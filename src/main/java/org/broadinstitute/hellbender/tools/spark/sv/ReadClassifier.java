package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Function;

/**
 * Figures out what kind of BreakpointEvidence, if any, a read represents.
 */
public class ReadClassifier implements Function<GATKRead, Iterator<BreakpointEvidence>> {
    @VisibleForTesting static final int MIN_SOFT_CLIP_LEN = 30; // minimum length of an interesting soft clip
    @VisibleForTesting static final int MIN_INDEL_LEN = 40; // minimum length of an interesting indel
    private static final byte MIN_QUALITY = 15; // minimum acceptable quality in a soft-clip window
    private static final float MAX_ZISH_SCORE = 6.f; // maximum fragment-length "z" score for a normal fragment
    private static final float MIN_CRAZY_ZISH_SCORE = 100.f; // "z" score that's probably associated with a mapping error
    private final ReadMetadata readMetadata;

    public ReadClassifier( final ReadMetadata readMetadata ) {
        this.readMetadata = readMetadata;
    }

    @Override
    public Iterator<BreakpointEvidence> apply( final GATKRead read ) {
        if ( read.isUnmapped() ) return Collections.emptyIterator();

        final List<BreakpointEvidence> evidenceList = new ArrayList<>();

        checkForSplitRead(read, evidenceList);
        checkDiscordantPair(read, evidenceList);

        return evidenceList.iterator();
    }

    @VisibleForTesting void checkForSplitRead( final GATKRead read,
                                               final List<BreakpointEvidence> evidenceList ) {
        final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
        if ( hasInitialSoftClip(cigarElements, read) ) {
            evidenceList.add(new BreakpointEvidence.SplitRead(read, readMetadata, true));
        }
        if ( hasFinalSoftClip(cigarElements, read) ) {
            evidenceList.add(new BreakpointEvidence.SplitRead(read, readMetadata, false));
        }
        checkBigIndel(cigarElements, read, evidenceList);
    }

    private static boolean hasInitialSoftClip( final List<CigarElement> cigarElements, final GATKRead read ) {
        final ListIterator<CigarElement> itr = cigarElements.listIterator();
        if ( !itr.hasNext() ) return false;

        CigarElement firstEle = itr.next();
        if ( firstEle.getOperator() == CigarOperator.HARD_CLIP && itr.hasNext() ) {
            firstEle = itr.next();
        }
        return firstEle.getOperator() == CigarOperator.SOFT_CLIP &&
                firstEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                isHighQualityRegion(read.getBaseQualities(), 0);
    }

    private static boolean hasFinalSoftClip( final List<CigarElement> cigarElements, final GATKRead read ) {
        final ListIterator<CigarElement> itr = cigarElements.listIterator(cigarElements.size());
        if ( !itr.hasPrevious() ) return false;

        CigarElement lastEle = itr.previous();
        if ( lastEle.getOperator() == CigarOperator.HARD_CLIP && itr.hasPrevious() ) {
            lastEle = itr.previous();
        }
        return lastEle.getOperator() == CigarOperator.SOFT_CLIP &&
                lastEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                isHighQualityRegion(read.getBaseQualities(), read.getLength() - lastEle.getLength());
    }

    private static boolean isHighQualityRegion( final byte[] quals, int idx ) {
        for ( final int end = idx+MIN_SOFT_CLIP_LEN; idx != end; ++idx ) {
            if ( quals[idx] < MIN_QUALITY ) return false;
        }
        return true;
    }

    private void checkBigIndel( final List<CigarElement> cigarElements,
                               final GATKRead read,
                               final List<BreakpointEvidence> evidenceList ) {
        int locus = read.getStart();
        for ( final CigarElement ele : cigarElements ) {
            final CigarOperator op = ele.getOperator();
            if ( ele.getLength() >= MIN_INDEL_LEN ) {
                if ( op == CigarOperator.INSERTION ) {
                    evidenceList.add(new BreakpointEvidence.LargeIndel(read, readMetadata, locus));
                } else if ( op == CigarOperator.DELETION ) {
                    evidenceList.add(new BreakpointEvidence.LargeIndel(read, readMetadata, locus+ele.getLength()/2));
                }
            }
            if ( op.consumesReferenceBases() ) {
                locus += ele.getLength();
            }
        }
    }

    @VisibleForTesting void checkDiscordantPair( final GATKRead read, final List<BreakpointEvidence> evidenceList ) {
        if ( read.mateIsUnmapped() ) {
            evidenceList.add(new BreakpointEvidence.MateUnmapped(read, readMetadata));
        } else if ( !read.getContig().equals(read.getMateContig()) ) {
            evidenceList.add(new BreakpointEvidence.InterContigPair(read, readMetadata));
        } else if ( read.isReverseStrand() == read.mateIsReverseStrand() ) {
            evidenceList.add(new BreakpointEvidence.SameStrandPair(read, readMetadata));
        } else if ( (read.getStart() < read.getMateStart() && read.isReverseStrand()) ||
                    (read.getStart() > read.getMateStart() && !read.isReverseStrand()) ) {
            evidenceList.add(new BreakpointEvidence.OutiesPair(read, readMetadata));
        } else {
            final float zIshScore =
                    readMetadata.getStatistics(read.getReadGroup()).getZIshScore(Math.abs(read.getFragmentLength()));
            if ( zIshScore > MAX_ZISH_SCORE && zIshScore < MIN_CRAZY_ZISH_SCORE ) {
                evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(read, readMetadata));
            }
            // TODO: see if there's anything we can do about anomalously short fragment sizes
            // (With current fragment sizes and read lengths there aren't enough bases to have a >=50bp insertion
            // between the mates.)
        }
    }
}
