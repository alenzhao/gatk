package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Created by davidben on 10/16/16.
 */
public class AnnotationArtificialData {
    private static final double MATCH_LIKELIHOOD = -1.0;

    public static ReadLikelihoods<Allele> makeLikelihoods(final String sample,
                                                          final int refDepth,
                                                          final int altDepth,
                                                          final int noninformativeDepth,
                                                          final double refReadAltLikelihood,
                                                          final double altReadRefLikelihood,
                                                          final double badReadAltLikelihood,
                                                          final Allele refAllele,
                                                          final Allele altAllele,
                                                          final byte baseQual,
                                                          final int mappingQual) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < altDepth; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "alt_" + i);
            read.setMappingQuality(mappingQual);
            read.setBaseQualities(Utils.dupBytes(baseQual, 10));
            reads.add(read);
        }

        for (int i = 0; i < refDepth; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "ref_" + i);
            read.setBaseQualities(Utils.dupBytes(baseQual, 10));
            read.setMappingQuality(mappingQual);
            reads.add(read);
        }

        for (int i = 0; i < noninformativeDepth; i++) {
            final GATKRead badRead = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "bad_" + i);
            badRead.setBaseQualities(Utils.dupBytes(baseQual, 10));
            badRead.setMappingQuality(mappingQual);
            reads.add(badRead);
        }

        final ReadLikelihoods<Allele> likelihoods = initializeReadLikelihoods(sample, refAllele, altAllele, reads);

        setLikelihoods(likelihoods, refDepth, altDepth, noninformativeDepth, refReadAltLikelihood, altReadRefLikelihood, badReadAltLikelihood);

        return likelihoods;
    }

    public static ReadLikelihoods<Allele> initializeReadLikelihoods(String sample, Allele refAllele, Allele altAllele, List<GATKRead> reads) {
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(sample));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(refAllele, altAllele));
        return new ReadLikelihoods<>(sampleList, alleleList, readsBySample);
    }

    public static void setLikelihoods(final ReadLikelihoods<Allele> likelihoods,
                                      final int refDepth,
                                      final int altDepth,
                                      final int noninformativeDepth,
                                      final double refReadAltLikelihood,
                                      final double altReadRefLikelihood,
                                      final double badReadAltLikelihood) {
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);

        int readIndex = 0;
        for (int i = 0; i < altDepth; i++) {
            matrix.set(0, readIndex, altReadRefLikelihood);
            matrix.set(1, readIndex, MATCH_LIKELIHOOD);
            readIndex++;
        }

        for (int i = 0; i < refDepth; i++) {
            matrix.set(0, readIndex, MATCH_LIKELIHOOD);
            matrix.set(1, readIndex, refReadAltLikelihood);
            readIndex++;
        }

        for (int i = 0; i < noninformativeDepth; i++) {
            matrix.set(0, readIndex, MATCH_LIKELIHOOD);
            matrix.set(1, readIndex, badReadAltLikelihood);
            readIndex++;
        }
    }
}
