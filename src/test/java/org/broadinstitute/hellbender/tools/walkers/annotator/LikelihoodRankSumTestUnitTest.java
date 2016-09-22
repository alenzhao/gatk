package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public final class LikelihoodRankSumTestUnitTest extends BaseTest {
    private final String sample1 = "NA1";
    private final String sample2 = "NA2";

    private VariantContext makeVC( final String contig, final long position, final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr(contig).start(position).stop(position).genotypes(testGC).make();
    }

    private GATKRead makeRead(final String contig, final int start, final int mq, final String name) {
        Cigar cigar = TextCigarCodec.decode("10M");
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar, name);
        read.setMappingQuality(mq);
        read.setPosition(contig, start);
        return read;
    }
    @Test
    public void testReadPos(){
        final String contig = "1";

        final Allele alleleRef = Allele.create("T", true);
        final Allele alleleAlt = Allele.create("A", false);

        final double[] altBestAlleleLL = {-1.0, -2.0};
        final double[] refBadAlleleLL =  {-100.0, -100.0};

        final double[] altBadAlleleLL =  {-100.0, -100.0};
        final double[] refBestAlleleLL = {-5.0, -7.0};
        final GATKRead read0 = makeRead(contig, 1,  30, "read1");
        final GATKRead read1 = makeRead(contig, 1, 30, "read2");
        final GATKRead read2 = makeRead(contig, 1, 30, "read3");
        final GATKRead read3 = makeRead(contig, 1, 30, "read4");

        final List<GATKRead> reads = Arrays.asList(read0, read1, read2, read3);
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample1, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(sample1, sample2));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(alleleRef, alleleAlt));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);
        matrix.set(1, 0, altBestAlleleLL[0]);
        matrix.set(0, 0, refBadAlleleLL[0]);
        matrix.set(1, 1, altBestAlleleLL[1]);
        matrix.set(0, 1, refBadAlleleLL[0]);
        matrix.set(1, 2, altBadAlleleLL[0]);
        matrix.set(0, 2, refBestAlleleLL[0]);
        matrix.set(1, 3, altBadAlleleLL[1]);
        matrix.set(0, 3, refBestAlleleLL[1]);


        final InfoFieldAnnotation ann = new LikelihoodRankSumTest();
        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), GATKVCFConstants.LIKELIHOOD_RANK_SUM_KEY);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), GATKVCFConstants.LIKELIHOOD_RANK_SUM_KEY);

        final ReferenceContext ref= null;

        final long position = 5L;  //middle of the read
        final VariantContext vc= makeVC(contig, position, alleleRef, alleleAlt);

        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);
        final double val= MannWhitneyU.runOneSidedTest(false,
                Arrays.asList(altBestAlleleLL[0], altBestAlleleLL[1]),
                Arrays.asList(refBestAlleleLL[0], refBestAlleleLL[1])).getLeft();
        final String valStr= String.format("%.3f", val);
        Assert.assertEquals(annotate.get(GATKVCFConstants.LIKELIHOOD_RANK_SUM_KEY), valStr);

    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testBlowupOnCallWithoutAllele(){
        final LikelihoodRankSumTest ann = new LikelihoodRankSumTest();
        final String contig = "1";
        final GATKRead read1 = makeRead(contig, 1,  30, "read1");
        ann.getElementForRead(read1, read1.getStart());
    }

}
