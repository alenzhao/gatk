package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_MappingQualityRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RankSumTest;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public final class MappingQualityRankSumTestUnitTest {

    private final String sample1 = "NA1";
    private final String sample2 = "NA2";

    private VariantContext makeVC( final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    private GATKRead makeRead(final int mq) {
        Cigar cigar = TextCigarCodec.decode("10M");
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(mq);
        return read;
    }

    @Test
    public void testMQ(){
        final InfoFieldAnnotation ann = new MappingQualityRankSumTest();
        final String key = GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY;

        final Allele alleleRef = Allele.create("T", true);
        final Allele alleleAlt = Allele.create("A", false);

        final int[] hardAlts = {10, 20};
        final int[] hardRefs = {100, 110};
        final GATKRead read0 = makeRead(hardAlts[0]);
        final GATKRead read1 = makeRead(hardAlts[1]);
        final GATKRead read2 = makeRead(hardRefs[0]);
        final GATKRead read3 = makeRead(hardRefs[1]);

        final List<GATKRead> reads = Arrays.asList(read0, read1, read2, read3);
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample1, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(sample1));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(alleleRef, alleleAlt));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);
        matrix.set(1, 0, -1.0);
        matrix.set(0, 0, -100.0);
        matrix.set(1, 1, -1.0);
        matrix.set(0, 1, -100.0);
        matrix.set(1, 2, -100.0);
        matrix.set(0, 2, -1.0);
        matrix.set(1, 3, -100.0);
        matrix.set(0, 3, -1.0);



        final ReferenceContext ref= null;
        final VariantContext vc= makeVC(alleleRef, alleleAlt);

        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);

        final double val= MannWhitneyU.runOneSidedTest(false, Arrays.asList(hardAlts[0], hardAlts[1]),
                                                              Arrays.asList(hardRefs[0], hardRefs[1])).getLeft();
        final String valStr= String.format("%.3f", val);
        Assert.assertEquals(annotate.get(key), valStr);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key);
    }

    @Test
    public void testAS_MQRaw(){
        final AS_RankSumTest ann = new AS_MappingQualityRankSumTest();
        final String key1 = GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY;
        final String key2 = GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY;

        final Allele alleleRef = Allele.create("T", true);
        final Allele alleleAlt = Allele.create("A", false);

        final int[] hardAlts = {10, 20};
        final int[] hardRefs = {100, 110};
        final GATKRead read0 = makeRead(hardAlts[0]);
        final GATKRead read1 = makeRead(hardAlts[1]);
        final GATKRead read2 = makeRead(hardRefs[0]);
        final GATKRead read3 = makeRead(hardRefs[1]);

        final List<GATKRead> reads = Arrays.asList(read0, read1, read2, read3);
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample1, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(sample1));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(alleleRef, alleleAlt));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);
        matrix.set(1, 0, -1.0);
        matrix.set(0, 0, -100.0);
        matrix.set(1, 1, -1.0);
        matrix.set(0, 1, -100.0);
        matrix.set(1, 2, -100.0);
        matrix.set(0, 2, -1.0);
        matrix.set(1, 3, -100.0);
        matrix.set(0, 3, -1.0);


        final ReferenceContext ref= null;
        final VariantContext vc= makeVC(alleleRef, alleleAlt);

        final Map<String, Object> annotateRaw = ann.annotateRawData(ref, vc, likelihoods);
        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);

        final String expected = hardRefs[0] + ",1," + hardRefs[1] + ",1" + AS_RankSumTest.PRINT_DELIM + hardAlts[0] + ",1," + hardAlts[1] + ",1";
        Assert.assertEquals(annotate.get(key1), expected);
        Assert.assertEquals(annotateRaw.get(key1), expected);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key1);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key2);
    }
}
