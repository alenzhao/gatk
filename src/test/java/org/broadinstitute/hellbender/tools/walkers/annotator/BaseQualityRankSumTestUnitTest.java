package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_BaseQualityRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RankSumTest;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public final class BaseQualityRankSumTestUnitTest {

    private final String sample1 = "NA1";
    private final String sample2 = "NA2";

    private VariantContext makeVC(final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30, 0, 190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    private GATKRead makeRead(final byte qual) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
        read.setMappingQuality(50);
        read.setBaseQualities(Utils.dupBytes(qual, 10));
        return read;
    }

    @Test
    public void testBaseQual() {
        final InfoFieldAnnotation ann = new BaseQualityRankSumTest();
        final String key = GATKVCFConstants.BASE_QUAL_RANK_SUM_KEY;

        final Allele alleleRef = Allele.create("T", true);
        final Allele alleleAlt = Allele.create("A", false);

        final byte[] hardAlts = {10, 20};
        final byte[] hardRefs = {50, 60};
        final GATKRead read0 = makeRead(hardAlts[0]);
        final GATKRead read1 = makeRead(hardAlts[1]);
        final GATKRead read2 = makeRead(hardRefs[0]);
        final GATKRead read3 = makeRead(hardRefs[1]);

        final List<GATKRead> reads = Arrays.asList(read0, read1, read2, read3);
        final ReadLikelihoods<Allele> likelihoods =
                AnnotationArtificialData.initializeReadLikelihoods(sample1, alleleRef, alleleAlt, reads);

        AnnotationArtificialData.setLikelihoods(likelihoods, 2, 2, 0, -100.0, -100.0, -1.0);

        final ReferenceContext ref = null;
        final VariantContext vc = makeVC(alleleRef, alleleAlt);

        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);

        final double val = MannWhitneyU.runOneSidedTest(false,
                Arrays.asList(hardAlts[0], hardAlts[1]),
                Arrays.asList(hardRefs[0], hardRefs[1])).getLeft();
        final String valStr = String.format("%.3f", val);
        Assert.assertEquals(annotate.get(key), valStr);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key);
    }

    @Test
    public void testBaseQualRawAnnotate() {
        final AS_RankSumTest ann =  new AS_BaseQualityRankSumTest();
        final String key1 = GATKVCFConstants.AS_RAW_BASE_QUAL_RANK_SUM_KEY;
        final String key2 = GATKVCFConstants.AS_BASE_QUAL_RANK_SUM_KEY;

        final Allele alleleRef = Allele.create("T", true);
        final Allele alleleAlt = Allele.create("A", false);

        final byte[] hardAlts = {10, 20};
        final byte[] hardRefs = {50, 60};
        final GATKRead read0 = makeRead(hardAlts[0]);
        final GATKRead read1 = makeRead(hardAlts[1]);
        final GATKRead read2 = makeRead(hardRefs[0]);
        final GATKRead read3 = makeRead(hardRefs[1]);

        final List<GATKRead> reads = Arrays.asList(read0, read1, read2, read3);
        final ReadLikelihoods<Allele> likelihoods =
                AnnotationArtificialData.initializeReadLikelihoods(sample1, alleleRef, alleleAlt, reads);

        AnnotationArtificialData.setLikelihoods(likelihoods, 2, 2, 0, -100.0, -100.0, -1.0);

        final ReferenceContext ref = null;
        final VariantContext vc = makeVC(alleleRef, alleleAlt);

        final Map<String, Object> annotateRaw = ann.annotateRawData(ref, vc, likelihoods);
        final Map<String, Object> annotateNonRaw = ann.annotate(ref, vc, likelihoods);

        final String expectedAnnotation = hardRefs[0] + ",1," + hardRefs[1] + ",1" + AS_RankSumTest.PRINT_DELIM + hardAlts[0] + ",1," + hardAlts[1] + ",1";
        Assert.assertEquals(annotateRaw.get(key1),    expectedAnnotation);
        Assert.assertEquals(annotateNonRaw.get(key1), expectedAnnotation);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key1);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key2);
    }

    @Test
    public void testEmptyIfNoGenotypes() throws Exception {
        final BaseQualityRankSumTest ann = new BaseQualityRankSumTest();

        final List<GATKRead> reads = Collections.emptyList();
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample1, reads);
        final SampleList sampleList = new IndexedSampleList(Arrays.asList(sample1));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(Allele.NO_CALL));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        final Map<String, Object> annotate = ann.annotate(null, when(mock(VariantContext.class).getGenotypesOrderedByName()).thenReturn(Collections.<Genotype>emptyList()).getMock(), likelihoods);
        Assert.assertTrue(annotate.isEmpty());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMapNotNull() {
        final BaseQualityRankSumTest ann = new BaseQualityRankSumTest();
        ann.annotate(null, mock(VariantContext.class), null);
    }
}