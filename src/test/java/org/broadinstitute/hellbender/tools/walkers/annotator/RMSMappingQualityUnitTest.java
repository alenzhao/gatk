package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Sets;
import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public final class RMSMappingQualityUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull() throws Exception {
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new RMSMappingQuality();
        cov.annotate(referenceContext, vc, null); //vc can't be null
    }

    @Test
    public void testDescriptions() throws Exception {
        final InfoFieldAnnotation cov = new RMSMappingQuality();
        Assert.assertEquals(cov.getDescriptions().size(), 2);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.RMS_MAPPING_QUALITY_KEY);
        Assert.assertEquals(cov.getDescriptions().get(1).getID(), GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY);
        Assert.assertEquals(new RMSMappingQuality().getRawKeyName(), GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY);
        Assert.assertEquals(new RMSMappingQuality().getKeyNames(), Sets.newHashSet(VCFConstants.RMS_MAPPING_QUALITY_KEY, GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY));
    }

    @Test
    public void testNullLikelihoods() throws Exception {
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new RMSMappingQuality();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null);
        Assert.assertTrue(annotate.isEmpty());

        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.RMS_MAPPING_QUALITY_KEY);
        Assert.assertEquals(cov.getDescriptions().get(1).getID(), GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.RMS_MAPPING_QUALITY_KEY);
    }

    private VariantContext makeVC() {
        final GenotypesContext testGC = GenotypesContext.create(2);
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    @Test
    public void testLikelihoods(){

        final Allele alleleA = Allele.create("A");
        final double lik= -1.0;  //ignored

        final int[] MQs = {1,2,3,4,5,6,7,8,9,10, QualityUtils.MAPPING_QUALITY_UNAVAILABLE};
        final List<Integer> MQsList = Arrays.asList(ArrayUtils.toObject(MQs));

        //MQ 255 are excluded from the calculations, we test it here.
        final List<Integer> MQsListOK = new ArrayList<>(MQsList);
        //NOTE: if we just call remove(i), Java thinks i is an index.
        //A workaround for this overloading bogosity to to call removeAll and pass a collection
        //(casting i to (Object) would work too but it's more error prone)
        MQsListOK.removeAll(Collections.singleton(QualityUtils.MAPPING_QUALITY_UNAVAILABLE));

        final int n1A= MQs.length;
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < n1A; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
            read.setMappingQuality(MQs[i]);
            reads.add(read);
        }

        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(alleleA));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);
        for (int n = 0; n < n1A; n++) {
            matrix.set(0, n, lik);
        }

        final VariantContext vc = makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new RMSMappingQuality().annotate(referenceContext, vc, likelihoods);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY), "annots");
        final double rms= MathUtils.sumOfSquares(MQsListOK); //only those are MQ0
        Assert.assertNull(annotate.get(VCFConstants.RMS_MAPPING_QUALITY_KEY));
        Assert.assertEquals(annotate.get(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY), String.format("%.2f", rms));
    }

    @Test
    public void testEmptyLikelihoods() throws Exception {
        final List<GATKRead> reads = Collections.emptyList();
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(Allele.create("A")));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);

        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new RMSMappingQuality().annotate(referenceContext, vc, likelihoods);
        Assert.assertTrue(annotate.isEmpty());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull_AS() throws Exception {
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new AS_RMSMappingQuality();
        cov.annotate(referenceContext, vc, null); //vc can't be null
    }

    @Test
    public void testDescriptions_AS() throws Exception {
        final InfoFieldAnnotation cov = new AS_RMSMappingQuality();
        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY);
    }

    @Test
    public void testNullLikelihoods_AS() throws Exception {
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new AS_RMSMappingQuality();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null);
        Assert.assertTrue(annotate.isEmpty());

        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY);
    }

    @Test
    public void testLikelihoods_AS(){

        final Allele alleleA = Allele.create("A", true);
        final double lik= -1.0;  //ignored

        final int[] MQs = {1,2,3,4,5,6,7,8,9,10, QualityUtils.MAPPING_QUALITY_UNAVAILABLE};
        final List<Integer> MQsList = Arrays.asList(ArrayUtils.toObject(MQs));

        //MQ 255 are excluded from the calculations, we test it here.
        final List<Integer> MQsListOK = new ArrayList<>(MQsList);
        //NOTE: if we just call remove(i), Java thinks i is an index.
        //A workaround for this overloading bogosity to to call removeAll and pass a collection
        //(casting i to (Object) would work too but it's more error prone)
        MQsListOK.removeAll(Collections.singleton(QualityUtils.MAPPING_QUALITY_UNAVAILABLE));

        final List<GATKRead> reads = new ArrayList<>();
        for (int MQ : MQs) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
            read.setMappingQuality(MQ);
            reads.add(read);
        }

        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(alleleA));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);
        for (int n = 0; n < reads.size(); n++) {
            matrix.set(0, n, lik);
        }

        final VariantContext vc = makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new AS_RMSMappingQuality().annotateRawData(referenceContext, vc, likelihoods);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY), "annots");
        final double rms= MathUtils.sumOfSquares(MQsListOK); //only those are MQ0
        final String[] split =((String)annotate.get(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY)).split(AS_RMSMappingQuality.SPLIT_DELIM);
        Assert.assertEquals(split.length, 2);
        Assert.assertEquals(split[0], String.format("%.2f", rms));
        Assert.assertEquals(split[1], String.format("%.2f", 0.0));
    }

    @Test
    public void testLikelihoodsEmpty_AS() throws Exception {
        final List<GATKRead> reads = Collections.emptyList();
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(Allele.create("A")));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new AS_RMSMappingQuality().annotate(referenceContext, vc, likelihoods);
        final String[] split =((String)annotate.get(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY)).split(AS_RMSMappingQuality.SPLIT_DELIM);
        Assert.assertEquals(split.length, 2);
        Assert.assertEquals(split[0], String.format("%.2f", 0.0));
        Assert.assertEquals(split[1], String.format("%.2f", 0.0));
    }
}
