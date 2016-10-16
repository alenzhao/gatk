package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public final class MappingQualityZeroUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull() throws Exception {
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new MappingQualityZero();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null);//vc can't be null
    }

    @Test
    public void testNullLikelihoods() throws Exception {
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new MappingQualityZero();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null);
        Assert.assertTrue(annotate.isEmpty());

        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.MAPPING_QUALITY_ZERO_KEY);
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

        final Allele alleleT = Allele.create("T");
        final Allele alleleA = Allele.create("A");
        final double lik= -1.0;

        final int n1A= 3;
        final int n1T= 5;
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < n1A; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "n1A_" + i);
            read.setMappingQuality(10);
            reads.add(read);
        }
        for (int i = 0; i < n1T; i++) {
            //try to fool it - add 2 alleles for same read
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "n1T_" + i);
            read.setMappingQuality(0);
            reads.add(read);
        }

        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(alleleT, alleleA));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);
        int n = 0;
        for (int i = 0; i < n1A; i++) {
            matrix.set(1, n, lik);
            n++;
        }
        for (int i = 0; i < n1T; i++) {
            matrix.set(1, n, lik);
            matrix.set(0, n, lik);

            n++;
        }

        final VariantContext vc = makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new MappingQualityZero().annotate(referenceContext, vc, likelihoods);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(VCFConstants.MAPPING_QUALITY_ZERO_KEY), "annots");
        Assert.assertEquals(annotate.get(VCFConstants.MAPPING_QUALITY_ZERO_KEY), String.valueOf(n1T));

    }

    @Test
    public void testEmptyLikelihoods() throws Exception {
        final List<GATKRead> reads = Collections.emptyList();
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(Allele.create("A")));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new MappingQualityZero().annotate(referenceContext, vc, likelihoods);

        final int n= 0; //strangely,  MappingQualityZero returns 0 if likelihoods is empty
        Assert.assertEquals(annotate.get(VCFConstants.MAPPING_QUALITY_ZERO_KEY), String.valueOf(n));
    }
}
