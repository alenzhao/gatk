package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public final class DepthPerAlleleBySampleUnitTest extends BaseTest {

    @Test
    public void testDescription(){
        Assert.assertEquals(new DepthPerAlleleBySample().getKeyNames(), Collections.singletonList(VCFConstants.GENOTYPE_ALLELE_DEPTHS));
        Assert.assertEquals(new DepthPerAlleleBySample().getDescriptions(), Collections.singletonList(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
    }

    @Test
    public void testUsingReads(){

        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AC = Arrays.asList(A, C);
        final int readDepthRef = 20;
        final int readDepthAlt = 17;
        final int[] extectedAD = {readDepthRef, readDepthAlt};

        final String sample1 = "sample1";
        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(sample1, AC).DP(dpDepth).make();

        final double log10PError = -5;

        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < readDepthAlt; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "readDepthAlt_" + i);
            read.setMappingQuality(20);
            reads.add(read);
        }
        for (int i = 0; i < readDepthRef; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "readDepthRef_" + i);
            read.setMappingQuality(20);
            reads.add(read);
        }

        //throw in one non-informative read
        final GATKRead badRead = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "non-informative");
        badRead.setMappingQuality(20);
        reads.add(badRead);

        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample1, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(sample1));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(AC);
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);

        int n = 0;
        for (int i = 0; i < readDepthAlt; i++) {
            matrix.set(0, n, -10.0);
            matrix.set(1, n, -1.0);
            n++;
        }
        for (int i = 0; i < readDepthRef; i++) {
            matrix.set(0, n, -1.0);
            matrix.set(1, n, -100.0);
            n++;
        }

        matrix.set(0, n, -1.0);
        matrix.set(1, n, -1.1);

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        new DepthPerAlleleBySample().annotate(null, vc, gAC, gb, likelihoods);
        final int[] ad = gb.make().getAD();
        Assert.assertEquals(ad, extectedAD);

        //now test a no-op
        final GenotypeBuilder gb1 = new GenotypeBuilder(gAC);
        new DepthPerAlleleBySample().annotate(null, vc, null, gb1, likelihoods);  //null genotype
        Assert.assertFalse(gb1.make().hasAD());
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testBlowUp(){

        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AC = Arrays.asList(A, C);
        final int readDepthRef = 20;
        final int readDepthAlt = 17;

        final String sample1 = "sample1";
        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(sample1, AC).DP(dpDepth).make();

        final double log10PError = -5;

        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < readDepthAlt; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
            read.setMappingQuality(20);
            reads.add(read);
        }
        for (int i = 0; i < readDepthRef; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
            read.setMappingQuality(20);
            reads.add(read);
        }

        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample1, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(sample1));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(A));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);

        int n = 0;
        for (int i = 0; i < readDepthAlt; i++) {
            matrix.set(0, n, -10.0);
            n++;
        }
        for (int i = 0; i < readDepthRef; i++) {
            matrix.set(0, n, -1.0);
            n++;
        }

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        //this blows up because there's no C allele in the map
        new DepthPerAlleleBySample().annotate(null, vc, gAC, gb, likelihoods);
    }


}
