package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class MultiVariantWalkerIntegrationTest extends CommandLineProgramTest {

    private static final String ENGINE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    @Override
    public String getTestedClassName() {
        return MultiVariantWalker.class.getSimpleName();
    }

    public static File getTestDataDir() {
        return new File(publicTestDir + "org/broadinstitute/hellbender/engine/MultiVariantDataSource/");
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithFeatures",
            oneLineSummary = "TestGATKToolWithFeatures",
            programGroup = TestProgramGroup.class
    )
    private static final class TestMultiVariantWalker extends MultiVariantWalker {

        int count = 0;

        SimpleInterval locus;

        public void apply(
                final VariantContext variant,
                final ReadsContext readsContext,
                final ReferenceContext referenceContext,
                final FeatureContext featureContext )
        {
            count++;

            //make sure we only move forward
            if (locus != null) {
                int locDiff = IntervalUtils.compareLocatables(
                        locus,
                        new SimpleInterval(variant),
                        this.getHeaderForVariants().getSequenceDictionary());
                Assert.assertTrue(locDiff == 0 || locDiff == -1);
            }
            locus = new SimpleInterval(variant);
        }
    }

    @Test(expectedExceptions = UserException.class)
    public void testQueryOverUnindexedFile() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(new File(ENGINE_TEST_DIRECTORY + "unindexed.vcf"), "unindexed");
        multiVariantSource.query(new SimpleInterval("1", 1, 1));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDuplicateSources() throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        File testFile = new File(getTestDataDir(), "baseVariants.vcf");
        args.add("--variant");
        args.add(testFile.getAbsolutePath());
        args.add("--variant");
        args.add(testFile.getAbsolutePath()); // add it again
        tool.instanceMain(args.toArray(new String[args.size()]));
    }

    @DataProvider(name="variantFiles")
    public Object[][] getVariantFiles() {
        return new Object[][]
            {
                { Arrays.asList(new File(getTestDataDir(), "baseVariants.vcf")), null, 26 },
                { Arrays.asList(new File(getTestDataDir(), "interleavedVariants_1.vcf")), null, 13 },
                { Arrays.asList(
                        new File(getTestDataDir(), "interleavedVariants_1.vcf"),
                        new File(getTestDataDir(), "interleavedVariants_2.vcf")), null, 26 },
                { Arrays.asList(
                        new File(getTestDataDir(), "splitVariants_1.vcf"),
                        new File(getTestDataDir(), "splitVariants_2.vcf")), null, 26 },

                // with intervals
                { Arrays.asList(new File(getTestDataDir(), "baseVariants.vcf")), "1", 14 },
                { Arrays.asList(new File(getTestDataDir(), "interleavedVariants_1.vcf")), "1", 7 },
                { Arrays.asList(
                        new File(getTestDataDir(), "interleavedVariants_1.vcf"),
                        new File(getTestDataDir(), "interleavedVariants_2.vcf")), "1", 14 },
                { Arrays.asList(
                        new File(getTestDataDir(), "interleavedVariants_1.vcf"),
                        new File(getTestDataDir(), "interleavedVariants_2.vcf")), "2:200-600", 3 },
            };
    }

    @Test(dataProvider = "variantFiles")
    public void testVariantOrder(final List<File> inputFiles, final String interval, final int expectedCount) throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>(inputFiles.size() * 2);
        inputFiles.forEach(f -> {args.add("--variant");
        args.add(f.getAbsolutePath());});
        if (interval != null) {
            args.add("-L");
            args.add(interval);
        }

        tool.instanceMain(args.toArray(new String[args.size()]));
        Assert.assertEquals(tool.count, expectedCount);
    }

    @Test
    public void testGetCompatibleHeader() throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        File testFile1 = new File(getTestDataDir(), "interleavedVariants_1.vcf");
        args.add("--variant");
        args.add(testFile1.getAbsolutePath());
        File testFile2 = new File(getTestDataDir(), "interleavedVariants_2.vcf");
        args.add("--variant");
        args.add(testFile2.getAbsolutePath());

        tool.instanceMain(args.toArray(new String[args.size()]));
        VCFHeader header = tool.getHeaderForVariants();
        Assert.assertEquals(header.getSequenceDictionary().getSequences().size(), 4 );
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testGetConflictingHeader() throws Exception {
        final TestMultiVariantWalker tool = new TestMultiVariantWalker();

        final List<String> args = new ArrayList<>();
        File testFile1 = new File(getTestDataDir(), "baseVariants.vcf");
        args.add("--variant");
        args.add(testFile1.getAbsolutePath());
        File testFile2 = new File(getTestDataDir(), "baseVariantsConflictingDictionary.vcf");
        args.add("--variant");
        args.add(testFile2.getAbsolutePath());

        tool.instanceMain(args.toArray(new String[args.size()]));
    }

}
