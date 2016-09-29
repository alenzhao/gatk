package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;

public final class MultiVariantDataSourceUnitTest extends BaseTest {
    private static final String ENGINE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String MULTI_VARIANT_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/MultiVariantDataSource/";

    private static final File baseVariants = new File(MULTI_VARIANT_TEST_DIRECTORY, "baseVariants.vcf");
    private static final File baseVariantsAlternateDictionary = new File(MULTI_VARIANT_TEST_DIRECTORY, "baseVariantsAlternateDictionary.vcf");
    private static final File baseVariantsConflictingDictionary = new File(MULTI_VARIANT_TEST_DIRECTORY, "baseVariantsConflictingDictionary.vcf");

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testRequireFeatureInput() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.iterator();
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testRejectNonExistentFile() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(BaseTest.getSafeNonExistentFile("nonexistent.vcf"), "nonexistent");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testRejectNullFile() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(null, "sourceName1");
    }

    @Test(expectedExceptions = UserException.class)
    public void testQueryOverUnindexedFile() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(new File(ENGINE_TEST_DIRECTORY + "unindexed.vcf"), "unindexed");
        multiVariantSource.query(new SimpleInterval("1", 1, 1));
    }

    @Test
    public void testGetName() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();

        multiVariantSource.addFeatureDataSource(baseVariants, "sourceName1");
        String name = multiVariantSource.getName();
        Assert.assertTrue(name.contains("sourceName1"));

        multiVariantSource.addFeatureDataSource(baseVariantsAlternateDictionary, "sourceName2");
        name = multiVariantSource.getName();
        Assert.assertTrue(name.contains("sourceName1"));
        Assert.assertTrue(name.contains("sourceName2"));
    }

    @Test
    public void testGetSequenceDictionaryCompatible() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_1.vcf"), "interleavedVariants_1");
        multiVariantSource.addFeatureDataSource(new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_2.vcf"), "interleavedVariants_2");

        Assert.assertEquals(multiVariantSource.getSequenceDictionary().getSequences().size(), 4);
    }

    @Test
    public void testGetSequenceDictionaryAlternate() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(baseVariants, "baseVariants");
        multiVariantSource.addFeatureDataSource(baseVariantsAlternateDictionary, "baseVariantsAlternateDictionary");

        multiVariantSource.getSequenceDictionary();
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testGetSequenceDictionaryIncompatible() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(baseVariants, "baseVariants");
        multiVariantSource.addFeatureDataSource(baseVariantsConflictingDictionary, "baseVariantsConflictingDictionary");

        multiVariantSource.getSequenceDictionary();
    }

    @Test
    public void testIterator() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_1.vcf"), "interleavedVariants_1.vcf");
        multiVariantSource.addFeatureDataSource(new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_2.vcf"), "interleavedVariants_2.vcf");

        int count = 0;
        for (final VariantContext vc: multiVariantSource) {
            count++;
        };
        Assert.assertEquals(count, 26);
    }

    @Test
    public void testSerialQueries() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_1.vcf"), "interleavedVariants_1.vcf");
        multiVariantSource.addFeatureDataSource(new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_2.vcf"), "interleavedVariants_2.vcf");

        int count = 0;
        Iterator<VariantContext> it = multiVariantSource.query(new SimpleInterval("1", 1, 1200));
        while (it.hasNext()) {
            it.next();
            count++;
        };
        Assert.assertEquals(count, 14);

        count = 0;
        it = multiVariantSource.query(new SimpleInterval("2", 200, 600));
        while (it.hasNext()) {
            it.next();
            count++;
        };
        Assert.assertEquals(count, 3);
    }

    @Test
    public void testSetIntervals() {
        MultiVariantDataSource multiVariantSource = new MultiVariantDataSource();
        multiVariantSource.addFeatureDataSource(new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_1.vcf"), "interleavedVariants_1.vcf");
        multiVariantSource.addFeatureDataSource(new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_2.vcf"), "interleavedVariants_2.vcf");

        int count = 0;
        multiVariantSource.setIntervalsForTraversal(
                Arrays.asList(
                    new SimpleInterval("1", 1, 1200),
                    new SimpleInterval("2", 200, 600)
                )
        );
        for (final VariantContext vc: multiVariantSource) {
            count++;
        };
        Assert.assertEquals(count, 17);
    }

}

