package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.util.ArrayList;
import java.util.List;
import java.util.Spliterator;

/**
 * A MultiVariantWalker is a tool that processes one variant at a time from multiple sources of variants, with
 * optional contextual information from a reference, sets of reads, and/or supplementary sources of Features.
 *
 * VariantWalker authors must implement the {@link #apply} method to process each variant, and may optionally implement
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class MultiVariantWalker extends VariantWalkerBase {

    // NOTE: using File rather than FeatureInput<VariantContext> here so that we can keep this driving source
    //       of variants separate from any other potential sources of Features
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
    public List<String> drivingVariantFiles = new ArrayList<String>();

    // NOTE: keeping the driving source of variants separate from other, supplementary FeatureInputs in our FeatureManager
    // in GATKTool we do add the driving source to the Feature manager but we do need to treat it differently and thus this
    // field.
    private MultiVariantDataSource drivingVariants = new MultiVariantDataSource();
    private List<FeatureInput<VariantContext>> drivingVariantsFeatureInputs = new ArrayList<>(2);

    @Override
    public boolean requiresFeatures() { return true; }

    protected SAMSequenceDictionary getSequenceDictionaryForDrivingVariants() { return drivingVariants.getSequenceDictionary(); }

    protected Spliterator<VariantContext> getSpliteratorForDrivingVariants() { return drivingVariants.spliterator(); }

    /**
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalStart} instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
        if ( hasIntervals() ) {
            drivingVariants.setIntervalsForTraversal(intervalsForTraversal);
        }
    }

    @SuppressWarnings("unchecked")
    protected void initializeDrivingVariants() {
        drivingVariantFiles.stream().forEach(
                f -> {
                    FeatureInput<VariantContext> featureInput = new FeatureInput<>(f, "drivingVariantFile");
                    if (drivingVariantsFeatureInputs.contains(featureInput)) {
                        //TODO: should we allow duplicates file names? If so we need to change the FeatureInput hash
                        //so FeatureManager handles them
                        throw new IllegalArgumentException("Feature inputs must be unique: " + featureInput.toString());
                    }
                    drivingVariants.addFeatureDataSource(featureInput, FEATURE_CACHE_LOOKAHEAD, VariantContext.class);
                    drivingVariantsFeatureInputs.add(featureInput);

                    //This is the data source for the driving source of variants, which uses a cache lookahead of FEATURE_CACHE_LOOKAHEAD
                    //drivingVariants = new FeatureDataSource<>(drivingVariantsFeatureInput, FEATURE_CACHE_LOOKAHEAD, VariantContext.class);

                    //Add the driving datasource to the feature manager too so that it can be queried. Setting lookahead to 0 to avoid caching.
                    //Note: we are disabling lookahead here because of windowed queries that need to "look behind" as well.
                    features.addToFeatureSources(0, featureInput, VariantContext.class);
                }
        );

        //Note: the intervals for the driving variants are set in onStartup
    }

    /**
     * Returns a list of feature inputs used for the driving variants for this source.
     */
    protected final List<FeatureInput<VariantContext>> getDrivingVariantsFeatureInputs() {
        return drivingVariantsFeatureInputs;
    }

    /**
     * Gets the header associated with our driving source of variants as a VCFHeader.
     *
     * @return VCFHeader for our driving source of variants
     */
    public final VCFHeader getHeaderForVariants() {
        return drivingVariants.getHeader();
    }

    /**
     * Close all data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalSuccess} and/or
     * {@link #closeTool} instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();
        drivingVariants.close();
    }
}
