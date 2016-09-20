package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.tribble.*;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * MultiVariantDataSource aggregates multiple FeatureDataSources of variants, and enables traversals and queries
 * over those sources through a single interface.
 *
 * Two basic operations are available on this data source:
 *
 * -Iteration over all Variants in the data sources, optionally restricted to Variants overlapping
 *  a set of intervals if intervals are provided via {@link #setIntervalsForTraversal(List)}. Traversal
 *  by a set of intervals requires the files to have been indexed using the bundled tool IndexFeatureFile.
 *  The set of intervals provided MUST be non-overlapping and sorted in increasing order of start position.
 *
 * -Targeted queries by one interval at a time. This also requires the files to have been indexed using
 *  the bundled tool IndexFeatureFile. Targeted queries by one interval at a time are unaffected by
 *  any intervals for full traversal set via {@link #setIntervalsForTraversal(List)}.
 */
public final class MultiVariantDataSource implements GATKDataSource<VariantContext>, AutoCloseable {
    private static final Logger logger = LogManager.getLogger(MultiVariantDataSource.class);

    /**
     * List of FeatureDataSource objects aggregated by this MultiVariantDataSource
     */
    private List<FeatureDataSource<VariantContext>> featureDataSources = new ArrayList<>();

    /**
     * Merged VCF header used for this (aggregate) source, derived from the individual soruces.
     */
    private VCFHeader mergedHeader;

    /**
     * SAMSequenceDictionary used for this (aggregate) source, derived from the individual sources.
     */
    private SAMSequenceDictionary mergedDictionary;

    /**
     * Iterator representing an open traversal over this data source initiated via a call to {@link #iterator}
     * (null if there is no open traversal). We need this to ensure that each iterator is properly closed,
     * and to enforce the constraint (required by Tribble) that we never have more than one iterator open
     * over our feature reader.
     */
    private MergingIterator<VariantContext> currentIterator;

    /**
     * Default value for queryLookaheadBases, if none is specified. This is designed to be large enough
     * so that in typical usage (ie., query intervals with gradually increasing start locations) there will
     * be a substantial number of cache hits between cache misses, reducing the number of times we need to
     * repopulate the cache from disk.
     */
    public static final int DEFAULT_QUERY_LOOKAHEAD_BASES = 1000;

    /**
     * Define the set of incompatibilities that are allowed for SequenceDictionaries in Variant input data sources
     */
    private final static Set<SequenceDictionaryUtils.SequenceDictionaryCompatibility> allowedSequenceDictionaryIncompatibilities =
            new HashSet<>(Arrays.asList(
                    SequenceDictionaryUtils.SequenceDictionaryCompatibility.COMMON_SUBSET,
                    SequenceDictionaryUtils.SequenceDictionaryCompatibility.SUPERSET,
                    SequenceDictionaryUtils.SequenceDictionaryCompatibility.NO_COMMON_CONTIGS));

    /**
     * Creates a FeatureDataSource backed by the provided File and assigns this data source the specified logical
     * name. We will look ahead the default number of bases ({@link #DEFAULT_QUERY_LOOKAHEAD_BASES}) during queries
     * that produce cache misses.
     *
     * @param featureFile file containing Features
     * @param name logical name for this data source (may be null)
     */
    public void addFeatureDataSource(final File featureFile, final String name) {
        featureDataSources.add(new FeatureDataSource<>(featureFile, name, DEFAULT_QUERY_LOOKAHEAD_BASES));
        invalidateCachedHeaderAndDictionary();
    }

    /**
     * Creates a FeatureDataSource backed by the provided File and assigns this data source the specified logical
     * name. We will look ahead the specified number of bases during queries that produce cache misses.
     *
     * @param featureFile file containing Features
     * @param name logical name for this data source (may be null)
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     */
    public void addFeatureDataSource(final File featureFile, final String name, final int queryLookaheadBases){
        addFeatureDataSource(Utils.nonNull(featureFile).getAbsolutePath(), name, queryLookaheadBases, null);
    }

    /**
     * Creates a FeatureDataSource backed by the resource at the provided path.
     *
     * @param featurePath path to file or GenomicsDB url containing features
     * @param name logical name for this data source (may be null)
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                          that produce this type of Feature. May be null, which results in an unrestricted search.
     */
    public void addFeatureDataSource(final String featurePath, final String name, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType) {
        addFeatureDataSource(new FeatureInput<>(featurePath, name != null ? name : featurePath), queryLookaheadBases, targetFeatureType);
    }

    /**
     * Creates a FeatureDataSource backed by the provided FeatureInput. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureInput a FeatureInput specifying a source of Features
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     * @param targetFeatureType When searching for a {@link FeatureCodec} for this data source, restrict the search to codecs
     *                          that produce this type of Feature. May be null, which results in an unrestricted search.
     */
    public void addFeatureDataSource(final FeatureInput<VariantContext> featureInput, final int queryLookaheadBases, final Class<? extends Feature> targetFeatureType) {
        Utils.validateArg( queryLookaheadBases >= 0, "Query lookahead bases must be >= 0");

        featureDataSources.add(new FeatureDataSource<>(featureInput, queryLookaheadBases, targetFeatureType));
        invalidateCachedHeaderAndDictionary();
    }

    /**
     * Returns the aggregate sequence dictionary for this source of Variants. Uses the dictionary resulting
     * from merging the individual VCF headers (if present) for variant inputs, otherwise attempts to create a
     * sequence dictionary from an index file (if present).
     *
     * @return the sequence dictionary derived from the input sources, or null if no dictionary could be created
     * from either the header or an index file.
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        validateFeatureSources();
        if (mergedDictionary == null && null != getHeader()) {
            validateSequenceDictionaries();
            final VCFHeader header = getHeader();
            mergedDictionary = header.getSequenceDictionary();
            if (mergedDictionary == null || mergedDictionary.isEmpty()) {
                // last resort - try to get a sequence dictionary from one of the feature data sources (which
                // will create one from the index if necessary)
                for (FeatureDataSource<VariantContext> fds : featureDataSources) {
                    mergedDictionary = fds.getSequenceDictionary();
                    //TODO: we're taking the first one we find
                    //should we even try to do this, or just return null ?
                    if (mergedDictionary != null) {
                        break;
                    }
                }
            }
        }
        return mergedDictionary;
    }

    /**
     * Restricts traversals of this data source via {@link #iterator} to only return Features that overlap the provided
     //* intervals. Calls to {@link #query(SimpleInterval)} and/or {@link # queryAndPrefetch(SimpleInterval)} are not
     * affected by these intervals.
     *
     * Intervals MUST be non-overlapping and sorted in order of increasing start position, otherwise traversal
     * results will be incorrect.
     *
     * Passing in a null or empty interval List clears the intervals for traversal, making future iterations
     * over this data source unrestricted by intervals.
     *
     * @param intervals Our next full traversal will return only Features overlapping these intervals
     */
    public void setIntervalsForTraversal( final List<SimpleInterval> intervals ) {
        validateFeatureSources();
        featureDataSources.forEach(ds -> ds.setIntervalsForTraversal(intervals));
    }

    /**
     * Gets an iterator over all Features in this data source, restricting traversal to Features
     * overlapping our intervals if intervals were provided via {@link #setIntervalsForTraversal(List)}
     *
     * Calling this method invalidates (closes) any previous iterator obtained from this method.
     *
     * @return an iterator over all Features in this data source, limited to Features that overlap the intervals supplied via {@link #setIntervalsForTraversal(List)} (if intervals were provided)
     */
    @Override
    public Iterator<VariantContext> iterator() {
        return getMergedIteratorFromDataSources(ds -> ds.iterator());
    }

    /**
     * Gets an iterator over all Variants in this data source that overlap the provided interval.
     *
     * This operation is not affected by intervals provided via {@link #setIntervalsForTraversal(List)}.
     *
     * Requires the backing files to have been indexed using the IndexFeatureFile tool, and to
     * be sorted in increasing order of start position for each contig.
     *
     * Query results are cached to improve the performance of future queries during typical access
     * patterns. See notes to the class as a whole for a description of the caching strategy.
     *
     * Calling this method potentially invalidates (closes) any other open iterator obtained
     * from this data source via a call to {@link #iterator}
     *
     * @param interval retrieve all Variants overlapping this interval
     * @return an iterator over all Variants in this data source that overlap the provided interval
     */
    @Override
    public Iterator<VariantContext> query( final SimpleInterval interval ) {
        return getMergedIteratorFromDataSources(ds -> ds.queryAndPrefetch(interval).iterator());
    }

    /**
     * Close any existing iterator, create a new iterator and update the local cached iterator reference.
     * @param iteratorFromSource function to retrieve individual iterator, to be applied to each data source
     * @return
     */
    private Iterator<VariantContext> getMergedIteratorFromDataSources(
            final Function<FeatureDataSource<VariantContext>, Iterator<VariantContext>> iteratorFromSource) {
        validateFeatureSources();

        // Tribble documentation states that having multiple iterators open simultaneously over the same FeatureReader
        // results in undefined behavior
        closeOpenIterationIfNecessary();

        List<CloseableIterator<VariantContext>> iterators = new ArrayList<>(1);
        featureDataSources.forEach(ds -> iterators.add(getCloseableIteratorWrapper(iteratorFromSource.apply((ds)))));

        VariantContextComparator varComparator = new VariantContextComparator(getSequenceDictionary());
        currentIterator = new MergingIterator<>(varComparator, iterators);
        return currentIterator;
    }

    /**
     * Get the logical name of this data source.
     *
     * @return the logical name of this data source
     */
    public String getName() {
        validateFeatureSources();
        return "MultiVariantDataSource: ("
                + Utils.join(", ", featureDataSources.stream().map(fds -> fds.getName()).collect(Collectors.toList()))
                + ")";
    }

    /**
     * Gets the header associated with this data source
     *
     * @return header associated with this data source as an Object
     */
    public VCFHeader getHeader() {
        validateFeatureSources();
        if (mergedHeader == null) { // merge and cache the resulting header
            List<VCFHeader> headers = featureDataSources
                    .stream()
                    .map(ds -> (VCFHeader) ds.getHeader())
                    .filter(h -> h != null)
                    .collect(Collectors.toList());

            // Now merge the headers using htsjdk, which is pretty promiscuous, and which only works properly
            // because of the cross-dictionary validation done in validateSequenceDictionaries.
            mergedHeader = headers.size() > 1 ?
                    new VCFHeader(VCFUtils.smartMergeHeaders(headers, true)) :
                    headers.get(0);
        }
        return mergedHeader;
    }

    /**
     * Permanently close this data source, invalidating any open iteration over it, and making it invalid for future
     * iterations and queries.
     */
    @Override
    public void close() {
        validateFeatureSources();
        featureDataSources.forEach(dataSource -> dataSource.close());
    }

    /**
     * Pairwise validate the sequence dictionary from each source with the sequence dictionary from each other source,
     * since its possible for any given pair to have conflicts that are unique to that pair. GATKTool only validates
     * individual feature dictionaries against the reference dictionary, so we need to cross-validate them against
     * each other here before we merge them together.
     */
    private void validateSequenceDictionaries() {
        featureDataSources.forEach(
                baseFDS -> featureDataSources.stream()
                    .filter(secondaryFDS -> baseFDS != secondaryFDS)
                    .forEach(
                        secondaryFDS -> SequenceDictionaryUtils.validateDictionaries(
                            baseFDS.getName(),
                            baseFDS.getSequenceDictionary(),
                            secondaryFDS.getName(),
                            secondaryFDS.getSequenceDictionary(),
                            allowedSequenceDictionaryIncompatibilities,
                            false,
                            false)
                    )
        );
    }

    /**
     * Close the iterator currently open over the data sources, if there is one.
     */
    private void closeOpenIterationIfNecessary() {
        if (currentIterator != null) {
            currentIterator.close();
            currentIterator = null;
        }
    }

    /**
     * Wrap the sourceIterator in a CloseableIterator to make it usable as a MergingIterator source.
     */
    private CloseableIterator<VariantContext> getCloseableIteratorWrapper(final Iterator<VariantContext> sourceIterator) {
        validateFeatureSources();
        Utils.nonNull(sourceIterator);

        return new CloseableIterator<VariantContext>() {
            Iterator<VariantContext> delegateIterator = sourceIterator;
            @Override
            public void close() { delegateIterator = null; }

            @Override
            public boolean hasNext() {
                return delegateIterator.hasNext();
            }

            @Override
            public VariantContext next() {
                return delegateIterator.next();
            }
        };
    }

    /**
     * Invalidate the cached VCFHeader and SAMSequenceDictionary whenever the set of backing FeatureDataSources
     * changed in order to force them to be recreated and cached on demand.
     */
    private void invalidateCachedHeaderAndDictionary() {
        mergedHeader = null;
        mergedDictionary = null;
    }

    /**
     * Validate the requirement that one or more variant data sources have been added
     */
    private void validateFeatureSources() {
        if (featureDataSources == null || featureDataSources.isEmpty()) {
            throw new IllegalArgumentException("MultiVariantDataSource requires one or more Variant sources to be added");
        }
    }

}
