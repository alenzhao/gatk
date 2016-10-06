package org.broadinstitute.hellbender.utils.spark;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple3;

import javax.annotation.Nullable;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

public abstract class ReadWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * This number controls the size of the cache for our FeatureInputs
     * (specifically, the number of additional bases worth of overlapping records to cache when querying feature sources).
     */
    public static final int FEATURE_CACHE_LOOKAHEAD = 1_000;

    private static final int SHARD_SIZE = 10000;
    private static final int SHARD_PADDING = 1000;

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    private FeatureManager features; // TODO: move up to GATKSparkTool?

    @Override
    protected void runPipeline(JavaSparkContext sparkContext) {
        initializeFeatures();
        super.runPipeline(sparkContext);
    }

    void initializeFeatures() {
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD);
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * Loads reads and the corresponding reference and features into a {@link JavaRDD} for the intervals specified.
     *
     * If no intervals were specified, returns all the reads.
     *
     * @return all reads from as a {@link JavaRDD}, bounded by intervals if specified.
     */
    public JavaRDD<Tuple3<GATKRead, ReferenceContext, FeatureContext>> getReads(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = hasIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        // use unpadded shards (padding is only needed for reference bases)
        final List<ShardBoundary> intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, SHARD_SIZE, SHARD_SIZE, 0, sequenceDictionary).stream())
                .collect(Collectors.toList());
        JavaRDD<Shard<GATKRead>> shardedReads = SparkUtils.shard(ctx, getReads(), GATKRead.class, sequenceDictionary, intervalShards, shuffle);
        Broadcast<ReferenceMultiSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedReads.flatMap(getReadsFunction(bReferenceSource, bFeatureManager, sequenceDictionary));
    }

    private static FlatMapFunction<Shard<GATKRead>, Tuple3<GATKRead, ReferenceContext, FeatureContext>> getReadsFunction(
            Broadcast<ReferenceMultiSource> bReferenceSource, Broadcast<FeatureManager> bFeatureManager,
            SAMSequenceDictionary sequenceDictionary) {
        return (FlatMapFunction<Shard<GATKRead>, Tuple3<GATKRead, ReferenceContext, FeatureContext>>) shard -> {
            // get reference bases for this shard (padded)
            SimpleInterval paddedInterval = shard.getInterval().expandWithinContig(SHARD_PADDING, sequenceDictionary);
            ReferenceDataSource reference = bReferenceSource == null ? null :
                    new ReferenceMemorySource(bReferenceSource.getValue().getReferenceBases(null, paddedInterval), sequenceDictionary);
            FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();

            Iterator<Tuple3<GATKRead, ReferenceContext, FeatureContext>> transform = Iterators.transform(shard.iterator(), new Function<GATKRead, Tuple3<GATKRead, ReferenceContext, FeatureContext>>() {
                @Nullable
                @Override
                public Tuple3<GATKRead, ReferenceContext, FeatureContext> apply(@Nullable GATKRead read) {
                    final SimpleInterval readInterval = getReadInterval(read);
                    return new Tuple3<>(read, new ReferenceContext(reference, readInterval), new FeatureContext(features, readInterval));
                }
            });
            // only include reads that start in the shard
            return () -> Iterators.filter(transform, r -> r._1().getStart() >= shard.getStart());
        };
    }

    /**
     * Returns an interval for the read.
     * Note: some walkers must be able to work on any read, including those whose coordinates do not form a valid SimpleInterval.
     * So here we check this condition and create null intervals for such reads.
     */
    static SimpleInterval getReadInterval(final GATKRead read) {
        return !read.isUnmapped() && SimpleInterval.isValid(read.getContig(), read.getStart(), read.getEnd()) ? new SimpleInterval(read) : null;
    }
}
