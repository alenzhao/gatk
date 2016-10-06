package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSpark;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.JoinStrategy;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.spark.ReadWalkerSpark;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import scala.Tuple2;
import scala.Tuple3;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        oneLineSummary = "BaseRecalibrator on Spark",
        programGroup = SparkProgramGroup.class
)
public class BaseRecalibratorSpark extends ReadWalkerSpark {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReference() { return true; }

    @Override
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return BaseRecalibrator.getStandardBQSRReadFilter(getHeaderForReads());
    }

    @Argument(fullName = "knownSites", shortName = "knownSites", doc = "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.", optional = false)
    private List<FeatureInput<Feature>> knownSites;

    @Argument(doc = "Path to save the final recalibration tables to.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputTablesPath = null;

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection bqsrArgs = new RecalibrationArgumentCollection();

    @Override
    protected void runTool( JavaSparkContext ctx ) {
        final RecalibrationReport bqsrReport = apply(getReads(ctx), getHeaderForReads(), knownSites, bqsrArgs);

        try ( final PrintStream reportStream = new PrintStream(BucketUtils.createFile(outputTablesPath, getAuthenticatedGCSOptions())) ) {
            RecalUtils.outputRecalibrationReport(reportStream, bqsrArgs, bqsrReport.getQuantizationInfo(), bqsrReport.getRecalibrationTables(), bqsrReport.getCovariates());
        }
    }

    public static RecalibrationReport apply(final JavaRDD<Tuple3<GATKRead, ReferenceContext, FeatureContext>> readsWithContext, final SAMFileHeader header, final List<FeatureInput<Feature>> knownSites, final RecalibrationArgumentCollection recalArgs ) {
        JavaRDD<RecalibrationTables> unmergedTables = readsWithContext.mapPartitions(readWithContextIterator -> {
            final BaseRecalibrationEngine bqsr = new BaseRecalibrationEngine(recalArgs, header);
            bqsr.logCovariatesUsed();

            while ( readWithContextIterator.hasNext() ) {
                final Tuple3<GATKRead, ReferenceContext, FeatureContext> readWithData = readWithContextIterator.next();
                final GATKRead read = readWithData._1();
                final ReferenceContext referenceContext = readWithData._2();
                final FeatureContext featureContext = readWithData._3();
                bqsr.processRead(read, referenceContext.getDataSource(), featureContext.getValues(knownSites));
            }
            // Need to wrap in ArrayList due to our current inability to serialize the return value of Arrays.asList() directly
            return new ArrayList<>(Arrays.asList(bqsr.getRecalibrationTables()));
        });

        final RecalibrationTables emptyRecalibrationTable = new RecalibrationTables(new StandardCovariateList(recalArgs, header));
        final RecalibrationTables combinedTables = unmergedTables.treeAggregate(emptyRecalibrationTable,
                RecalibrationTables::inPlaceCombine,
                RecalibrationTables::inPlaceCombine,
                Math.max(1, (int)(Math.log(unmergedTables.partitions().size()) / Math.log(2))));

        BaseRecalibrationEngine.finalizeRecalibrationTables(combinedTables);

        final QuantizationInfo quantizationInfo = new QuantizationInfo(combinedTables, recalArgs.QUANTIZING_LEVELS);

        final StandardCovariateList covariates = new StandardCovariateList(recalArgs, header);
        return RecalUtils.createRecalibrationReport(recalArgs.generateReportTable(covariates.covariateNames()), quantizationInfo.generateReportTable(), RecalUtils.generateReportTables(combinedTables, covariates));
    }
}
