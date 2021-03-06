package org.broadinstitute.hellbender.utils.nio;

import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * NioBam holds a reference to a BAM file on Google Cloud Storage, and can give you
 * an RDD with the reads inside of it.
 *
 * The constructor will open the file to make sure it exists and we have access to it.
 * This is preferable to waiting until we're running on the cloud.
 *
 * Although we don't expect you to move these objects across computers, you can.
 */
public class NioBam implements Serializable {

    private static final long serialVersionUID = 1L;
    private final String bam;
    private final String index;
    private transient byte[] indexCache;

    /** Checks the files exists, then stores them. **/
    public NioBam(String gcsFilename, String indexGcsFilename) throws IOException {
        this.bam = gcsFilename;
        this.index = indexGcsFilename;
        init();
    }

    /** Finds the index file, then calls NioBam(bam, index). **/
    public NioBam(String gcsFilename) throws IOException {
        String indexFilename = gcsFilename + ".bai";
        if (!Files.exists(BucketUtils.getPathOnGcs(indexFilename))) {
            int i = gcsFilename.lastIndexOf('.');
            if (i>=0) {
                indexFilename = gcsFilename.substring(0, i) + ".bai";
            }
        }
        this.bam = gcsFilename;
        this.index = indexFilename;
        init();
    }

    private void init() throws IOException {
        Path bamPath = BucketUtils.getPathOnGcs(bam);
        Path bamIndexPath = BucketUtils.getPathOnGcs(index);
        if (!Files.exists(bamPath)) {
            throw new FileNotFoundException(bamPath.toString());
        }
        if (!Files.exists(bamIndexPath)) {
            throw new FileNotFoundException(bamIndexPath.toString());
        }
    }

    /** Parses the BAM file into SAMRecords. Will be distributed onto at least 'numPartitions' partitions. **/
    public JavaRDD<SAMRecord> getReads(JavaSparkContext ctx, int numPartitions) throws IOException {
        Path bamPath = BucketUtils.getPathOnGcs(bam);
        ChannelAsSeekableStream bamOverNIO = new ChannelAsSeekableStream(Files.newByteChannel(bamPath), bamPath.toString());
        final byte[] index  = getIndex();
        SeekableStream indexInMemory = new ByteArraySeekableStream(index);

        SamReader bam3 = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.LENIENT)
                .enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
                .open(SamInputResource.of(bamOverNIO).index(indexInMemory));
        List<QueryInterval> chunks = getAllChunksBalanced(bam3, numPartitions);

        // Ideally we'd get exactly the number of chunks the user is asking for, but until then...
        System.out.println("We got: " + chunks.size() + " chunks.");

        return ctx.parallelize(chunks, chunks.size()).flatMap(qi -> new ReadsIterable(bam, index, qi));
    }

    private synchronized byte[] getIndex() throws IOException {
        if (null!= indexCache) {
            return indexCache;
        }
        indexCache = Files.readAllBytes(BucketUtils.getPathOnGcs(index));
        return indexCache;
    }

    // this isn't very good yet, ideally we want just this number of query intervals, not per-contig.
    private static List<QueryInterval> getAllChunksBalanced(SamReader bam, int countPerContig) {
        List<QueryInterval> ret = new ArrayList<>();
        SAMFileHeader header = bam.getFileHeader();
        for (SAMSequenceRecord s : header.getSequenceDictionary().getSequences()) {
            ret.addAll(getChunksBalanced(bam, s.getSequenceIndex(), countPerContig));
        }
        return ret;
    }

    private static List<QueryInterval> getChunksBalanced(SamReader bam, int sequenceIndex, int retCount) {
        List<QueryInterval> ret = new ArrayList<>();
        BAMIndex index = bam.indexing().getIndex();
        SAMFileHeader header = bam.getFileHeader();
        SAMSequenceRecord s = header.getSequence(sequenceIndex);
        long totalLength = chunksLength(getChunks(index, sequenceIndex, 1, s.getSequenceLength() + 1));
        if (totalLength == 0) {
            return ret;
        }
        int sofar = 0;
        long targetLength = totalLength / retCount;
        int end = s.getSequenceLength();
        int step = s.getSequenceLength() / (100 * retCount);
        if (step < 1) step = 1;
        int start = 1;
        for (int j = step; j < end; j += step) {
            if (j > end) j = end;
            List<Chunk> candidate = getChunks(index, sequenceIndex, start, j);
            long size = chunksLength(candidate);
            if (size < targetLength) {
                // not big enough yet
                continue;
            }
            if (size > targetLength * 2) {
                // too large, search for a good separation point
                // TODO
            }
            // good, emit.
            ret.add(new QueryInterval(sequenceIndex, start, j + 1));
            start = j;
            sofar += size;
            if (ret.size() < retCount) {
                targetLength = (totalLength - sofar) / (retCount - ret.size());
            } else {
                targetLength = totalLength / retCount;
            }

        }
        return ret;
    }


    private static List<Chunk> getChunks(BAMIndex index, int sequenceIndex, int start, int endExcluded) {
        if (endExcluded <= start) return new ArrayList<>();
        BAMFileSpan span = index.getSpanOverlapping(sequenceIndex, start, endExcluded - 1);
        if (null == span) return new ArrayList<>();
        return span.getChunks();
    }

    private static long chunksLength(List<Chunk> chunks) {
        long totalLength = 0;
        for (Chunk c : chunks) {
            totalLength += chunkSize(c);
        }
        return totalLength;
    }

    private static long chunkSize(Chunk c) {
        long start = c.getChunkStart() >> 16;
        long end = (c.getChunkEnd() >> 16) + (c.getChunkEnd() & 0xffff);
        return (end - start);
    }

}
