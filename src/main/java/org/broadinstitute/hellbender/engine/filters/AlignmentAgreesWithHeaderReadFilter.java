package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;

/**
 * Checks to ensure that the alignment of each read makes sense based on the contents of the header.
 */
public final class AlignmentAgreesWithHeaderReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1l;

    public AlignmentAgreesWithHeaderReadFilter( ) {}

    public AlignmentAgreesWithHeaderReadFilter( final SAMFileHeader header ) {
        super.setHeader(header);
    }

    @Override
    public boolean test( GATKRead read ) {
        return ReadUtils.alignmentAgreesWithHeader(samHeader, read);
    }
}
