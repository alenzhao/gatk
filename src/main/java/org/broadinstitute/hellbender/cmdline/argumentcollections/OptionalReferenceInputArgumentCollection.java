package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;

/**
 * An argument collection for use with tools that optionally accept a reference file as input.
 */
public final class OptionalReferenceInputArgumentCollection extends ReferenceInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence", optional = true)
    private String referenceFileName;

    private File referenceFile = null;

    @Override
    public File getReferenceFile() {
        if (null!=referenceFile) return referenceFile;
        if (null==referenceFileName) return null;
        referenceFile = new File(referenceFileName);
        return referenceFile;
    }

    @Override
    public String getReferenceFileName() {
        return referenceFileName;
    }
}
