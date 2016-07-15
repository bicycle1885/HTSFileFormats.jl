module HTSFileFormats

export
    SAM,
    SAMHeader,
    SAMRecord,
    BAM,
    BAMRecord,
    BAI,
    Tabix,
    ismapped,
    refname,
    next_refname,
    refid,
    next_refid,
    position,
    next_position,
    bin,
    mapping_quality,
    flag,
    template_length,
    seqname,
    cigar,
    sequence,
    qualities,
    header,
    # SAM flags
    SAM_FLAG_PAIRED,
    SAM_FLAG_PROPER_PAIR,
    SAM_FLAG_UNMAP,
    SAM_FLAG_MUNMAP,
    SAM_FLAG_REVERSE,
    SAM_FLAG_MREVERSE,
    SAM_FLAG_READ1,
    SAM_FLAG_READ2,
    SAM_FLAG_SECONDARY,
    SAM_FLAG_QCFAIL,
    SAM_FLAG_DUP,
    SAM_FLAG_SUPPLEMENTARY


import Bio
import Bio: Ragel
import Bio.Seq: DNASequence
import BGZFStreams: BGZFStream, VirtualOffset, virtualoffset
import BufferedStreams   # this is necessary though I don't know why.
import BufferedStreams: BufferedInputStream

include("sam/sam.jl")
include("bam/bam.jl")
include("tabix/tabix.jl")

end # module
