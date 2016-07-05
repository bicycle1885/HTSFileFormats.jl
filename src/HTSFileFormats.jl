module HTSFileFormats

export
    KeyTag,
    @tag_str,
    BAM,
    BAMRecord,
    SAM,
    SAMRecord,
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
    header

import Bio
import Bio: Ragel
import Bio.Seq: DNASequence
import BGZFStreams: BGZFStream
import BufferedStreams   # this is necessary though I don't know why.
import BufferedStreams: BufferedInputStream

include("keytag.jl")
include("sam/sam.jl")
include("bam/bam.jl")

end # module
