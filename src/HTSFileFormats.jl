module HTSFileFormats

export
    KeyTag,
    @tag_str,
    BAM,
    BAMRecord,
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
import BGZFStreams: BGZFStream
import DataStructures: OrderedDict

include("keytag.jl")
include("sam/sam.jl")
include("bam/bam.jl")

end # module
