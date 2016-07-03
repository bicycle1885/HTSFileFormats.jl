module HTSFileFormats

export
    @tag_str,
    BAM,
    BAMRecord,
    header

import Bio
import BGZFStreams: BGZFStream
import DataStructures: OrderedDict

include("keytag.jl")
include("sam/sam.jl")
include("bam/bam.jl")

end # module
