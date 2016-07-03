# SAM
# ===

"""
The SAM file format.
"""
immutable SAM <: Bio.IO.FileFormat end

include("flags.jl")
include("header.jl")
include("record.jl")
include("reader.jl")
include("parser.jl")
