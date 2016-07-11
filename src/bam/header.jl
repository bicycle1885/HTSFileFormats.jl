# BAM Header
# ==========

type BAMHeader
    samheader::Dict{String,Any}
    refseqnames::Vector{String}
    refseqlens::Vector{Int}
end

function BAMHeader()
    return BAMHeader(Dict{String,Any}(), String[], Int[])
end
