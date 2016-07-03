# BAM Header
# ==========

type BAMHeader
    samheader::OrderedDict{String,Any}
    refseqnames::Vector{String}
    refseqlens::Vector{Int}
end

function BAMHeader()
    return BAMHeader(OrderedDict{String,Any}())
end
