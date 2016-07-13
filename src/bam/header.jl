# BAM Header
# ==========

type BAMHeader
    samheader::SAMHeader
    refseqnames::Vector{String}
    refseqlens::Vector{Int}
end

function BAMHeader()
    return BAMHeader(SAMHeader(), String[], Int[])
end

function Base.show(io::IO, header::BAMHeader)
    print(io, summary(header), ":")
    @assert length(header.refseqnames) == length(header.refseqlens)
    for (i, (name, len)) in enumerate(zip(header.refseqnames, header.refseqlens))
        println(io)
        print("  [", i, "]: ", name, " (length: ", len, ")")
    end
end

function Base.getindex(header::BAMHeader, i::Integer)
    return header.refseqnames[i], header.refseqlens[i]
end
