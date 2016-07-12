# BAI
# ===

typealias BinIndex Dict{UInt32,Vector{Chunk}}
typealias LinearIndex Vector{VirtualOffset}

const LinearWindowSize = 16 * 1024

type BAI
    indexes::Vector{Tuple{BinIndex,LinearIndex}}
    n_no_coors::Nullable{UInt64}
end

function chunks(index::BAI, refid::Integer, grange::UnitRange)
    if isempty(grange)
        return Chunk[]
    end
    binindex, linindex = index.indexes[refid]
    bins = reg2bins(first(grange), last(grange))
    offset = linindex[cld(first(grange), LinearWindowSize)]
    ret = Chunk[]
    for bin in bins
        if haskey(binindex, bin)
            for chunk in binindex[bin]
                if chunk.stop > offset
                    push!(ret, chunk)
                end
            end
        end
    end
    return sort!(ret)
end

function Base.read(input::IO, ::Type{BAI})
    # magic bytes
    B = read(input, UInt8)
    A = read(input, UInt8)
    I = read(input, UInt8)
    x = read(input, UInt8)
    if B != UInt8('B') || A != UInt8('A') || I != UInt8('I') || x != 0x01
        error("input is not a valid BAI file")
    end

    indexes = Tuple{BinIndex,LinearIndex}[]
    n_refs = read(input, Int32)
    for _ in 1:n_refs
        binindex = BinIndex()
        n_bins = read(input, Int32)
        for _ in 1:n_bins
            bin = read(input, UInt32)
            n_chunks = read(input, Int32)
            chunks = Vector{Chunk}()
            for i in 1:n_chunks
                chunk_beg::VirtualOffset = read(input, UInt64)
                chunk_end::VirtualOffset = read(input, UInt64)
                push!(chunks, Chunk(chunk_beg, chunk_end))
            end
            binindex[bin] = chunks
        end

        n_intvs = read(input, Int32)
        linindex = LinearIndex()
        for _ in 1:n_intvs
            push!(linindex, read(input, UInt64))
        end

        push!(indexes, (binindex, linindex))
    end

    if !eof(input)
        n_no_coors = Nullable(read(input, UInt64))
    else
        n_no_coors = Nullable{UInt64}()
    end

    return BAI(indexes, n_no_coors)
end

# Calculate bins overlapping a region [from, to] (one-based).
function reg2bins(from, to)
    bins = UInt32[]
    bin_start = 0
    for scale in 29:-3:14
        for k in ((from - 1) >> scale):((to - 1) >> scale)
            push!(bins, bin_start + k)
        end
        bin_start = 8 * bin_start + 1
    end
    return bins
end
