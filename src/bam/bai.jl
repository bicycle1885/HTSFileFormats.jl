# BAI
# ===

typealias BinIndex Dict{UInt32,Vector{NTuple{2,VirtualOffset}}}
typealias LinearIndex Vector{VirtualOffset}

type BAI
    indexes::Vector{Tuple{BinIndex,LinearIndex}}
    n_no_coors::Nullable{UInt64}
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
            chunks = Vector{NTuple{2,VirtualOffset}}()
            for i in 1:n_chunks
                chunk_beg = read(input, UInt64)
                chunk_end = read(input, UInt64)
                push!(chunks, (chunk_beg, chunk_end))
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