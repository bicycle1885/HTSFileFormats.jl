# AuxDataDict
# ===========
#
# Auxiliary data dictionary for BAM.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# This is not designed for very large dictionaries: time complexities in lookup
# and update operations are O(N).
immutable AuxDataDict <: Associative{String,Any}
    data::Vector{UInt8}
end

function Base.getindex(dict::AuxDataDict, tag)
    checkkeytag(tag)
    return _auxiliary(dict.data, 1, UInt8(tag[1]), UInt8(tag[2]))
end

#function Base.eltype(::Type{AuxDataDict})
#    return Tuple{String,Any}
#end

function Base.length(dict::AuxDataDict)
    return count_auxtags(dict.data, 1)
end

function Base.start(dict::AuxDataDict)
    return 1
end

function Base.done(dict::AuxDataDict, pos)
    return pos > length(dict.data)
end

function Base.next(dict::AuxDataDict, pos)
    data = dict.data
    tag = String([data[pos], data[pos+1]])
    pos, typ = getauxtype(data, pos + 2)
    pos, value = getauxdata(data, pos, typ)
    return (tag, value), pos
end


# Internals
# ---------

function _auxiliary(data::Vector{UInt8}, pos::Integer, t1::UInt8, t2::UInt8)
    p::Int = pos

    while p ≤ length(data) && (data[p] != t1 || data[p+1] != t2)
        p = next_tag_position(data, p)
    end

    if p > length(data)
        throw(KeyError(KeyTag(t1, t2)))
    else
        p, typ = getauxtype(data, p + 2)
        _, value = getauxdata(data, p, typ)
        return value
    end
end

function getauxtype(data::Vector{UInt8}, p::Int)
    t = data[p]
    if t == UInt8('B')
        return p + 2, Vector{auxtype[data[p+1]]}
    else
        return p + 1, auxtype[t]
    end
end

function getauxdata{T}(data::Vector{UInt8}, p::Int, ::Type{T})
    return p + sizeof(T), unsafe_load(Ptr{T}(pointer(data, p)))
end

function getauxdata(data::Vector{UInt8}, p::Int, ::Type{Char})
    return p + 1, Char(unsafe_load(pointer(data, p)))
end

function getauxdata{T}(data::Vector{UInt8}, p::Int, ::Type{Vector{T}})
    n = unsafe_load(Ptr{Int32}(pointer(data, p)))
    p += 4
    xs = Array(T, n)
    unsafe_copy!(pointer(xs), Ptr{T}(pointer(data, p)), n)
    return p + n * sizeof(T), xs
end

function getauxdata(data::Vector{UInt8}, p::Int, ::Type{String})
    dataptr = pointer(data, p)
    endptr = ccall(:memchr, Ptr{Void}, (Ptr{Void}, Cint, Csize_t),
                   dataptr, '\0', length(data) - p + 1)
    q = p + (endptr - dataptr) - 1
    return q + 2, String(data[p:q])
end

function count_auxtags(data::Vector{UInt8}, p::Int)
    count = 0
    while p ≤ length(data)
        count += 1
        p = next_tag_position(data, p)
    end
    return count
end

# Find the starting position of a next tag in `data` after `p`.
# `(data[p], data[p+1])` is supposed to be a current tag.
function next_tag_position(data::Vector{UInt8}, p::Int)
    typ = Char(data[p+2])
    p += 3
    if typ == 'A'
        p += 1
    elseif typ == 'c' || typ == 'C'
        p += 1
    elseif typ == 's' || typ == 'S'
        p += 2
    elseif typ == 'i' || typ == 'I'
        p += 4
    elseif typ == 'f'
        p += 4
    elseif typ == 'd'
        p += 8
    elseif typ == 'Z' || typ == 'H'
        while data[p] != 0x00  # NULL-terminalted string
            p += 1
        end
        p += 1
    elseif typ == 'B'
        eltyp = Char(data[p])
        elsize = eltyp == 'c' || eltyp == 'C'                 ? 1 :
                 eltyp == 's' || eltyp == 'S'                 ? 2 :
                 eltyp == 'i' || eltye == 'I' || eltyp == 'f' ? 4 :
                 error("unrecognized auxiliary type: ", eltyp)
        p += 1
        n = unsafe_load(Ptr{Int32}(pointer(data, p)))
        p += elsize * n
    else
        error("unrecognized auxiliary type: ", typ)
    end
    return p
end
