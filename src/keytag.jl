# KeyTag
# ======
#
# Two-byte tag used as keys for mapping.

immutable KeyTag <: AbstractString
    data::NTuple{2,UInt8}
end

"""
    KeyTag(s::AbstractString)
    KeyTag(t1::Char, t2::Char)
    KeyTag(t1::UInt8, t2::UInt8)

Create a tag string used as keys for mapping.
"""
function KeyTag(s::AbstractString)
    return convert(KeyTag, s)
end

function KeyTag(t1::Char, t2::Char)
    return KeyTag(UInt8(t1), UInt8(t2))
end

function KeyTag(t1::UInt8, t2::UInt8)
    return KeyTag((t1, t2))
end

function Base.convert(::Type{KeyTag}, str::AbstractString)
    if !isascii(str) || length(str) != 2
        throw(InexactError())
    end
    return KeyTag((UInt8(str[1]), UInt8(str[2])))
end

macro tag_str(str)
    return convert(KeyTag, str)
end

function Base.:(==)(x::KeyTag, y::KeyTag)
    return x.data == y.data
end

function Base.isless(x::KeyTag, y::KeyTag)
    return isless(x.data, y.data)
end

function Base.hash(tag::KeyTag, h::UInt)
    return hash(tag.data, h)
end

function Base.length(::KeyTag)
    return 2
end

function Base.endof(::KeyTag)
    return 2
end

function Base.start(::KeyTag)
    return 1
end

function Base.done(::KeyTag, i)
    return i > 2
end

function Base.next(tag::KeyTag, i)
    return Char(tag.data[i]), i + 1
end

function Base.checkbounds(tag::KeyTag, i::Integer)
    if 1 ≤ i ≤ endof(tag)
        return true
    end
    throw(BoundsError(i))
end

function Base.getindex(tag::KeyTag, i::Int)
    checkbounds(tag, i)
    return Char(tag.data[i])
end

function Base.show(io::IO, tag::KeyTag)
    write(io, "tag\"", tag.data[1], tag.data[2], '"')
    return
end
