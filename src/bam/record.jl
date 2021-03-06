# BAM Record
# ==========

# NOTE: The order and size of fields are important; the data will be loaded
# using memory copy. See specs for the details.
type BAMRecord
    refid::Int32
    pos::Int32
    bin_mq_nl::UInt32
    flag_nc::UInt32
    l_seq::Int32
    next_refid::Int32
    next_pos::Int32
    tlen::Int32

    # variable length data
    data::Vector{UInt8}

    # filled bytes of data (≤ length(.data))
    datasize::Int

    # reference sequence names (shared)
    refseqnames::Vector{String}
end

function BAMRecord()
    return BAMRecord(-1, -1, 0, 0, 0, -1, -1, 0, UInt8[], 0, String[])
end

# NOTE: this does not copy `refseqnames`.
function Base.copy(rec::BAMRecord)
    return BAMRecord(
        rec.refid,
        rec.pos,
        rec.bin_mq_nl,
        rec.flag_nc,
        rec.l_seq,
        rec.next_refid,
        rec.next_pos,
        rec.tlen,
        copy(rec.data),
        rec.datasize,
        rec.refseqnames)
end

function Base.show(io::IO, rec::BAMRecord)
    println(summary(rec), ":")
    println(io, "reference name: ", refname(rec))
    println(io, "next reference name: ", nextrefname(rec))
    println(io, "position: ", position(rec))
    println(io, "next position: ", nextposition(rec))
    println(io, "mapping quality: ", mappingquality(rec))
    println(io, "flag: ", flag(rec))
    println(io, "template length: ", templatelength(rec))
    println(io, "sequence name: ", seqname(rec))
    println(io, "CIGAR string: ", cigar(rec))
    println(io, "sequence: ", sequence(rec))
    println(io, "base qualities: ", qualities(rec))
      print(io, "optional fields: ", optinal_fields(rec))
end

# the data size of fixed-length fields (.refid-.tlen)
const BAM_FIXED_FIELDS_BYTES = 32

"""
    ismapped(rec::BAMRecord)

Return `true` if and only if `rec` is mapped to a reference sequence.
"""
function ismapped(rec::BAMRecord)
    return rec.pos != -1
end

"""
    refid(rec::BAMRecord)

Return the index of a reference sequence that `rec` is mapped onto.

The index is 1-based and will be 0 for an alignment without mapping position.
"""
function refid(rec::BAMRecord)
    return rec.refid + 1
end

function nextrefid(rec::BAMRecord)
    return rec.next_refid + 1
end

"""
    refname(rec::BAMRecord)

Return the name of a reference sequence that `rec` is mapped onto.

If `rec` is unmapped, it returns `"*"` like SAM records.
"""
function refname(r::BAMRecord)
    id = refid(r)
    if id == 0
        return "*"
    else
        return r.refseqnames[id]
    end
end

function nextrefname(r::BAMRecord)
    id = nextrefid(r)
    if id == 0
        return "*"
    else
        return r.refseqnames[id]
    end
end

"""
    position(rec::BAMRecord)

Return the position of a mapped read.

The index is 1-based and will be 0 for an alignment without mapping position.
"""
function Base.position(rec::BAMRecord)
    return rec.pos + 1
end

function nextposition(rec::BAMRecord)
    return rec.next_pos + 1
end

"""
    mappingquality(rec::BAMRecord)

Return the mapping quality of the alignment `rec`.
"""
function mappingquality(rec::BAMRecord)
    return UInt8((rec.bin_mq_nl >> 8) & 0xff)
end

"""
    flag(rec::BAMRecord)

Return the flag of the alignment `rec`.
"""
function flag(rec::BAMRecord)
    return UInt16(rec.flag_nc >> 16)
end

"""
    templatelength(rec::BAMRecord)

Return the template length of the alignment `rec`.
"""
function templatelength(rec::BAMRecord)
    return rec.tlen
end

"""
    seqname(rec::BAMRecord)

Return the read name of the alignment `rec`.
"""
function seqname(rec::BAMRecord)
    # drop the last NUL character
    return unsafe_string(pointer(rec.data), max(seqname_length(rec) - 1, 0))
end

"""
    cigar_rle(rec::BAMRecord)

Return a run-length encoded tuple `(ops, lens)` of the CIGAR string.
See also `cigar`.
"""
function cigar_rle(rec::BAMRecord)
    offset = seqname_length(rec)
    ops = Bio.Align.Operation[]
    lens = Int[]
    for i in offset+1:4:offset+n_cigar_op(rec)*4
        x = unsafe_load(Ptr{UInt32}(pointer(rec.data, i)))
        op = Bio.Align.Operation(x & 0x0f)
        push!(ops, op)
        push!(lens, x >> 4)
    end
    return ops, lens
end

"""
    cigar(rec::BAMRecord)

Return a CIGAR string of the alignment `rec`. See also `cigar_rle`.
"""
function cigar(rec::BAMRecord)
    buf = IOBuffer()
    for (op, len) in zip(cigar_rle(rec)...)
        print(buf, len, Char(op))
    end
    return takebuf_string(buf)
end

"""
    sequence(rec::BAMRecord)

Return a DNA sequence of the alignment `rec`.
"""
function sequence(rec::BAMRecord)
    seqlen = sequence_length(rec)
    data = Vector{UInt64}(cld(seqlen, 16))
    src::Ptr{UInt64} = pointer(rec.data, seqname_length(rec) + n_cigar_op(rec) * 4 + 1)
    for i in 1:endof(data)
        # copy data flipping high and low nybble
        x = unsafe_load(src, i)
        data[i] = (x & 0x0f0f0f0f0f0f0f0f) << 4 | (x & 0xf0f0f0f0f0f0f0f0) >> 4
    end
    return Bio.Seq.DNASequence(data, 1:seqlen, false)
end

"""
    qualities(rec::BAMRecord)

Return base qualities of the alignment `rec`.
"""
function qualities(rec::BAMRecord)
    seqlen = sequence_length(rec)
    offset = seqname_length(rec) + n_cigar_op(rec) * 4 + cld(seqlen, 2)
    return [rec.data[i+offset] for i in 1:seqlen]
end

function Base.getindex(rec::BAMRecord, tag::AbstractString)
    checkkeytag(tag)
    return getvalue(rec.data, auxdata_position(rec), UInt8(tag[1]), UInt8(tag[2]))
end

function Base.setindex!(rec::BAMRecord, val, tag::AbstractString)
    checkkeytag(tag)
    setvalue!(rec.data, auxdata_position(rec), val, UInt8(tag[1]), UInt8(tag[2]))
    return rec
end

function Base.delete!(rec::BAMRecord, tag::AbstractString)
    checkkeytag(tag)
    deletevalue!(rec.data, auxdata_position(rec), UInt8(tag[1]), UInt8(tag[2]))
    return rec
end

function Base.haskey(rec::BAMRecord, tag::AbstractString)
    checkkeytag(tag)
    return findtag(rec.data, auxdata_position(rec), UInt8(tag[1]), UInt8(tag[2])) > 0
end

function optinal_fields(rec::BAMRecord)
    return AuxDataDict(rec.data[auxdata_position(rec):rec.datasize])
end

function auxdata_position(rec)
    seqlen = sequence_length(rec)
    return seqname_length(rec) + n_cigar_op(rec) * 4 + cld(seqlen, 2) + seqlen + 1
end

# Return the right-most position of alignment.
function rightmost_position(rec::BAMRecord)
    return Int32(position(rec) + alignment_length(rec) - 1)
end

# Return the length of alignment.
function alignment_length(rec::BAMRecord)
    offset = seqname_length(rec)
    length::Int = 0
    for i in offset+1:4:offset+n_cigar_op(rec)*4
        x = unsafe_load(Ptr{UInt32}(pointer(rec.data, i)))
        op = Bio.Align.Operation(x & 0x0f)
        if Bio.Align.ismatchop(op) || Bio.Align.isdeleteop(op)
            length += x >> 4
        end
    end
    return length
end

# Return the length of the read name.
function seqname_length(rec)
    return rec.bin_mq_nl & 0xff
end

# Return the number of CIGAR operations.
function n_cigar_op(rec)
    return rec.flag_nc & 0xffff
end

# Return the length of the DNA sequence.
function sequence_length(rec)
    return rec.l_seq
end
