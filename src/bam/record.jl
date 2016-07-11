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

    # file header (shared)
    header::BAMHeader
end

function BAMRecord()
    return BAMRecord(
        -1, -1, 0, 0, 0, -1, -1, 0, UInt8[], 0, BAMHeader())
end

function Base.show(io::IO, rec::BAMRecord)
    println(summary(rec), ":")
    println(io, "reference name: ", refname(rec))
    println(io, "next reference name: ", next_refname(rec))
    println(io, "position: ", position(rec))
    println(io, "next position: ", next_position(rec))
    println(io, "bin: ", bits(bin(rec)))
    println(io, "mapping quality: ", mapping_quality(rec))
    println(io, "flag: ", flag(rec))
    println(io, "template length: ", template_length(rec))
    println(io, "sequence name: ", seqname(rec))
    println(io, "CIGAR string: ", cigar(rec))
    println(io, "sequence: ", sequence(rec))
    println(io, "base qualities: ", qualities(rec))
      print(io, "auxiliary: ", auxiliary(rec))
end

# the data size of fixed-length fields (.refid-.tlen)
const BAM_FIXED_FIELDS_BYTES = 32

"""
    refid(aln::BAMRecord)

Return the index of a reference sequence that `aln` is mapped onto.

The index is 1-based and will be 0 for an alignment without mapping position.
"""
function refid(aln::BAMRecord)
    return aln.refid + 1
end

function next_refid(aln::BAMRecord)
    return aln.next_refid + 1
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
        return r.header.refseqnames[id]
    end
end

function next_refname(r::BAMRecord)
    id = next_refid(r)
    if id == 0
        return "*"
    else
        return r.header.refseqnames[id]
    end
end

"""
    position(aln::BAMRecord)

Return the position of a mapped read.

The index is 1-based and will be 0 for an alignment without mapping position.
"""
function Base.position(aln::BAMRecord)
    return aln.pos + 1
end

function next_position(aln::BAMRecord)
    return aln.next_pos + 1
end

"""
    bin(aln::BAMRecord)

Return the bin of the alignment `aln`.
"""
function Base.bin(aln::BAMRecord)
    return UInt16(aln.bin_mq_nl >> 16)
end

"""
    mapping_quality(aln::BAMRecord)

Return the mapping quality of the alignment `aln`.
"""
function mapping_quality(aln::BAMRecord)
    return UInt8((aln.bin_mq_nl >> 8) & 0xff)
end

"""
    flag(aln::BAMRecord)

Return the flag of the alignment `aln`.
"""
function flag(aln::BAMRecord)
    return UInt16(aln.flag_nc >> 16)
end

"""
    template_length(aln::BAMRecord)

Return the template length of the alignment `aln`.
"""
function template_length(aln::BAMRecord)
    return aln.tlen
end

"""
    seqname(aln::BAMRecord)

Return the read name of the alignment `aln`.
"""
function seqname(aln::BAMRecord)
    # drop the last NUL character
    return unsafe_string(pointer(aln.data), max(seqname_length(aln) - 1, 0))
end

"""
    cigar_rle(aln::BAMRecord)

Return a run-length encoded tuple `(ops, lens)` of the CIGAR string.
See also `cigar`.
"""
function cigar_rle(aln::BAMRecord)
    offset = seqname_length(aln)
    ops = Bio.Align.Operation[]
    lens = Int[]
    for i in offset+1:4:offset+n_cigar_op(aln)*4
        x = UInt32(aln.data[i])         |
            UInt32(aln.data[i+1]) <<  8 |
            UInt32(aln.data[i+2]) << 16 |
            UInt32(aln.data[i+3]) << 24
        op = Bio.Align.Operation(x & 0x0f)
        len = x >> 4
        push!(ops, op)
        push!(lens, len)
    end
    return ops, lens
end

"""
    cigar(aln::BAMRecord)

Return a CIGAR string of the alignment `aln`. See also `cigar_rle`.
"""
function cigar(aln::BAMRecord)
    buf = IOBuffer()
    for (op, len) in zip(cigar_rle(aln)...)
        print(buf, len, Char(op))
    end
    return takebuf_string(buf)
end

"""
    sequence(aln::BAMRecord)

Return a DNA sequence of the alignment `aln`.
"""
function sequence(aln::BAMRecord)
    return decode_bamseq!(Bio.Seq.DNASequence(sequence_length(aln)), aln)
end

"""
    qualities(aln::BAMRecord)

Return base qualities of the alignment `aln`.
"""
function qualities(aln::BAMRecord)
    seqlen = sequence_length(aln)
    offset = seqname_length(aln) + n_cigar_op(aln) * 4 + cld(seqlen, 2)
    return [aln.data[i+offset] for i in 1:seqlen]
end

"""
    auxiliary(aln::BAMRecord)

Return a auxiliary data dictionary of the alignment `aln`.
"""
function auxiliary(aln::BAMRecord)
    seqlen = sequence_length(aln)
    offset = seqname_length(aln) + n_cigar_op(aln) * 4 + cld(seqlen, 2) + seqlen
    return AuxDataDict(aln.data[offset+1:end])
end

function Base.getindex(aln::BAMRecord, field::AbstractString)
    seqlen = sequence_length(aln)
    offset = seqname_length(aln) + n_cigar_op(aln) * 4 + cld(seqlen, 2) + seqlen
    return _auxiliary(aln.data, offset + 1, UInt8(field[1]), UInt8(field[2]))
end

# Return the length of the read name.
function seqname_length(aln)
    return aln.bin_mq_nl & 0xff
end

# Return the number of CIGAR operations.
function n_cigar_op(aln)
    return aln.flag_nc & 0xffff
end

# Return the length of the DNA sequence.
function sequence_length(aln)
    return aln.l_seq
end

# "=ACMGRSVTWYHKDBN" -> [0,16)
const bam_nucs = [
    Bio.Seq.DNA_Gap, Bio.Seq.DNA_A, Bio.Seq.DNA_C, Bio.Seq.DNA_M,
    Bio.Seq.DNA_G,   Bio.Seq.DNA_R, Bio.Seq.DNA_S, Bio.Seq.DNA_V,
    Bio.Seq.DNA_T,   Bio.Seq.DNA_W, Bio.Seq.DNA_Y, Bio.Seq.DNA_H,
    Bio.Seq.DNA_K,   Bio.Seq.DNA_D, Bio.Seq.DNA_B, Bio.Seq.DNA_N
]

# Decode the DNA sequence in a BAM alignment into a DNASequence.
function decode_bamseq!(seq, aln)
    seqlen = sequence_length(aln)
    @assert length(seq) == seqlen

    i = 2
    j = seqname_length(aln) + n_cigar_op(aln) * 4 + 1
    while i ≤ seqlen
        x = aln.data[j]
        seq[i-1] = bam_nucs[(x >> 4) + 1]
        seq[i  ] = bam_nucs[(x & 0x0f) + 1]
        i += 2
        j += 1
    end
    if isodd(seqlen)
        x = aln.data[j]
        seq[i-1] = bam_nucs[(x >> 4) + 1]
    end

    return seq
end
