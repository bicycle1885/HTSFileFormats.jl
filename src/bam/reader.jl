# BAM Reader
# ==========

type BAMReader <: Bio.IO.AbstractParser
    stream::BGZFStream
    header::BAMHeader
end

function header(reader::BAMReader)
    return reader.header
end

function Base.eof(reader::BAMReader)
    return eof(reader.stream)
end

function Base.close(reader::BAMReader)
    close(reader.stream)
end

function Base.seek(reader::BAMReader, voffset::VirtualOffset)
    seek(reader.stream, voffset)
end

function Base.open(filename::AbstractString, ::Type{BAM})
    stream = BGZFStream(filename)

    # magic bytes
    B = read(stream, UInt8)
    A = read(stream, UInt8)
    M = read(stream, UInt8)
    x = read(stream, UInt8)
    if B != UInt8('B') || A != UInt8('A') || M != UInt8('M') || x != 0x01
        error("input was not a valid BAM file")
    end

    # SAM header
    textlen = read(stream, Int32)
    samheader = parse_samheader(read(stream, UInt8, textlen))

    # reference sequences
    refseqnames = String[]
    refseqlens = Int[]
    n_refs = read(stream, Int32)
    for _ in 1:n_refs
        namelen = read(stream, Int32)
        data = read(stream, UInt8, namelen)
        seqname = unsafe_string(pointer(data))
        seqlen = read(stream, Int32)
        push!(refseqnames, seqname)
        push!(refseqlens, seqlen)
    end

    return BAMReader(stream, BAMHeader(samheader, refseqnames, refseqlens))
end

function Base.read!(reader::BAMReader, aln::BAMRecord)
    datasize = read(reader.stream, Int32) - BAM_FIXED_FIELDS_BYTES
    unsafe_read(reader.stream, pointer_from_objref(aln), BAM_FIXED_FIELDS_BYTES)
    if length(aln.data) < datasize
        resize!(aln.data, datasize)
    end
    unsafe_read(reader.stream, pointer(aln.data), datasize)
    aln.datasize = datasize
    aln.header = reader.header
    return aln
end
