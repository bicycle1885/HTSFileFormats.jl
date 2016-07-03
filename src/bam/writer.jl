# BAM Writer
# ==========

type BAMWriter <: Bio.IO.AbstractWriter
    stream::BGZFStream
end

function Base.close(writer::BAMWriter)
    close(writer.stream)
end

function Base.write(writer::BAMWriter, header::BAMHeader)
    stream = writer.stream
    n = 0

    # magic bytes
    n += write(stream, "BAM\1")

    # SAM header
    buf = IOBuffer()
    l = write_samheader(buf, header.samheader)
    n += write(stream, Int32(l))
    n += write(stream, takebuf_array(buf))

    # reference sequences
    n += write(stream, Int32(length(header.refseqnames)))
    for (seqname, seqlen) in zip(header.refseqnames, header.refseqlens)
        namelen = length(seqname)
        n += write(stream, Int32(namelen + 1))
        n += write(stream, seqname, '\0')
        n += write(stream, Int32(seqlen))
    end

    return n
end

function Base.write(writer::BAMWriter, aln::BAMRecord)
    n = 0
    n += write(writer.stream, Int32(BAM_FIXED_FIELDS_BYTES + aln.datasize))
    n += unsafe_write(writer.stream, pointer_from_objref(aln), BAM_FIXED_FIELDS_BYTES)
    n += unsafe_write(writer.stream, pointer(aln.data), aln.datasize)
    return n
end
