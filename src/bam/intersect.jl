# Intersect
# =========

type BAMIntersectionIterator
    reader::BAMReader
    chunks::Vector{Chunk}
    refid::Int
    interval::UnitRange{Int}
end

function Base.start(iter::BAMIntersectionIterator)
    if !isempty(iter.chunks)
        seek(iter.reader, first(iter.chunks).start)
    end
    rec = BAMRecord()
    return advance!(iter, rec, 1)
end

function Base.done(iter::BAMIntersectionIterator, i_rec)
    i, _ = i_rec
    return i > endof(iter.chunks)
end

function Base.next(iter::BAMIntersectionIterator, i_rec)
    i, rec = i_rec
    retrec = deepcopy(rec)
    return retrec, advance!(iter, rec, i)
end

function advance!(iter, rec, i)
    while i ≤ endof(iter.chunks)
        chunk = iter.chunks[i]
        while virtualoffset(iter.reader.stream) < chunk.stop
            read!(iter.reader, rec)
            if isoverlapping(rec, iter.refid, iter.interval)
                return i, rec
            end
        end
        i += 1
        if i ≤ endof(iter.chunks)
            seek(iter.reader, iter.chunks[i].start)
        end
    end
    return i, rec
end

function isoverlapping(rec, refid_, interval)
    leftmost = position(rec)

    if !ismapped(rec) || refid(rec) != refid_ || leftmost > last(interval) || isempty(interval)
        return false
    end

    # TODO: this might be slow because allocations are needed.
    rightmost = leftmost - 1
    ops, lens = cigar_rle(rec)
    for (op, len) in zip(ops, lens)
        if Bio.Align.ismatchop(op) || Bio.Align.isdeleteop(op)
            rightmost += len
        end
    end

    return rightmost ≥ first(interval)
end

function Base.intersect(reader::BAMReader, refid::Integer, interval::UnitRange)
    if isnull(reader.index)
        error("no index")
    end
    chunks = overlapchunks(get(reader.index), refid, interval)
    return BAMIntersectionIterator(reader, chunks, refid, interval)
end

function Base.intersect(reader::BAMReader, refname::AbstractString, interval::UnitRange)
    refid = findfirst(reader.refseqnames, refname)
    if refid == 0
        error("sequence name $refname is not in the header")
    end
    return intersect(reader, refid, interval)
end
