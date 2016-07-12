# Chunk
# =====

# BGZF file chunk [.start, .stop).
immutable Chunk
    start::VirtualOffset
    stop::VirtualOffset
end

function Base.:(==)(chunk1::Chunk, chunk2::Chunk)
    return chunk1.start == chunk2.start && chunk1.stop == chunk2.stop
end

function Base.isless(chunk1::Chunk, chunk2::Chunk)
    return isless(chunk1.start, chunk2.start) ||
        (chunk1.start == chunk2.start &&
         isless(chunk1.stop, chunk2.stop))
end
