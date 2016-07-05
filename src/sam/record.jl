# SAM Record
# ==========

type SAMRecord
    name::Bio.StringFields.StringField
    flag::UInt16
    refname::Bio.StringFields.StringField
    pos::Int64
    mapq::Int8
    cigar::Bio.StringFields.StringField
    next_refname::Bio.StringFields.StringField
    next_pos::Int64
    tlen::Int32
    seq::DNASequence
    qual::Vector{UInt8}
    optional_fields::Dict{String,Any}
end

function SAMRecord()
    return SAMRecord("", 0x0000, "*", 0, 0, "*", "*", 0, 0, "", UInt8[], Dict())
end

function Base.copy(rec::SAMRecord)
    return deepcopy(rec)
end

function seqname(r::SAMRecord)
    return r.name
end

function flag(r::SAMRecord)
    return r.flag
end

function refname(r::SAMRecord)
    return r.refname
end

function next_refname(r::SAMRecord)
    return r.next_refname
end

function Base.position(r::SAMRecord)
    return r.pos
end

function next_position(r::SAMRecord)
    return r.next_pos
end

function mapping_quality(r::SAMRecord)
    return d.mapq
end

function template_length(r::SAMRecord)
    return r.tlen
end

function cigar(r::SAMRecord)
    return r.cigar
end

function sequence(r::SAMRecord)
    return r.seq
end

function qualities(r::SAMRecord)
    return r.qual
end
