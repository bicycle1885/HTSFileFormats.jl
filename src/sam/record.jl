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
    return SAMRecord("", 0x0000, "", 0, 0, "*", "", 0, 0, "", UInt8[], Dict())
end

function Base.copy(rec::SAMRecord)
    return deepcopy(rec)
end
