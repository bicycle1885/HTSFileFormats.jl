# SAM Header
# ==========

# Write the SAM header into the output and return the number of written bytes.
function write_samheader(io::IO, header::OrderedDict{String})
    first = tag"HD"
    last = tag"CO"
    n = 0

    if haskey(header, first)
        for values in header[first]
            n += write_samheader_values(io, first, values)
        end
    end

    for (tag, values) in header; if tag != first && tag != last
        for values in header[tag]
            n += write_samheader_values(io, tag, values)
        end
    end; end

    if haskey(header, last)
        for values in header[last]
            n += write_samheader_values(io, last, values)
        end
    end

    return n
end

function write_samheader_values(io, headtag, values)
    n = 0
    n += write(io, '@', headtag)
    for (tag, value) in values
        n += write(io, '\t', tag, ':', string(value))
    end
    n += write(io, '\n')
    return n
end
