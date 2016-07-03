# SAM Header
# ==========

# Read a SAM header from the input and return the result as a dictionary.
function read_samheader(io::IO)
    d = OrderedDict{String,Vector}()
    while !eof(io)
        line = readline(io)
        if line[1] != '@'
            error("'@' is expected")
        end
        tag = line[2:3]
        if line[4] != '\t'
            error("tab is expected")
        end
        rest = chomp(line[5:end])
        if !haskey(d, tag)
            d[tag] = []
        end
        if tag == "CO"
            push!(d[tag], rest)
        else
            push!(d[tag], read_samheader_values(rest))
        end
    end
    return d
end

function read_samheader_values(line)
    # keep the order of values
    ret = OrderedDict{String,String}()
    for pair in split(line, '\t')
        tag = pair[1:2]
        if pair[3] != ':'
            error("':' is expected")
        end
        ret[tag] = pair[4:end]
    end
    return ret
end

# Write the SAM header into the output and return the number of written bytes.
function write_samheader(io::IO, header::OrderedDict{String})
    first = "HD"
    last = "CO"
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

function write_samheader_values(io::IO, headtag::String, values::OrderedDict)
    n = 0
    n += write(io, '@', headtag)
    for (tag, value) in values
        n += write(io, '\t', tag, ':', string(value))
    end
    n += write(io, '\n')
    return n
end
