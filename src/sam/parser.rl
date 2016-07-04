%%{
    machine samparser;

    action count_line {
        input.state.linenum += 1
    }

    action anchor {
        Ragel.@anchor!
    }

    action qname {
        Ragel.@copy_from_anchor!(output.name)
    }

    action flag {
        output.flag = Ragel.@int64_from_anchor!
    }

    action rname {
        Ragel.@copy_from_anchor!(output.refname)
    }

    action pos {
        output.pos = Ragel.@int64_from_anchor!
    }

    action mapq {
        output.mapq = Ragel.@int64_from_anchor!
    }

    action cigar {
        Ragel.@copy_from_anchor!(output.cigar)
    }

    action tlen {
        output.tlen = Ragel.@int64_from_anchor!
    }

    action seq {
        seqstr = Ragel.@ascii_from_anchor!()
        resize!(output.seq, length(seqstr))
        Bio.Seq.encode_copy!(output.seq, seqstr)
    }

    action qual {
        qualstr = Ragel.@ascii_from_anchor!()
        resize!(output.qual, length(qualstr))
        for i in 1:endof(qualstr)
            output.qual[i] = UInt8(qualstr[i]) - 33
        end
        empty!(output.optional_fields)
    }

    action optfield {
        optfieldstr = Ragel.@ascii_from_anchor!()
        tag = KeyTag(optfieldstr[1], optfieldstr[2])
        typ = optfieldstr[4]
        if typ == 'A'
            value = optfieldstr[6]
        elseif typ == 'Z'
            value = optfieldstr[6:end]
        elseif typ == 'H'
            value = parse_hexbytearray(optfieldstr[6:end])
        elseif typ == 'B'
            eltyp = auxtype[UInt8[optfieldstr[6]]]
            value = [parse(eltyp, x) for x in split(optfieldstr[8:end], ',')]
        else
            value = parse(auxtype[UInt8(typ)], optfieldstr[6:end])
        end
        output.optional_fields[tag] = value
    }

    action record {
        # need to anchor here for the next record
        Ragel.@anchor!
        Ragel.@yield ftargs
    }

    newline = '\n' >count_line;
    sign    = ('-' | '+')?;
    float   = sign? digit* '.'? digit+ ([eE] sign? digit+)?;

    qname = (graph - '@')+ >anchor %qname;
    flag  = digit+ >anchor %flag;
    rname = graph+ >anchor %rname;
    pos   = digit+ >anchor %pos;
    mapq  = digit+ >anchor %mapq;
    cigar = ('*' | (digit+ [MIDNSHPX=])+) >anchor %cigar;
    tlen  = ('-'? digit+) >anchor %tlen;
    seq   = ('*' | [A-Za-z=.]+) >anchor %seq;
    qual  = graph+ >anchor %qual;

    tag      = alpha alnum;
    typ      = [A-Za-z];
    value    = print+;
    optfield = (tag ':' (
        ("A:" graph) |
        ("Z:" print+) |
        ("H:" [0-9A-F]+) |
        ("B:" (',' float)+) |
        ("f:" float) |
        ([cCsSiI] ':' sign? digit+))) >anchor %optfield;

    record = (
        qname '\t'
        flag  '\t'
        rname '\t'
        pos   '\t'
        mapq  '\t'
        cigar '\t'
        rname '\t'
        pos   '\t'
        tlen  '\t'
        seq   '\t'
        qual ('\t' optfield)*
        newline) %record;

    main := record*;
}%%

%% write data;

Ragel.@generate_read!_function(
    "samparser",
    SAMReader,
    SAMRecord,
    begin
        %% write exec;
    end)
