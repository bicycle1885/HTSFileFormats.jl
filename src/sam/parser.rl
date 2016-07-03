%%{
    machine samparser;

    action count_line {
        input.state.linenum += 1
    }

    action anchor {
        Ragel.@anchor!
    }

    action init {
        empty!(output.optional_fields)
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
    }

    action optfield {
        optfieldstr = Ragel.@ascii_from_anchor!()
        tag = KeyTag(optfieldstr[1], optfieldstr[2])
        typ = optfieldstr[4]
        output.optional_fields[tag] = optfieldstr[6:end]
    }

    action record {
        # need to anchor here for the next record
        Ragel.@anchor!
        Ragel.@yield ftargs
    }

    newline = '\n' >count_line;

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
    optfield = (tag ':' typ ':' value) >anchor %optfield;

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
        newline) >init %record;

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
