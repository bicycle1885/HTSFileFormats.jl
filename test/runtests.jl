using HTSFileFormats
using Base.Test
using Bio.Seq
using BGZFStreams

function testfile(filename)
    return Pkg.dir("HTSFileFormats", "test", filename)
end

function rand_interval(intval)
    x = rand(intval)
    y = rand(intval)
    if x < y
        return x:y
    else
        return y:x
    end
end

@testset "AuxDataDict" begin
    AuxDataDict = HTSFileFormats.AuxDataDict

    dict = AuxDataDict()
    @test length(dict) == 0
    @test isempty(dict)
    @test_throws KeyError dict["NM"]

    dict = AuxDataDict(
        "X1" => UInt8(1),
        "X2" => UInt16(2),
        "X3" => UInt32(3),
        "X4" => Int8(4),
        "X5" => Int16(5),
        "X6" => Int32(6),
        "X7" => Float32(7),
        "X8" => "eight",
        "X9" => Int32[9])
    @test length(dict) == 9
    @test !isempty(dict)
    @test dict["X1"] === UInt8(1)
    @test dict["X2"] === UInt16(2)
    @test dict["X3"] === UInt32(3)
    @test dict["X4"] === Int8(4)
    @test dict["X5"] === Int16(5)
    @test dict["X6"] === Int32(6)
    @test dict["X7"] === Float32(7)
    @test dict["X8"] == "eight"
    @test typeof(dict["X8"]) == String
    @test dict["X9"] == Int32[9]
    @test typeof(dict["X9"]) == Vector{Int32}

    dict = AuxDataDict("NM" => 0x01, "XY" => Int32(100), "XZ" => [0x11, 0x23])
    @test length(dict) == 3
    @test dict["NM"] === 0x01
    @test dict["XY"] === Int32(100)
    @test dict["XZ"] == [0x11, 0x23]
    @test eltype(dict["XZ"]) == UInt8

    dict = AuxDataDict("NM" => 0x01, "MD" => "8T1T39")
    @test length(dict) == 2
    @test dict["NM"] === 0x01
    @test dict["MD"] == "8T1T39"
    dict["NM"] = 0x00
    @test dict["NM"] === 0x00
    dict["MD"] = "50"
    @test dict["MD"] == "50"
    @test collect(dict) == ["NM" => 0x00, "MD" => "50"]
    dict["XY"] = "foobar"
    @test dict["XY"] == "foobar"
    @test collect(dict) == ["NM" => 0x00, "MD" => "50", "XY" => "foobar"]
    delete!(dict, "NM")
    @test length(dict) == 2
    @test collect(dict) == ["MD" => "50", "XY" => "foobar"]
end

@testset "SAM" begin
    @testset "SAMHeader" begin
        h = SAMHeader()
        @test isa(h, Associative)
        @test isempty(h)
        h["HD"] = Dict("VN" => "100.100", "SO" => "unknown")
        @test length(h) == 1
        @test h["HD"]["VN"] == "100.100"
        h["CO"] = ["comment1", "comment2"]
        @test length(h) == 2
        @test h["CO"] == ["comment1", "comment2"]
        delete!(h, "CO")
        @test length(h) == 1
    end

    @testset "Record" begin
        rec = SAMRecord()
        @test !ismapped(rec)
        # default values
        @test seqname(rec) == ""
        @test flag(rec) == 0x0000
        @test refname(rec) == "*"
        @test nextrefname(rec) == "*"
        @test position(rec) == 0
        @test nextposition(rec) == 0
        @test templatelength(rec) == 0
        @test cigar(rec) == "*"
        @test sequence(rec) == dna""
        @test qualities(rec) == UInt8[]
    end

    @testset "Reader" begin
        reader = open(testfile("sam1.sam"), SAM)
        @test isa(reader, HTSFileFormats.SAMReader)
        n = 0
        aln = SAMRecord()
        read!(reader, aln)
        n += 1
        read!(reader, aln)
        n += 1
        @test aln["XT"] === 'U'
        @test aln["NM"] === Int32(5)
        while !eof(reader)
            read!(reader, aln)
            n += 1
        end
        close(reader)
        @test n == 200
    end

    @testset "Round trip" begin
        mktemp() do path, io
            # copy
            reader = open(testfile("sam1.sam"), SAM)
            writer = HTSFileFormats.SAMWriter(io)
            write(writer, header(reader))
            aln = SAMRecord()
            while !eof(reader)
                read!(reader, aln)
                write(writer, aln)
            end
            close(reader)
            close(writer)

            # read again
            reader = open(path, SAM)
            n = 0
            aln = SAMRecord()
            while !eof(reader)
                read!(reader, aln)
                n += 1
            end
            close(reader)
            @test n == 200
        end
    end
end

@testset "BAM" begin
    @testset "Record" begin
        rec = BAMRecord()
        @test !ismapped(rec)
        # default values
        @test refname(rec) == "*"
        @test refid(rec) == 0
        @test position(rec) == 0
        @test mappingquality(rec) == 0
        @test flag(rec) == 0
        @test nextrefname(rec) == "*"
        @test nextrefid(rec) == 0
        @test nextposition(rec) == 0
        @test templatelength(rec) == 0
        @test seqname(rec) == ""
        @test cigar(rec) == ""
        @test sequence(rec) == dna""
        @test qualities(rec) == UInt8[]

        @test HTSFileFormats.rightmost_position(rec) === Int32(-1)
        @test HTSFileFormats.alignment_length(rec) === 0

        rec = BAMRecord()
        @test !haskey(rec, "MN")
        rec["MN"] = 0x01
        @test rec["MN"] === 0x01
        @test haskey(rec, "MN")
        delete!(rec, "MN")
        @test !haskey(rec, "MN")
    end

    @testset "Reader" begin
        reader = open(testfile("bam1.bam"), BAM)
        h = header(reader)
        @test h["SQ"] == [Dict("SN" => "1", "LN" => "239940")]
        @test h["PG"] == [Dict("ID" => "bwa", "PN" => "bwa", "VN" => "0.6.2-r126")]
        @test isa(reader, HTSFileFormats.BAMReader)
        n = 0
        aln = BAMRecord()
        read!(reader, aln)
        n += 1
        read!(reader, aln)
        n += 1
        @test aln["XT"] === 'U'
        # TODO: check type
        @test aln["NM"] == 5
        while !eof(reader)
            read!(reader, aln)
            n += 1
        end
        close(reader)
        @test n == 200
    end

    @testset "Round trip" begin
        mktemp() do path, _
            # copy
            reader = open(testfile("bam1.bam"), BAM)
            writer = HTSFileFormats.BAMWriter(
                BGZFStream(path, "w"),
                header(reader, true))
            aln = BAMRecord()
            while !eof(reader)
                read!(reader, aln)
                write(writer, aln)
            end
            close(reader)
            close(writer)

            # read again
            reader = open(path, BAM)
            n = 0
            aln = BAMRecord()
            while !eof(reader)
                read!(reader, aln)
                n += 1
            end
            close(reader)
            @test n == 200
        end
    end

    @testset "random access" begin
        reader = open(testfile("test.bam"), BAM)
        for (name, interval, expected) in [
                ("chr1", 1000:10000, 21),
                ("chr1", 8000:10000, 20),
                ("chr1", 766_000:800_000, 142),
                ("chr1", 786_000:800_000, 1),
                ("chr1", 796_000:800_000, 0)]
            n = 0
            for rec in intersect(reader, name, interval)
                n += 1
            end
            @test n == expected
        end

        for n in 1:50
            name = "chr1"
            refid = 1
            interval = rand_interval(1:1_000_000)

            expected = BAMRecord[]
            seekstart(reader)
            for rec in reader
                if HTSFileFormats.isoverlapping(rec, refid, interval)
                    push!(expected, rec)
                end
            end

            actual = BAMRecord[]
            for rec in intersect(reader, name, interval)
                push!(actual, rec)
            end

            @test map(seqname, actual) == map(seqname, expected)
        end
    end
end

@testset "Tabix" begin
    bedfile = testfile("knownGene.part.bed.gz")
    indexfile = bedfile * ".tbi"

    @testset "load" begin
        index = Tabix(indexfile)
        @test isa(index, Tabix)
        @test index.format === Int32(0x10000)
        @test index.columns === (1, 2, 3)
        @test index.meta === '#'
        @test index.skip === 0
        @test index.names == ["chr1", "chr2", "chr3"]
        @test length(index.indexes) == 3
    end

    @testset "chunks" begin
        stream = BGZFStream(bedfile)
        index = Tabix(indexfile)

        # Expected values were counted using the tabix tool as follows:
        #   $ tabix test/knownGene.part.bed.gz chr1:5,000,000-10,000,000 | wc -l
        #   $     64
        for (seqname, interval, expected) in [
                ("chr1", 5_000_000:10_000_000, 64),
                ("chr1", 6_000_000:10_000_000, 54),
                ("chr2", 5_000_000:10_000_000, 81),
                ("chr3", 9_000_000:10_000_000, 15),
                ("chr3", 9_800_000:10_000_000,  3),
                ("chr3", 9_900_000:10_000_000,  0)]
            n = 0
            for chunk in overlapchunks(index, seqname, interval)
                seek(stream, chunk)
                while virtualoffset(stream) in chunk
                    line = readline(stream)
                    values = split(chomp(line), '\t')
                    nm = values[index.columns[1]]
                    @test nm == seqname
                    int = parse(Int, values[index.columns[2]]):parse(Int, values[index.columns[3]])-1
                    if !isempty(intersect(int, interval))
                        n += 1
                    elseif first(int) > last(interval)
                        break
                    end
                end
            end
            @test n == expected
        end
    end
end
