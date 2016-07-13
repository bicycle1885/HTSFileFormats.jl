using HTSFileFormats
using Base.Test
using Bio.Seq
import BGZFStreams: BGZFStream

function testfile(filename)
    return Pkg.dir("HTSFileFormats", "test", filename)
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
        @test next_refname(rec) == "*"
        @test position(rec) == 0
        @test next_position(rec) == 0
        @test template_length(rec) == 0
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
        @test bin(rec) == 0
        @test mapping_quality(rec) == 0
        @test flag(rec) == 0
        @test next_refname(rec) == "*"
        @test next_refid(rec) == 0
        @test next_position(rec) == 0
        @test template_length(rec) == 0
        @test seqname(rec) == ""
        @test cigar(rec) == ""
        @test sequence(rec) == dna""
        @test qualities(rec) == UInt8[]
    end

    @testset "Reader" begin
        reader = open(testfile("bam1.bam"), BAM)
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
end
