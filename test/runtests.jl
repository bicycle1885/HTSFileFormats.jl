using HTSFileFormats
using Base.Test
using Bio.Seq
import BGZFStreams: BGZFStream

function testfile(filename)
    return Pkg.dir("HTSFileFormats", "test", filename)
end

@testset "KeyTag" begin
    @test KeyTag("BF") == KeyTag("BF")
    @test KeyTag("BF") != KeyTag("AP")
    @test tag"BF" == KeyTag("BF") == KeyTag('B', 'F') == KeyTag(0x42, 0x46)
    t = tag"BF"
    @test length(t) == endof(t) == 2
    @test t[1] == 'B'
    @test t[2] == 'F'
    @test collect(t) == ['B', 'F']
end

@testset "SAM" begin
    @testset "Reader" begin
        reader = open(testfile("sam1.sam"), SAM)
        @test isa(reader, HTSFileFormats.SAMReader)
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

@testset "BAM" begin
    @testset "Record" begin
        rec = BAMRecord()
        # default values
        @test refid(rec) == 0
        @test position(rec) == 0
        @test bin(rec) == 0
        @test mapping_quality(rec) == 0
        @test flag(rec) == 0
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
            writer = HTSFileFormats.BAMWriter(BGZFStream(path, "w"))
            write(writer, header(reader))
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
