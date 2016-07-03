using HTSFileFormats
using Base.Test
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

@testset "BAM" begin
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
