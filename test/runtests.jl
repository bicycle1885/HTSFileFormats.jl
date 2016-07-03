using HTSFileFormats
using Base.Test
import BGZFStreams: BGZFStream

function testfile(filename)
    return Pkg.dir("HTSFileFormats", "test", filename)
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
