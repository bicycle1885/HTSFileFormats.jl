using HTSFileFormats
using Base.Test
using Bio.Seq
using BGZFStreams

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

    @testset "random access" begin
        reader = open(testfile("test.bam"), BAM)
        for (seqname, interval, expected) in [("chr1", 8000:10000, 20)]
            n = 0
            for rec in intersect(reader, seqname, interval)
                n += 1
            end
            @test n == expected
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
