# HTSFileFormats

**NOTE**: This package is merged in [Bio.jl](https://github.com/BioJulia/Bio.jl).

[![Build Status](https://travis-ci.org/bicycle1885/HTSFileFormats.jl.svg?branch=master)](https://travis-ci.org/bicycle1885/HTSFileFormats.jl)

## Synopsis

Getting all records:
```julia
using HTSFileFormats

reader = open("data.bam", BAM)
for record in reader
    # do something
end
close(reader)
```

Getting records overlapping a specific region:
```julia
for record in intersect(reader, "chr2", 10000:15000)
    # do something
end
```


## Performance Tips

### In-place reading

Iterating over records will copy each record for safety. Instead, you can read
data into a preallocated record object to improve the performance of reading:
```julia
record = BAMRecord()
while !eof(reader)
    read!(reader, record)
    # do something
end
```

### Type stability for optional fields

Optional fields in a record is not type-stable: the value type depends on the
tag key and hence its type is unpredictable for the Julia compiler. If you
know the type of a value in advance, you can specify it as a type annotation:
```julia
sum_nm = 0
for record in open(filename, BAM)
    sum_nm += record["NM"]::UInt8  # the type of `sum_nm` is stable now
end
```
