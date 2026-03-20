# Nucleotide 2-bit codec
C functions to encode and decode nucleotides between ASCII character strings and packed 2-bit values, using scalar, BMI2, SSE4.1 and AVX2 x86 instruction.

## Features
- [`nucleotide_2bit_codec.h`](nucleotide_2bit_codec.h) / [`nucleotide_2bit_codec.c`](nucleotide_2bit_codec.c):
  - Functions encoding ASCII bases and decoding 2-bit bases, that are tied to specific instruction sets (e.g. `encode_bases_scalar()`, `encode_bases_bmi2()`, `decode_bases_avx2()`, etc.)
  - Functions encoding a fixed number ASCII bases and decoding a fixed number of 2-bit bases, that are not tied to any specific instruction set
  - Utility functions to convert between base count and encoded/unencoded byte count
- [`decoded_tables.h`](decoded_tables.h) / [`decoded_tables.c`](decoded_tables.c):
  - Tables of all 256 combinations of 4 nucleotides represented by ASCII characters "ACGT", "acgt", "ACGU" and "acgu", to be passed to decoding functions
- Aimed at 64-bit x86 architectures
- Compiles with GCC, clang and MSVC
- Benchmarked and tested

## What is encoded and decoded
This library can encode ASCII characters representing non-ambiguous IUPAC DNA/RNA nucleotides to the following 2-bit values and decode them back:
Nucleotide base    | 2-bit encoding
-------------------|---------------
`A`, `a`           | `0b00`
`C`, `c`           | `0b01`
`G`, `g`           | `0b10`
`T`, `t`, `U`, `u` | `0b11`

Other ASCII characters between `' '` (32) and `DEL` (127) are also encoded to one of these four 2-bit values, but will be decoded to the matching nucleotides in the table above, not to the original ASCII character. See [`ascii_to_encoded` in `nucleotide_2bit_codec.c`](https://github.com/badsami/nucleotide_2bit_codec/blob/main/nucleotide_2bit_codec.c#L56-L135) for more details.

Here are the 2-bit encodings of some non-nucleotide characters encountered in sequences found in FASTQ and FASTA files:
ASCII char | 2-bit encoding
-----------|---------------
`N`        | `0b00`
`.`        | `0b00`
`-`        | `0b01`
`_`        | `0b00`
`~`        | `0b00`

## Using
1. Copy, download or clone parts of the code in this repository
2. Compile with `-msse4.1`, `-mavx2`, `-mbmi2`, `-arch:AVX2`, etc. depending on your compiler and target hardware

Although the BMI2, SSE4.1 and AVX2 instructions are available on most x86 CPUs nowadays, be aware that they may not always be. One possible way to deal with this case is to perform runtime [CPU dispatching](https://johnnysswlab.com/cpu-dispatching-make-your-code-both-portable-and-fast/) based on the [`CPUID` x86 instruction](https://en.wikipedia.org/wiki/CPUID).

## Benchmarks & tests
The benchmarks and tests code (which is Windows-specific), along with instructions to build and run it, is available in [nucleotide_2bit_codec_benchmarks_and_tests](https://github.com/badsami/nucleotide_2bit_codec_benchmarks_and_tests).

Benchmarks encode/decode 1.00 MiB-worth of unencoded/encoded bases 1024 times, on a single core.  

The numbers reported below are from benchmarks compiled with MSVC using compilers flags `-O2 -arch:AVX2` and with clang using compiler flags `-O2 -mavx2 -mbmi2`, and run on an AMD Zen 3 5800H CPU (released January 2021) on Windows. For all 4 benchmark scenarii, measurements were within ±5% of each others. So the results reported below apply benchmarks compiled with both MSVC and clang.  
   
The benchmarks were ran 10 times each, without any open application running in the background. The best of the 10 runs are reported below.

#### Encoding
- CPU clock boost disabled
  Type   | Minimum bandwidth | Maximum bandwidth | Average bandwidth 
  -------|-------------------|-------------------|------------------
  Scalar |        3.43 GiB/s |        4.32 GiB/s |    **4.26 GiB/s**
  BMI2   |        7.71 GiB/s |       15.09 GiB/s |   **14.54 GiB/s**
  SSE4.1 |        9.10 GiB/s |       31.00 GiB/s |   **29.67 GiB/s**
  AVX2   |       13.37 GiB/s |       38.44 GiB/s |   **36.78 GiB/s**

- CPU clock boost enabled
  Type   | Minimum bandwidth | Maximum bandwidth | Average bandwidth
  -------|-------------------|-------------------|------------------
  Scalar |        4.81 GiB/s |        6.05 GiB/s |    **5.94 GiB/s**
  BMI2   |        9.29 GiB/s |       20.86 GiB/s |   **19.84 GiB/s**
  SSE4.1 |       15.62 GiB/s |       42.64 GiB/s |   **39.98 GiB/s**
  AVX2   |       17.25 GiB/s |       53.65 GiB/s |   **49.65 GiB/s**

#### Decoding
- CPU clock boost disabled
  Type   | Minimum bandwidth | Maximum bandwidth | Average bandwidth
  -------|-------------------|-------------------|------------------
  Scalar |        9.72 GiB/s |       13.45 GiB/s |   **13.30 GiB/s**
  SSE4.1 |       13.16 GiB/s |       22.24 GiB/s |   **21.87 GiB/s**
  AVX2   |       19.14 GiB/s |       41.37 GiB/s |   **40.16 GiB/s**

- CPU clock boost enabled
  Type   | Minimum bandwidth | Maximum bandwidth | Average bandwidth
  -------|-------------------|-------------------|------------------
  Scalar |       10.26 GiB/s |       18.74 GiB/s |   **17.24 GiB/s**
  SSE4.1 |       16.41 GiB/s |       30.14 GiB/s |   **28.75 GiB/s**
  AVX2   |       21.51 GiB/s |       53.65 GiB/s |   **51.42 GiB/s**

## Goals & motivation
- Promoting the 2-bit representation of nucleotides for operations repeatedly using/reading the same sequences (e.g. sequence alignment, indexing & searching reference genomes), shifting away from text
- Encouraging work on more compact data structures and accelerated algorithms for determined sequences
- Providing a reasonably-efficient and straightforward way to convert between ASCII character strings and packed 2-bit elements to build upon and improve
- Scaling sequence analysis for 3rd generation high-throughput sequencing technologies, which output reads longer than 2nd generation sequencing technologies
- Helping reduce processing times, memory usage, storage used space, and overall computing energy consumption

## License
The code in this repository is released in the public domain. You are free to use the code in this repository for any purpose.  

I only ask that you do not misrepresent the origin of this code: acknowledging or disclosing its origin is not required, but please do not claim that someone other than me wrote the original software, if inquired.