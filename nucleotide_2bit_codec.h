#pragma once
#include "types.h"

// ASCII characters representing DNA/RNA nucleotides are always encoded to and decoded from the
// following 2-bit values:
// Nucleotide | 2-bit encoding
// -----------|---------------
// A, a       | 0b00
// C, c       | 0b01
// G, g       | 0b10
// T, t, U, u | 0b11
//
// Non-encoding ASCII characters, including ambiguous IUPAC nucleotide codes such as R (A or G) or B
// (C, G or T), between ' ' (32) and DEL (127) are also always converted to one of these 2-bit
// values. See "ascii_to_encoded" in nucleotide_2bit_codec.c for detailed encodings.
//
// That means, for instance, that "OFIL" will be encoded to 0b01101000 which might be decoded to
// AGGC or aggc, or that "...." will be encoded to 0b00000000 which might be decoded to AAAA or aaaa

////////////////////////////////////////////////////////////////////////////////////////////////////
//// Helpers
u64 base_count_to_encoded_byte_count(u64 base_count);
u64 base_count_to_decoded_byte_count(u64 base_count);
u64 encoded_byte_count_to_base_count(u64 base_count);
u64 decoded_byte_count_to_base_count(u64 base_count);


////////////////////////////////////////////////////////////////////////////////////////////////////
//// Encoding (8-bit ASCII chars to 2-bit elements)
u8 encode_base(char base);

void encode_bases_scalar(const char* restrict unencoded_bases,
                         u8*         restrict encoded_bases,
                         u64                  base_count);

void encode_bases_bmi2(const char* restrict unencoded_bases,
                       u8*         restrict encoded_bases,
                       u64                  base_count);

void encode_bases_sse4_1(const char* restrict unencoded_bases,
                         u8*         restrict encoded_bases,
                         u64                  base_count);

void encode_bases_avx2(const char* restrict unencoded_bases,
                       u8*         restrict encoded_bases,
                       u64                  base_count);


////////////////////////////////////////////////////////////////////////////////////////////////////
//// Decoding with specific reference bases (2-bit elements to 8-bit ASCII chars)
// decoded_table is a table of all 256 combinations of 4 unencoded bases, whatever those bases are.
// Tables for ACGT, acgt, ACGU, acgu are provided in decoded_tables.h.
char decode_base(u64 encoded_base, const char decoded_table[256 * 4]);

void decode_bases_x4(u64 encoded_bases_x4, char* decoded_bases, const char decoded_table[256 * 4]);

void decode_bases_scalar(const u8* restrict encoded_bases,
                         char*     restrict decoded_bases,
                         u64                base_count,
                         const char         decoded_table[256 * 4]);

void decode_bases_sse4_1(const u8* restrict encoded_bases,
                         char*     restrict decoded_bases,
                         u64                base_count,
                         const char         decoded_table[256 * 4]);

void decode_bases_avx2(const u8* restrict encoded_bases,
                       char*     restrict decoded_bases,
                       u64                base_count,
                       const char         decoded_table[256 * 4]);