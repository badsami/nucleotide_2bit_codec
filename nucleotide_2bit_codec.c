#include "types.h"

#include <immintrin.h>

//  ______________________________________________
// | SIMD instruction set | Max char per register |
// | ---------------------|-----------------------|
// | None        (x86)    | 8  (u64)              |
// | SSE1 - SSE4 (x86)    | 16 (__m128i)          |
// | AVX1 - AVX2 (x86)    | 32 (__m256i)          |
//  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

// Encoding ASCII characters
//              _____________________________________________________
//             |                   ASCII                   | Encoded |
//  ___________|‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾|‾‾‾‾‾‾‾‾‾|
// | Character | Hexadecimal | Decimal | Octal | Binary    | Binary  |
// |‾‾‾‾‾‾‾‾‾‾‾|‾‾‾‾‾‾‾‾‾‾‾‾‾|‾‾‾‾‾‾‾‾‾|‾‾‾‾‾‾‾|‾‾‾‾‾‾‾‾‾‾‾|‾‾‾‾‾‾‾‾‾|
// | A         | 0x41        | 65      | 0101  | 0b1000001 | 0b00    |
// | C         | 0x43        | 67      | 0103  | 0b1000011 | 0b01    |
// | G         | 0x47        | 71      | 0107  | 0b1000111 | 0b10    |
// | T         | 0x54        | 84      | 0124  | 0b1010100 | 0b11    |
// | U         | 0x55        | 85      | 0125  | 0b1010101 | 0b11    |
//  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

////////////////////////////////////////////////////////////////////////////////////////////////////
//// Helpers
u64 base_count_to_encoded_byte_count(u64 base_count)
{
  // 8-bit (1 byte) char are encoded into packed 2-bit (0.25 byte) elements, so base_count
  // unencoded chars need ceiled(base_count / 4) encoded bytes
  return (base_count + 3) / 4;
}


u64 base_count_to_decoded_byte_count(u64 base_count)
{
  return base_count;
}


u64 encoded_byte_count_to_base_count(u64 byte_count)
{
  return byte_count * 4;
}


u64 decoded_byte_count_to_base_count(u64 byte_count)
{
  return byte_count;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//// Encoding (8-bit ASCII chars to 2-bit elements)
#define LOWEST_VALID_ASCII_CHAR ' '
#define HIGHEST_VALID_ASCII_CHAR '_'
static const u8 ascii_to_encoded[64] =
{
  // Encoding DNA nucleotides
  ['A'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['C'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['G'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['T'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['U'  - LOWEST_VALID_ASCII_CHAR] = 0b11,

  // Non-encoding
  //
  // Both scalar and vectorized versions of functions encoding more than one base from text to
  // binary encode bits 2 and 1 of each 8-bit element with
  // 
  //   base ^ (base >> 1)
  //
  // then pack pairs of encoded bits together. For actual nucleotides, that always yields the binary
  // value specified above. Otherwise, that always yields the binary value below
  [' '  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['!'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['"'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['#'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['$'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['%'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['&'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['\'' - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['('  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  [')'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['*'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['+'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  [','  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['-'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['.'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['/'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['0'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['1'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['2'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['3'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['4'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['5'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['6'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['7'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['8'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['9'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  [':'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  [';'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['<'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['='  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['>'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['?'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['@'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['B'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['D'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['E'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['F'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['H'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['I'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['J'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['K'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['L'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['M'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['N'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['O'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['P'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['Q'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['R'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['S'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['V'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['W'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['X'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['Y'  - LOWEST_VALID_ASCII_CHAR] = 0b10,
  ['Z'  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['['  - LOWEST_VALID_ASCII_CHAR] = 0b11,
  ['\\' - LOWEST_VALID_ASCII_CHAR] = 0b01,
  [']'  - LOWEST_VALID_ASCII_CHAR] = 0b01,
  ['^'  - LOWEST_VALID_ASCII_CHAR] = 0b00,
  ['_'  - LOWEST_VALID_ASCII_CHAR] = 0b00
};

u8 encode_base(char base)
{
  // Handles ASCII characters from ' ' (32) to DEL (127)
  u64 idx = (base & 0b01011111) - LOWEST_VALID_ASCII_CHAR;
  return ascii_to_encoded[idx];
}


void encode_bases_scalar(const char* restrict unencoded_bases,
                         u8*         restrict encoded_bases,
                         u64                  base_count)
{
  const char*       unencoded_end          = unencoded_bases + base_count;
  const char* const unencoded_bases_x4_end = unencoded_bases + (base_count & ~3);
  const char* const unencoded_bases_x8_end = unencoded_bases + (base_count & ~7);
  
  // The following loops will likely be vectorized using the largest SIMD instruction set allowed,
  // unless a pragma directive (such as #pragma loop(no_vector) for MSVC), attribute or compiler
  // flag explicitly forbids it

  while (unencoded_bases < unencoded_bases_x8_end)
  {
    u64 bases_x8        = *(u64*)unencoded_bases;
    u64 bits_21_encoded = bases_x8 ^ (bases_x8 >> 1);

    u16 base_0 = (bits_21_encoded >>  1) &        3;
    u16 base_1 = (bits_21_encoded >>  7) & (3 <<  2);
    u16 base_2 = (bits_21_encoded >> 13) & (3 <<  4);
    u16 base_3 = (bits_21_encoded >> 19) & (3 <<  6);
    u16 base_4 = (bits_21_encoded >> 25) & (3 <<  8);
    u16 base_5 = (bits_21_encoded >> 31) & (3 << 10);
    u16 base_6 = (bits_21_encoded >> 37) & (3 << 12);
    u16 base_7 = (bits_21_encoded >> 43) & (3 << 14);

    *(u16*)encoded_bases = base_0 | base_1 | base_2 | base_3 | base_4 | base_5 | base_6 | base_7;

    unencoded_bases += 8;
    encoded_bases   += 2;
  }

  while (unencoded_bases < unencoded_bases_x4_end)
  {
    u32 bases_x4 = *(u32*)unencoded_bases;
    u32 bits_21_encoded = bases_x4 ^ (bases_x4 >> 1);
    u8 base_0 = (bits_21_encoded >>  1) &       3;
    u8 base_1 = (bits_21_encoded >>  7) & (3 << 2);
    u8 base_2 = (bits_21_encoded >> 13) & (3 << 4);
    u8 base_3 = (bits_21_encoded >> 19) & (3 << 6);

    *encoded_bases = base_0 | base_1 | base_2 | base_3;

    unencoded_bases += 4;
    encoded_bases   += 1;
  }

  if (unencoded_bases < unencoded_end)
  {
    // At most 3 unencoded bases remain
    u8 packed = 0;
    u8 shift  = 0;
    do
    {
      const u8 encoded_base = encode_base(*unencoded_bases);
      packed |= encoded_base << shift;

      unencoded_bases += 1;
      shift           += 2;
    }
    while (unencoded_bases < unencoded_end);

    *encoded_bases = packed;
  }
}


void encode_bases_bmi2(const char* restrict unencoded_bases,
                       u8*         restrict encoded_bases,
                       u64                  base_count)
{
  const char* const unencoded_bases_x8_end  = unencoded_bases + (base_count & ~7);
  const char* const unencoded_bases_x32_end = unencoded_bases + (base_count & ~31);

  while (unencoded_bases < unencoded_bases_x32_end)
  {
    // Encode bases in bits 2 (most significant bit) and 1 (least significant bit) of each byte.
    // Then extract and pack them in 16 bits. Example with "AAACAGAT":
    //
    //                     < most significant bit/byte                 least significant bit/byte >
    //  ____________________________________________________________________________________________
    // | Base (x)         | 'T'      'A'      'G'      'A'      'C'      'A'      'A'      'A'      |
    // | Binary           | 01010100 01000001 01000111 01000001 01000011 01000001 01000001 01000001 |
    // | x ^ (x >> 1)     | 01111110 01100001 01100100 01100001 01100010 01100001 01100001 01100001 |
    // | _pext_u64() mask | 00000110 00000110 00000110 00000110 00000110 00000110 00000110 00000110 |
    // | _pext_u64()      |                                                       11001000 01000000 |
    //  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
    u64 bases_x8_0         = *(u64*)unencoded_bases;
    u64 bases_x8_1         = *(u64*)(unencoded_bases + 8);
    u64 bases_x8_2         = *(u64*)(unencoded_bases + 16);
    u64 bases_x8_3         = *(u64*)(unencoded_bases + 24);
    u64 bits_2_1_encoded_0 = bases_x8_0 ^ (bases_x8_0 >> 1);
    u64 bits_2_1_encoded_1 = bases_x8_1 ^ (bases_x8_1 >> 1);
    u64 bits_2_1_encoded_2 = bases_x8_2 ^ (bases_x8_2 >> 1);
    u64 bits_2_1_encoded_3 = bases_x8_3 ^ (bases_x8_3 >> 1);
    u64 encoded_bases_x8_0 = _pext_u64(bits_2_1_encoded_0, 0x0606060606060606);
    u64 encoded_bases_x8_1 = _pext_u64(bits_2_1_encoded_1, 0x0606060606060606);
    u64 encoded_bases_x8_2 = _pext_u64(bits_2_1_encoded_2, 0x0606060606060606);
    u64 encoded_bases_x8_3 = _pext_u64(bits_2_1_encoded_3, 0x0606060606060606);

    *(u16*)encoded_bases       = (u16)encoded_bases_x8_0;
    *(u16*)(encoded_bases + 2) = (u16)encoded_bases_x8_1;
    *(u16*)(encoded_bases + 4) = (u16)encoded_bases_x8_2;
    *(u16*)(encoded_bases + 6) = (u16)encoded_bases_x8_3;

    encoded_bases   += 8;
    unencoded_bases += 32;
  }

  while (unencoded_bases < unencoded_bases_x8_end)
  {
    u64 bases_x8_0         = *(u64*)unencoded_bases;
    u64 bits_2_1_encoded_0 = bases_x8_0 ^ (bases_x8_0 >> 1);
    u64 encoded_bases_x8_0 = _pext_u64(bits_2_1_encoded_0, 0x0606060606060606);

    *(u16*)encoded_bases = (u16)encoded_bases_x8_0;
    encoded_bases   += 2;
    unencoded_bases += 8;
  }

  u64 remaining_base_count = base_count & 7;
  encode_bases_scalar(unencoded_bases, encoded_bases, remaining_base_count);
}


void encode_bases_sse4_1(const char* restrict unencoded_bases,
                         u8*         restrict encoded_bases,
                         u64                  base_count)
{
  const char* restrict const unencoded_bases_end = unencoded_bases + (base_count & ~15);

  const __m128i unpacked_encoded_mask = _mm_set1_epi8(0x06);
  const __m128i sum_of_shifts_mul     = _mm_set1_epi32(0x820820);
  const __m128i shuffle_mask          = _mm_set_epi64x(0x8080808080808080, 0x808080800F0B0703);

  // Neither clang nor GCC unroll the loop, even after it has been unrolled manually.
  // MSVC unrolls the loop by default since MSVC v19.27 (Visual Studio 2019 16.7, 5th August 2020)
  #if defined(__clang__)
  #  pragma clang loop unroll_count(4)
  #elif defined(__GNUC__)
  #  pragma GCC unroll 4
  #endif
  while (unencoded_bases < unencoded_bases_end)
  {
    __m128i bases_x16 = _mm_loadu_si128((__m128i*)unencoded_bases);

    // Encode bases in bits 2 (most significant bit) and 1 (least significant bit) of each byte:
    //
    //   bit_21_encoded = bases_x16 ^ (bases_x16 >> 1)
    //                  = xxxxxabx xxxxxcdx xxxxxefx xxxxxghx xxxxxijx xxxxxklx ...
    __m128i base_x16_shr_1  = _mm_srli_epi64(bases_x16, 1);
    __m128i bits_21_encoded = _mm_xor_si128(bases_x16, base_x16_shr_1);

    // Mask out unused bits (bits 7, 6, 5, 4, 3 and 0 of each byte):
    //
    //   unpacked_encoded = bits_21_encoded & 0x0606060606060606
    //                    = 00000ab0 00000cd0 00000ef0 00000gh0 00000ij0 00000kl0 ...
    __m128i unpacked_encoded = _mm_and_si128(bits_21_encoded, unpacked_encoded_mask);

    // Pack encoded bits (2 and 1) of consecutive groups of 4 bytes into 8 bits (4 x 2), in the
    // highest byte of each group of 4 bytes:
    //
    //   packed_8x8_encoded = unpacked_encoded * 0x820820
    //                      = abcdefgh 00000000 00000000 00000000 ijklmnop 00000000 ...
    //
    // The multiplication works as a series of shifts applied to unpacked_encoded, ORed together:
    //
    //   x       = unpacked_encoded
    //   x       = 00000ab0 00000cd0 00000ef0 00000gh0 00000ij0 00000kl0 ...
    //
    //   x << 5  = ab000000 cd000000 ef000000 gh000000 ij000000 kl000000 ...
    //   x << 11 = 00cd0000 00ef0000 00gh0000 00ij0000 00kl0000 00mn0000 ...
    //   x << 17 = 0000ef00 0000gh00 0000ij00 0000kl00 0000mn00 0000op00 ...
    //   x << 23 = 000000gh 000000ij 000000kl 000000mn 000000op 00000000 ...
    //   ORed    = abcdefgh cdefghij efghijkl ghijklmn ijklmnop klmnop00 ...
    //
    //   packed_8x8_encoded = (x << 5) | (x << 11) | (x << 17) | (x << 23)
    //   packed_8x8_encoded = x * (1 << 5) | x * (1 << 11) | x * (1 << 17) | x * (1 << 23)
    //   packed_8x8_encoded = x * ((1 << 5) | (1 << 11) | (1 << 17) | (1 << 23))
    //   packed_8x8_encoded = x * 0x820820
    //
    //   packed_8x8_encoded = unpacked_encoded * 0x820820
    //                      = abcdefgh xxxxxxxx xxxxxxxx xxxxxxxx ijklmnop xxxxxxxx ...
    __m128i packed_8x4_encoded = _mm_mullo_epi32(unpacked_encoded, sum_of_shifts_mul);

    // Further pack the high 8-bit of each 4 bytes (packed encoded bits) into the low 32 bits. Put
    // all other bits to 0 (shuffle mask byte 0x80):
    //
    //   packed_encoded = abcdefgh ijkl.... ........ ....wxyz ABCDEFGH IJKL.... ........ ....WXYZ
    __m128i packed_encoded = _mm_shuffle_epi8(packed_8x4_encoded, shuffle_mask);

    // Finally, extract the packed low 32 bits off the vector register (with the floating-point
    // version of the mov instruction since its throughput is better)
    f32 packed = _mm_cvtss_f32(_mm_castsi128_ps(packed_encoded));
    *(u32*)encoded_bases = *(u32*)&packed;
    
    encoded_bases   += 4;
    unencoded_bases += 16;
  }

  // Some x86 architectures, such as AMD's Bulldozer, support SSE 4.1 but not BMI2 (nor BMI1)
  u64 remaining_base_count = base_count & 15;
  encode_bases_scalar(unencoded_bases, encoded_bases, remaining_base_count);
}


void encode_bases_avx2(const char* restrict unencoded_bases,
                       u8*         restrict encoded_bases,
                       u64                  base_count)
{
  // Load 32 chars per 256 bits, encode and pack them into 64 bits
  const char* const unencoded_bases_end = unencoded_bases + (base_count & ~31);

  const __m256i unpacked_encoded_mask = _mm256_set1_epi8(0x06);
  const __m256i sum_of_shifts_mul     = _mm256_set1_epi32(0x820820);
  const __m256i shuffle_mask          = _mm256_set_epi64x(0x8080808080808080,
                                                          0x0F0B070380808080,
                                                          0x8080808080808080,
                                                          0x808080800F0B0703);
  const __m256i packed_encoded_mask   = _mm256_set_epi64x(0x0, 0x0, 0x0, (5ull << 32));

  while (unencoded_bases < unencoded_bases_end)
  {
    __m256i bases_x32 = _mm256_loadu_si256((__m256i*)unencoded_bases);

    // Encode bases in bits 2 (most significant bit) and 1 (least significant bit) of each byte:
    //
    //   bit_21_encoded = bases_x32 ^ (bases_x32 >> 1)
    //                  = xxxxxabx xxxxxcdx xxxxxefx xxxxxghx xxxxxijx xxxxxklx ...
    __m256i base_x32_shr_1  = _mm256_srli_epi64(bases_x32, 1);
    __m256i bits_21_encoded = _mm256_xor_si256(bases_x32, base_x32_shr_1);

    // Mask out unused bits (bits 7, 6, 5, 4, 3 and 0 of each byte):
    //
    //   unpacked_encoded = bits_21_encoded & 0x0606060606060606
    //                    = 00000ab0 00000cd0 00000ef0 00000gh0 00000ij0 00000kl0 ...
    __m256i unpacked_encoded = _mm256_and_si256(bits_21_encoded, unpacked_encoded_mask);

    // Pack encoded bits (2 and 1) of consecutive groups of 4 bytes into 8 bits (4 x 2), in the
    // highest byte of each group of 4 bytes:
    //
    //   packed_8x8_encoded = unpacked_encoded * 0x820820
    //                      = abcdefgh 00000000 00000000 00000000 ijklmnop 00000000 ...
    //
    // The multiplication works as a series of shifts applied to unpacked_encoded, ORed together:
    //
    //   x       = unpacked_encoded
    //   x       = 00000ab0 00000cd0 00000ef0 00000gh0 00000ij0 00000kl0 ...
    //
    //   x << 5  = ab000000 cd000000 ef000000 gh000000 ij000000 kl000000 ...
    //   x << 11 = 00cd0000 00ef0000 00gh0000 00ij0000 00kl0000 00mn0000 ...
    //   x << 17 = 0000ef00 0000gh00 0000ij00 0000kl00 0000mn00 0000op00 ...
    //   x << 23 = 000000gh 000000ij 000000kl 000000mn 000000op 00000000 ...
    //   ORed    = abcdefgh cdefghij efghijkl ghijklmn ijklmnop klmnop00 ...
    //
    //   packed_8x8_encoded = (x << 5) | (x << 11) | (x << 17) | (x << 23)
    //   packed_8x8_encoded = x * (1 << 5) | x * (1 << 11) | x * (1 << 17) | x * (1 << 23)
    //   packed_8x8_encoded = x * ((1 << 5) | (1 << 11) | (1 << 17) | (1 << 23))
    //   packed_8x8_encoded = x * 0x820820
    //
    //   packed_8x8_encoded = unpacked_encoded * 0x820820
    //                      = abcdefgh xxxxxxxx xxxxxxxx xxxxxxxx ijklmnop xxxxxxxx ...
    __m256i packed_8x8_encoded = _mm256_mullo_epi32(unpacked_encoded, sum_of_shifts_mul);

    // In the two groups of 128 bits ("lanes"), separately, further pack the high 8-bit of each
    // 4 bytes (packed encoded bits) into 32 bits.
    // 
    // In the low  128 bits: put the packed 32 bits in bits 0  to 31 (4 lowest bytes) of the lane
    // In the high 128 bits: put the packed 32 bits in bits 32 to 63 of the lane
    //
    // (Shuffles on 8-bit elements cannot move bits from a lane to another, so the high lane is setup
    // to be later ORed with the low lane)
    //
    // Put all other bits to 0 (shuffle mask byte 0x80) in each lane:
    //
    //   packed_2x32_encoded = abcdefgh ijkl.... ........ ....wxyz 00000000 00000000 00000000 00000000 ...
    //                         00000000 00000000 00000000 00000000 ABCDEFGH IJKL.... ........ ....WXYZ
    __m256i packed_2x32_encoded = _mm256_shuffle_epi8(packed_8x8_encoded, shuffle_mask);

    // Move the 32 packed bits located in the high lane in bits 32 to 63 of the low lane, while
    // keeping the 32 packed bits located in the low lane in bits 0 to 31
    //
    //   packed_encoded = ... abcdefgh ijkl.... ........ ....wxyz ABCDEFGH IJKL.... ........ ....WXYZ

    // Shuffling the high and low lane the same way then using _mm_unpacklo_epi32() may yield a
    // slightly better throughput (~1 GB/s) on some processors
    //
    //   __m128i hi_packed_1x32_encoded = _mm256_extracti128_si256(packed_2x32_encoded, 1);
    //   __m128i unpacked = _mm_unpacklo_epi32(_mm256_castsi256_si128(packed_2x32_encoded), hi_packed_1x32_encoded);
    //   f64 packed = _mm_cvtsd_f64(_mm_castsi128_pd(unpacked));
    __m256i packed_encoded = _mm256_permutevar8x32_epi32(packed_2x32_encoded, packed_encoded_mask);

    // Finally, extract the packed low 64 bits off the vector register (with the floating-point
    // version of the vmov instruction since there's not equivalent for u64)
    f64 packed = _mm256_cvtsd_f64(_mm256_castsi256_pd(packed_encoded));
    *(u64*)encoded_bases = *(u64*)&packed;

    encoded_bases   += 8;
    unencoded_bases += 32;
  }

  u64 remaining_base_count = base_count & 31;
  encode_bases_bmi2(unencoded_bases, encoded_bases, remaining_base_count);
}




////////////////////////////////////////////////////////////////////////////////////////////////////
//// Decoding with specific reference bases (2-bit elements to 8-bit ASCII chars)
char decode_base(u64 encoded_base, const char decoded_table[256 * 4])
{
  u64 idx = encoded_base << 2;
  return decoded_table[idx];
}


// Only the low 8 bits of encoded_bases_x4 are expected to be set, at most
void decode_bases_x4(u64 encoded_bases_x4, char decoded_bases[4], const char decoded_table[256 * 4])
{
  // decoded_table is a table of all 256 combinations of 4 unencoded bases
  *(u32*)decoded_bases = ((u32*)decoded_table)[encoded_bases_x4];
}


void decode_bases_scalar(const u8* restrict encoded_bases,
                         char*     restrict decoded_bases,
                         u64                base_count,
                         const char         decoded_table[256 * 4])
{
  const char* const decoded_bases_end     = decoded_bases + base_count;
  const char* const decoded_bases_x4_end  = decoded_bases + (base_count & ~3);
  const char* const decoded_bases_x64_end = decoded_bases + (base_count & ~63);

  while (decoded_bases < decoded_bases_x64_end)
  {
    u32 bases_x16_0 = *(u32*)encoded_bases;
    u32 bases_x16_1 = *(u32*)(encoded_bases + 4);
    u32 bases_x16_2 = *(u32*)(encoded_bases + 8);
    u32 bases_x16_3 = *(u32*)(encoded_bases + 12);

    u8 bases_x4_0  = (u8)bases_x16_0;
    u8 bases_x4_1  = (u8)(bases_x16_0 >>  8);
    u8 bases_x4_2  = (u8)(bases_x16_0 >> 16);
    u8 bases_x4_3  = (u8)(bases_x16_0 >> 24);
    u8 bases_x4_4  = (u8)bases_x16_1;
    u8 bases_x4_5  = (u8)(bases_x16_1 >>  8);
    u8 bases_x4_6  = (u8)(bases_x16_1 >> 16);
    u8 bases_x4_7  = (u8)(bases_x16_1 >> 24);
    u8 bases_x4_8  = (u8)bases_x16_2;
    u8 bases_x4_9  = (u8)(bases_x16_2 >>  8);
    u8 bases_x4_10 = (u8)(bases_x16_2 >> 16);
    u8 bases_x4_11 = (u8)(bases_x16_2 >> 24);
    u8 bases_x4_12 = (u8)bases_x16_3;
    u8 bases_x4_13 = (u8)(bases_x16_3 >>  8);
    u8 bases_x4_14 = (u8)(bases_x16_3 >> 16);
    u8 bases_x4_15 = (u8)(bases_x16_3 >> 24);

    decode_bases_x4(bases_x4_0,  decoded_bases,      decoded_table);
    decode_bases_x4(bases_x4_1,  decoded_bases +  4, decoded_table);
    decode_bases_x4(bases_x4_2,  decoded_bases +  8, decoded_table);
    decode_bases_x4(bases_x4_3,  decoded_bases + 12, decoded_table);
    decode_bases_x4(bases_x4_4,  decoded_bases + 16, decoded_table);
    decode_bases_x4(bases_x4_5,  decoded_bases + 20, decoded_table);
    decode_bases_x4(bases_x4_6,  decoded_bases + 24, decoded_table);
    decode_bases_x4(bases_x4_7,  decoded_bases + 28, decoded_table);
    decode_bases_x4(bases_x4_8,  decoded_bases + 32, decoded_table);
    decode_bases_x4(bases_x4_9,  decoded_bases + 36, decoded_table);
    decode_bases_x4(bases_x4_10, decoded_bases + 40, decoded_table);
    decode_bases_x4(bases_x4_11, decoded_bases + 44, decoded_table);
    decode_bases_x4(bases_x4_12, decoded_bases + 48, decoded_table);
    decode_bases_x4(bases_x4_13, decoded_bases + 52, decoded_table);
    decode_bases_x4(bases_x4_14, decoded_bases + 56, decoded_table);
    decode_bases_x4(bases_x4_15, decoded_bases + 60, decoded_table);

    encoded_bases += 16;
    decoded_bases += 64;
  }

  while (decoded_bases < decoded_bases_x4_end)
  {
    u8 bases_x4 = *encoded_bases;
    decode_bases_x4(bases_x4, decoded_bases, decoded_table);

    encoded_bases += 1;
    decoded_bases += 4;
  }
  
  if (decoded_bases < decoded_bases_end)
  {
    // At most 3 bases remain
    u64 encoded_bases_xn = *encoded_bases;
    u8  shift            = 0;
    do
    {
      u64 encoded_base = (encoded_bases_xn >> shift) & 0b11;
      *decoded_bases = decode_base(encoded_base, decoded_table);
      decoded_bases += 1;
      shift         += 2;
    }
    while (decoded_bases < decoded_bases_end);
  }
}


void decode_bases_sse4_1(const u8* restrict encoded_bases,
                         char*     restrict decoded_bases,
                         u64                base_count,
                         const char         decoded_table[256 * 4])
{
  const char* const decoded_bases_x32_end = decoded_bases + (base_count & ~15);

  // The 229th element in decoded_table represents the four ordered decoded bases
  // (e.g ACGT, ACGU, acgt, acgu)
  const __m128i spaced_bases_4x4_mask = _mm_set_epi64x(0x8003800380028002, 0x8001800180008000);
  const __m128i reference_bases       = _mm_set1_epi32(*(u32*)(decoded_table + (228 * 4)));
  const __m128i sum_of_shifts_mul     = _mm_set1_epi32(0x40100401);
  const __m128i spaced_bases_x16_mask = _mm_set_epi64x(0x0C0D0E0F08090A0B, 0x0405060700010203);
  const __m128i unset_8bit_msb_mask   = _mm_set1_epi8(0x7F);

  while (decoded_bases < decoded_bases_x32_end)
  {
    // Load 16 packed encoded bases 4 times in the same vector register
    //
    //  bases_x16 = abcdefgh ijklmnop ........ ........ ........ ........ ........ stuvwxyz 
    __m128i bases_x16 = _mm_set1_epi32(*(u32*)encoded_bases);

    // Pad 8-bit sets of 4 encoded bases with 0s so there is only 1 such set per 32 bits (4 bytes)
    //
    //   spaced_bases_4x4 = 00000000 00000000 00000000 abcdefgh 00000000 00000000 00000000 ijklmnop 00..
    __m128i spaced_bases_4x4 = _mm_shuffle_epi8(bases_x16, spaced_bases_4x4_mask);

    // Move each individual encoded base in the high 2 bits of its own 8-bit element:
    //
    //  unordered_spaced_bases_x32_msb = spaced_bases_4x8 * 0x40100401
    //                                   ghxxxxxx efxxxxxx cdxxxxxx abxxxxxx opxxxxxx mnxxxxxx ...
    //
    // The multiplication works as a series of shifts on 2-byte elements in spaced_bases_4x8,
    // added together:
    //
    //   x            = spaced_bases_4x8
    //
    //   x            = 00000000 abcdefgh   00000000 abcdefgh   00000000 ijklmnop   00000000 ijklmnop ...
    //   x_odd        = 00000000 abcdefgh   ________ ________   00000000 ijklmnop   ________ ________ ...
    //   x_even       = ________ ________   00000000 abcdefgh   ________ ________   00000000 ijklmnop ...
    //
    //   x_odd  << 14 = gh000000 00000000   ________ ________   op000000 00000000   ________ ________ ...
    //   x_odd  <<  4 = 0000abcd efgh0000   ________ ________   0000ijkl mnop0000   ________ ________ ...
    //   x_even << 10 = ________ ________   cdefgh00 00000000   ________ ________   klmnop00 00000000 ...
    //   x_even <<  0 = ________ ________   00000000 abcdefgh   ________ ________   00000000 ijklmnop ...
    //   Added        = gh00abcd efgh0000   cdefgh00 abcdefgh   op00ijkl mnop0000   klmnop00 ijklmnop ...
    //
    //   unordered_spaced_bases_x32_msb = (x_even << 0) + (x_even << 10) + (x_odd << 4) + (x_odd << 14)
    //                                  = x_even * ((1 << 0) + (1 << 10)) + x_odd * ((1 << 4) + (1 << 14))
    //                                  = _mm256_mullo_epi16(spaced_bases_4x8, sum_of_shifts_mul)
    //
    // unordered_spaced_bases_x32_lsb = unordered_spaced_bases_x32_msb >> 6
    //                                = 000000gh 00abcdef gh0000cd efgh00ab  000000op 00ijklmn ...
    //
    // ordered_spaced_bases_x32_lsb = efgh00ab gh00abcd 00abcdef 000000gh mnop00ij op00ijkl ...
    __m128i unordered_spaced_bases_x16_msb = _mm_mullo_epi16(spaced_bases_4x4, sum_of_shifts_mul);
    __m128i unordered_spaced_bases_x16_lsb = _mm_srli_epi64(unordered_spaced_bases_x16_msb, 6);
    __m128i ordered_spaced_bases_x16_lsb   = _mm_shuffle_epi8(unordered_spaced_bases_x16_lsb,
                                                              spaced_bases_x16_mask);

    // Mask out the top bits of each 8-bit element, so 0s aren't introduced in the next shuffle:
    //
    //   shuffle_mask = 0xxxxxab 0xxxxxcd 0xxxxxef 0xxxxxgh 0xxxxxij 0xxxxxkl ...
    __m128i shuffle_mask = _mm_and_si128(ordered_spaced_bases_x16_lsb, unset_8bit_msb_mask);

    // shuffle_mask now contains 8-bit elements whose low 4 bits are used as indices on
    // reference_bases's content, which holds decoded bases. Decoded bases in reference_bases are
    // stored as ACGT ACGT ACGT (or U)... Two values with different upper 2 bits but identical
    // 2 bits are indices toward an identical decoded base:
    //
    //   decoded = (01000001 for A, 01000011 for C, 01000111 for G, 01010100 for T) x 32
    __m128i decoded = _mm_shuffle_epi8(reference_bases, shuffle_mask);

    _mm_storeu_si128((__m128i*)decoded_bases, decoded);

    encoded_bases += 4;
    decoded_bases += 16;
  }

  // Process the remaining bases with scalar operations
  u64 remaining_base_count = base_count & 15;
  decode_bases_scalar(encoded_bases, decoded_bases, remaining_base_count, decoded_table);
}


void decode_bases_avx2(const u8* restrict encoded_bases,
                       char*     restrict decoded_bases,
                       u64                base_count,
                       const char         decoded_table[256 * 4])
{
  const __m256i spaced_bases_4x8_mask = _mm256_set_epi64x(0x8007800780068006,
                                                          0x8005800580048004,
                                                          0x8003800380028002,
                                                          0x8001800180008000);
  const __m256i sum_of_shifts_mul     = _mm256_set1_epi32(0x40100401);
  const __m256i spaced_bases_x32_mask = _mm256_set_epi64x(0x0C0D0E0F08090A0B,
                                                          0x0405060700010203,
                                                          0x0C0D0E0F08090A0B,
                                                          0x0405060700010203);
  const __m256i unset_8bit_msb_mask   = _mm256_set1_epi8(0x7F);
  const __m256i reference_bases       = _mm256_set1_epi32(*(u32*)(decoded_table + (228 * 4)));

  const char* const decoded_bases_x32_end = decoded_bases + (base_count & ~31);
  while (decoded_bases < decoded_bases_x32_end)
  {
    // Load 32 packed encoded bases 4 times in the same vector register
    //
    //  bases_x32 = abcdefgh ijklmnop ........ ........ ........ ........ ........ stuvwxyz
    __m256i bases_x32 = _mm256_set1_epi64x(*(u64*)encoded_bases);

    // Pad 8-bit sets of 4 encoded bases with 0s so there is only 1 such set per 32 bits (4 bytes)
    //
    //   spaced_bases_4x8 = 00000000 abcdefgh 00000000 abcdefgh 00000000 ijklmnop 00000000 ijklmnop 00..
    __m256i spaced_bases_4x8 = _mm256_shuffle_epi8(bases_x32, spaced_bases_4x8_mask);

    // Move each individual encoded base in the high 2 bits of its own 8-bit element:
    //
    //  unordered_spaced_bases_x32_msb = spaced_bases_4x8 * 0x40100401
    //                                   ghxxxxxx efxxxxxx cdxxxxxx abxxxxxx opxxxxxx mnxxxxxx ...
    //
    // The multiplication works as a series of shifts on 2-byte elements in spaced_bases_4x8,
    // added together:
    //
    //   x            = spaced_bases_4x8
    //
    //   x            = 00000000 abcdefgh   00000000 abcdefgh   00000000 ijklmnop   00000000 ijklmnop ...
    //   x_odd        = 00000000 abcdefgh   ________ ________   00000000 ijklmnop   ________ ________ ...
    //   x_even       = ________ ________   00000000 abcdefgh   ________ ________   00000000 ijklmnop ...
    //
    //   x_odd  << 14 = gh000000 00000000   ________ ________   op000000 00000000   ________ ________ ...
    //   x_odd  <<  4 = 0000abcd efgh0000   ________ ________   0000ijkl mnop0000   ________ ________ ...
    //   x_even << 10 = ________ ________   cdefgh00 00000000   ________ ________   klmnop00 00000000 ...
    //   x_even <<  0 = ________ ________   00000000 abcdefgh   ________ ________   00000000 ijklmnop ...
    //   Added        = gh00abcd efgh0000   cdefgh00 abcdefgh   op00ijkl mnop0000   klmnop00 ijklmnop ...
    //
    //   unordered_spaced_bases_x32_msb = (x_even << 0) + (x_even << 10) + (x_odd << 4) + (x_odd << 14)
    //                                  = x_even * ((1 << 0) + (1 << 10)) + x_odd * ((1 << 4) + (1 << 14))
    //                                  = _mm256_mullo_epi16(spaced_bases_4x8, sum_of_shifts_mul)
    //
    // unordered_spaced_bases_x32_lsb = unordered_spaced_bases_x32_msb >> 6
    //                                = 000000gh 00abcdef gh0000cd efgh00ab  000000op 00ijklmn ...
    //
    // ordered_spaced_bases_x32_lsb = efgh00ab gh00abcd 00abcdef 000000gh mnop00ij op00ijkl ...
    __m256i unordered_spaced_bases_x32_msb = _mm256_mullo_epi16(spaced_bases_4x8, sum_of_shifts_mul);
    __m256i unordered_spaced_bases_x32_lsb = _mm256_srli_epi64(unordered_spaced_bases_x32_msb, 6);
    __m256i ordered_spaced_bases_x32_lsb   = _mm256_shuffle_epi8(unordered_spaced_bases_x32_lsb,
                                                                 spaced_bases_x32_mask);

    // Mask out the top bits of each 8-bit element, so 0s aren't introduced in the next shuffle:
    //
    //   shuffle_mask = 0xxxxxab 0xxxxxcd 0xxxxxef 0xxxxxgh 0xxxxxij 0xxxxxkl ...
    __m256i shuffle_mask = _mm256_and_si256(ordered_spaced_bases_x32_lsb, unset_8bit_msb_mask);

    // shuffle_mask now contains 8-bit elements whose low 4 bits are used as indices on
    // reference_bases's content, which holds decoded bases. Decoded bases in reference_bases are
    // stored as ACGT ACGT ACGT (or U)... Two values with different upper 2 bits but identical
    // 2 bits represent an index toward identical decoded base:
    //
    //   decoded = (01000001 for A, 01000011 for C, 01000111 for G, 01010100 for T) x 32
    __m256i decoded = _mm256_shuffle_epi8(reference_bases, shuffle_mask);

    _mm256_storeu_si256((__m256i*)decoded_bases, decoded);

    encoded_bases += 8;
    decoded_bases += 32;
  }

  u64 remaining_base_count = base_count & 31;
  decode_bases_scalar(encoded_bases, decoded_bases, remaining_base_count, decoded_table);
}