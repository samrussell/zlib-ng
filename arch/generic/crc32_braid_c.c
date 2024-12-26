/* crc32_braid.c -- compute the CRC-32 of a data stream
 * Copyright (C) 1995-2022 Mark Adler
 * For conditions of distribution and use, see copyright notice in zlib.h
 *
 * This interleaved implementation of a CRC makes use of pipelined multiple
 * arithmetic-logic units, commonly found in modern CPU cores. It is due to
 * Kadatch and Jenkins (2010). See doc/crc-doc.1.0.pdf in this distribution.
 */

#include "zbuild.h"
#include "crc32_braid_p.h"
#include "crc32_braid_tbl.h"

/*
  A CRC of a message is computed on N braids of words in the message, where
  each word consists of W bytes (4 or 8). If N is 3, for example, then three
  running sparse CRCs are calculated respectively on each braid, at these
  indices in the array of words: 0, 3, 6, ..., 1, 4, 7, ..., and 2, 5, 8, ...
  This is done starting at a word boundary, and continues until as many blocks
  of N * W bytes as are available have been processed. The results are combined
  into a single CRC at the end. For this code, N must be in the range 1..6 and
  W must be 4 or 8. The upper limit on N can be increased if desired by adding
  more #if blocks, extending the patterns apparent in the code. In addition,
  crc32 tables would need to be regenerated, if the maximum N value is increased.

  N and W are chosen empirically by benchmarking the execution time on a given
  processor. The choices for N and W below were based on testing on Intel Kaby
  Lake i7, AMD Ryzen 7, ARM Cortex-A57, Sparc64-VII, PowerPC POWER9, and MIPS64
  Octeon II processors. The Intel, AMD, and ARM processors were all fastest
  with N=5, W=8. The Sparc, PowerPC, and MIPS64 were all fastest at N=5, W=4.
  They were all tested with either gcc or clang, all using the -O3 optimization
  level. Your mileage may vary.
*/

/* ========================================================================= */
#ifdef W
/*
  Return the CRC of the W bytes in the word_t data, taking the
  least-significant byte of the word as the first byte of data, without any pre
  or post conditioning. This is used to combine the CRCs of each braid.
 */
#if BYTE_ORDER == LITTLE_ENDIAN
static uint32_t crc_word(z_word_t data) {
    int k;
    for (k = 0; k < W; k++)
        data = (data >> 8) ^ crc_table[data & 0xff];
    return (uint32_t)data;
}
#elif BYTE_ORDER == BIG_ENDIAN
static z_word_t crc_word(z_word_t data) {
    int k;
    for (k = 0; k < W; k++)
        data = (data << 8) ^
            crc_big_table[(data >> ((W - 1) << 3)) & 0xff];
    return data;
}
#endif /* BYTE_ORDER */

#endif /* W */

/* ========================================================================= */
Z_INTERNAL uint32_t crc32_braid_base(uint32_t c, const uint8_t *buf, size_t len) {

#ifdef W
    /* If provided enough bytes, do a braided CRC calculation. */
    if (len >= N * W + W - 1) {
        size_t blks;
        z_word_t const *words;
        int k;

        /* Compute the CRC up to a z_word_t boundary. */
        while (len && ((uintptr_t)buf & (W - 1)) != 0) {
            len--;
            DO1;
        }

        /* Compute the CRC on as many N z_word_t blocks as are available. */
        blks = len / (N * W);
        len -= blks * N * W;
        words = (z_word_t const *)buf;

        z_word_t crc0, word0, comb;
#if N > 1
        z_word_t crc1, word1;
#if N > 2
        z_word_t crc2, word2;
#if N > 3
        z_word_t crc3, word3;
#if N > 4
        z_word_t crc4, word4;
#if N > 5
        z_word_t crc5, word5;
#endif
#endif
#endif
#endif
#endif
        /* Initialize the CRC for each braid. */
        crc0 = ZSWAPWORD(c);
#if N > 1
        crc1 = 0;
#if N > 2
        crc2 = 0;
#if N > 3
        crc3 = 0;
#if N > 4
        crc4 = 0;
#if N > 5
        crc5 = 0;
#endif
#endif
#endif
#endif
#endif
        /* Process the first blks-1 blocks, computing the CRCs on each braid independently. */
        while (--blks) {
            /* Load the word for each braid into registers. */
            word0 = crc0 ^ words[0];
#if N > 1
            word1 = crc1 ^ words[1];
#if N > 2
            word2 = crc2 ^ words[2];
#if N > 3
            word3 = crc3 ^ words[3];
#if N > 4
            word4 = crc4 ^ words[4];
#if N > 5
            word5 = crc5 ^ words[5];
#endif
#endif
#endif
#endif
#endif
            words += N;

            /* Compute and update the CRC for each word. The loop should get unrolled. */
            crc0 = BRAID_TABLE[0][word0 & 0xff];
#if N > 1
            crc1 = BRAID_TABLE[0][word1 & 0xff];
#if N > 2
            crc2 = BRAID_TABLE[0][word2 & 0xff];
#if N > 3
            crc3 = BRAID_TABLE[0][word3 & 0xff];
#if N > 4
            crc4 = BRAID_TABLE[0][word4 & 0xff];
#if N > 5
            crc5 = BRAID_TABLE[0][word5 & 0xff];
#endif
#endif
#endif
#endif
#endif
            for (k = 1; k < W; k++) {
                crc0 ^= BRAID_TABLE[k][(word0 >> (k << 3)) & 0xff];
#if N > 1
                crc1 ^= BRAID_TABLE[k][(word1 >> (k << 3)) & 0xff];
#if N > 2
                crc2 ^= BRAID_TABLE[k][(word2 >> (k << 3)) & 0xff];
#if N > 3
                crc3 ^= BRAID_TABLE[k][(word3 >> (k << 3)) & 0xff];
#if N > 4
                crc4 ^= BRAID_TABLE[k][(word4 >> (k << 3)) & 0xff];
#if N > 5
                crc5 ^= BRAID_TABLE[k][(word5 >> (k << 3)) & 0xff];
#endif
#endif
#endif
#endif
#endif
            }
        }

        /* Process the last block, combining the CRCs of the N braids at the same time. */
        comb = crc_word(crc0 ^ words[0]);
#if N > 1
        comb = crc_word(crc1 ^ words[1] ^ comb);
#if N > 2
        comb = crc_word(crc2 ^ words[2] ^ comb);
#if N > 3
        comb = crc_word(crc3 ^ words[3] ^ comb);
#if N > 4
        comb = crc_word(crc4 ^ words[4] ^ comb);
#if N > 5
        comb = crc_word(crc5 ^ words[5] ^ comb);
#endif
#endif
#endif
#endif
#endif
        words += N;
        Assert(comb <= UINT32_MAX, "comb should fit in uint32_t");
        c = (uint32_t)ZSWAPWORD(comb);

        /* Update the pointer to the remaining bytes to process. */
        buf = (const unsigned char *)words;
    }

#endif /* W */

    /* Complete the computation of the CRC on any remaining bytes. */
    while (len >= 8) {
        len -= 8;
        DO8;
    }
    while (len) {
        len--;
        DO1;
    }

    /* Return the CRC, post-conditioned. */
    return c;
}

/* Implement Chorba algorithm from https://arxiv.org/abs/2412.16398 */
#define bitbuffersizebytes (16 * 1024 * sizeof(z_word_t))
#define bitbuffersizezwords (bitbuffersizebytes / sizeof(z_word_t))
#define bitbuffersizeqwords (bitbuffersizebytes / sizeof(uint64_t))

uint32_t chorba_118960_nondestructive (uint32_t crc, const uint8_t* input, size_t len) {
    z_word_t bitbuffer[bitbuffersizezwords];
    const uint8_t* bitbufferbytes = (const uint8_t*) bitbuffer;

    int i = 0;

    z_word_t next1 = crc;
    z_word_t next2 = 0;
    z_word_t next3 = 0;
    z_word_t next4 = 0;
    z_word_t next5 = 0;
    z_word_t next6 = 0;
    z_word_t next7 = 0;
    z_word_t next8 = 0;
    z_word_t next9 = 0;
    z_word_t next10 = 0;
    z_word_t next11 = 0;
    z_word_t next12 = 0;
    z_word_t next13 = 0;
    z_word_t next14 = 0;
    z_word_t next15 = 0;
    z_word_t next16 = 0;
    z_word_t next17 = 0;
    z_word_t next18 = 0;
    z_word_t next19 = 0;
    z_word_t next20 = 0;
    z_word_t next21 = 0;
    z_word_t next22 = 0;
    crc = 0;

    // do a first pass to zero out bitbuffer
    for(; i < (14848 * sizeof(z_word_t)); i += (32 * sizeof(z_word_t))) {
        z_word_t in1, in2, in3, in4, in5, in6, in7, in8;
        z_word_t in9, in10, in11, in12, in13, in14, in15, in16;
        z_word_t in17, in18, in19, in20, in21, in22, in23, in24;
        z_word_t in25, in26, in27, in28, in29, in30, in31, in32;
        int outoffset1 = ((i / sizeof(z_word_t)) + 14848) % bitbuffersizezwords;
        int outoffset2 = ((i / sizeof(z_word_t)) + 14880) % bitbuffersizezwords;

        in1 = *((z_word_t*) (input + i + (0 * sizeof(z_word_t)))) ^ next1;
        in2 = *((z_word_t*) (input + i + (1 * sizeof(z_word_t)))) ^ next2;
        in3 = *((z_word_t*) (input + i + (2 * sizeof(z_word_t)))) ^ next3;
        in4 = *((z_word_t*) (input + i + (3 * sizeof(z_word_t)))) ^ next4;
        in5 = *((z_word_t*) (input + i + (4 * sizeof(z_word_t)))) ^ next5;
        in6 = *((z_word_t*) (input + i + (5 * sizeof(z_word_t)))) ^ next6;
        in7 = *((z_word_t*) (input + i + (6 * sizeof(z_word_t)))) ^ next7;
        in8 = *((z_word_t*) (input + i + (7 * sizeof(z_word_t)))) ^ next8 ^ in1;
        in9 = *((z_word_t*) (input + i + (8 * sizeof(z_word_t)))) ^ next9 ^ in2;
        in10 = *((z_word_t*) (input + i + (9 * sizeof(z_word_t)))) ^ next10 ^ in3;
        in11 = *((z_word_t*) (input + i + (10 * sizeof(z_word_t)))) ^ next11 ^ in4;
        in12 = *((z_word_t*) (input + i + (11 * sizeof(z_word_t)))) ^ next12 ^ in1 ^ in5;
        in13 = *((z_word_t*) (input + i + (12 * sizeof(z_word_t)))) ^ next13 ^ in2 ^ in6;
        in14 = *((z_word_t*) (input + i + (13 * sizeof(z_word_t)))) ^ next14 ^ in3 ^ in7;
        in15 = *((z_word_t*) (input + i + (14 * sizeof(z_word_t)))) ^ next15 ^ in4 ^ in8;
        in16 = *((z_word_t*) (input + i + (15 * sizeof(z_word_t)))) ^ next16 ^ in5 ^ in9;
        in17 = *((z_word_t*) (input + i + (16 * sizeof(z_word_t)))) ^ next17 ^ in6 ^ in10;
        in18 = *((z_word_t*) (input + i + (17 * sizeof(z_word_t)))) ^ next18 ^ in7 ^ in11;
        in19 = *((z_word_t*) (input + i + (18 * sizeof(z_word_t)))) ^ next19 ^ in8 ^ in12;
        in20 = *((z_word_t*) (input + i + (19 * sizeof(z_word_t)))) ^ next20 ^ in9 ^ in13;
        in21 = *((z_word_t*) (input + i + (20 * sizeof(z_word_t)))) ^ next21 ^ in10 ^ in14;
        in22 = *((z_word_t*) (input + i + (21 * sizeof(z_word_t)))) ^ next22 ^ in11 ^ in15;
        in23 = *((z_word_t*) (input + i + (22 * sizeof(z_word_t)))) ^ in1 ^ in12 ^ in16;
        in24 = *((z_word_t*) (input + i + (23 * sizeof(z_word_t)))) ^ in2 ^ in13 ^ in17;
        in25 = *((z_word_t*) (input + i + (24 * sizeof(z_word_t)))) ^ in3 ^ in14 ^ in18;
        in26 = *((z_word_t*) (input + i + (25 * sizeof(z_word_t)))) ^ in4 ^ in15 ^ in19;
        in27 = *((z_word_t*) (input + i + (26 * sizeof(z_word_t)))) ^ in5 ^ in16 ^ in20;
        in28 = *((z_word_t*) (input + i + (27 * sizeof(z_word_t)))) ^ in6 ^ in17 ^ in21;
        in29 = *((z_word_t*) (input + i + (28 * sizeof(z_word_t)))) ^ in7 ^ in18 ^ in22;
        in30 = *((z_word_t*) (input + i + (29 * sizeof(z_word_t)))) ^ in8 ^ in19 ^ in23;
        in31 = *((z_word_t*) (input + i + (30 * sizeof(z_word_t)))) ^ in9 ^ in20 ^ in24;
        in32 = *((z_word_t*) (input + i + (31 * sizeof(z_word_t)))) ^ in10 ^ in21 ^ in25;

        next1 = in11 ^ in22 ^ in26;
        next2 = in12 ^ in23 ^ in27;
        next3 = in13 ^ in24 ^ in28;
        next4 = in14 ^ in25 ^ in29;
        next5 = in15 ^ in26 ^ in30;
        next6 = in16 ^ in27 ^ in31;
        next7 = in17 ^ in28 ^ in32;
        next8 = in18 ^ in29;
        next9 = in19 ^ in30;
        next10 = in20 ^ in31;
        next11 = in21 ^ in32;
        next12 = in22;
        next13 = in23;
        next14 = in24;
        next15 = in25;
        next16 = in26;
        next17 = in27;
        next18 = in28;
        next19 = in29;
        next20 = in30;
        next21 = in31;
        next22 = in32;

        bitbuffer[outoffset1 + 22] = in1;
        bitbuffer[outoffset1 + 23] = in2;
        bitbuffer[outoffset1 + 24] = in3;
        bitbuffer[outoffset1 + 25] = in4;
        bitbuffer[outoffset1 + 26] = in5;
        bitbuffer[outoffset1 + 27] = in6;
        bitbuffer[outoffset1 + 28] = in7;
        bitbuffer[outoffset1 + 29] = in8;
        bitbuffer[outoffset1 + 30] = in9;
        bitbuffer[outoffset1 + 31] = in10;
        bitbuffer[outoffset2 + 0] = in11;
        bitbuffer[outoffset2 + 1] = in12;
        bitbuffer[outoffset2 + 2] = in13;
        bitbuffer[outoffset2 + 3] = in14;
        bitbuffer[outoffset2 + 4] = in15;
        bitbuffer[outoffset2 + 5] = in16;
        bitbuffer[outoffset2 + 6] = in17;
        bitbuffer[outoffset2 + 7] = in18;
        bitbuffer[outoffset2 + 8] = in19;
        bitbuffer[outoffset2 + 9] = in20;
        bitbuffer[outoffset2 + 10] = in21;
        bitbuffer[outoffset2 + 11] = in22;
        bitbuffer[outoffset2 + 12] = in23;
        bitbuffer[outoffset2 + 13] = in24;
        bitbuffer[outoffset2 + 14] = in25;
        bitbuffer[outoffset2 + 15] = in26;
        bitbuffer[outoffset2 + 16] = in27;
        bitbuffer[outoffset2 + 17] = in28;
        bitbuffer[outoffset2 + 18] = in29;
        bitbuffer[outoffset2 + 19] = in30;
        bitbuffer[outoffset2 + 20] = in31;
        bitbuffer[outoffset2 + 21] = in32;
    }

    // one intermediate pass where we pull half the values
    for(; i < (14880 * sizeof(z_word_t)); i += (32 * sizeof(z_word_t))) {
        z_word_t in1, in2, in3, in4, in5, in6, in7, in8;
        z_word_t in9, in10, in11, in12, in13, in14, in15, in16;
        z_word_t in17, in18, in19, in20, in21, in22, in23, in24;
        z_word_t in25, in26, in27, in28, in29, in30, in31, in32;
        int inoffset = (i / sizeof(z_word_t)) % bitbuffersizezwords;
        int outoffset1 = ((i / sizeof(z_word_t)) + 14848) % bitbuffersizezwords;
        int outoffset2 = ((i / sizeof(z_word_t)) + 14880) % bitbuffersizezwords;

        in1 = *((z_word_t*) (input + i + (0 * sizeof(z_word_t)))) ^ next1;
        in2 = *((z_word_t*) (input + i + (1 * sizeof(z_word_t)))) ^ next2;
        in3 = *((z_word_t*) (input + i + (2 * sizeof(z_word_t)))) ^ next3;
        in4 = *((z_word_t*) (input + i + (3 * sizeof(z_word_t)))) ^ next4;
        in5 = *((z_word_t*) (input + i + (4 * sizeof(z_word_t)))) ^ next5;
        in6 = *((z_word_t*) (input + i + (5 * sizeof(z_word_t)))) ^ next6;
        in7 = *((z_word_t*) (input + i + (6 * sizeof(z_word_t)))) ^ next7;
        in8 = *((z_word_t*) (input + i + (7 * sizeof(z_word_t)))) ^ next8 ^ in1;
        in9 = *((z_word_t*) (input + i + (8 * sizeof(z_word_t)))) ^ next9 ^ in2;
        in10 = *((z_word_t*) (input + i + (9 * sizeof(z_word_t)))) ^ next10 ^ in3;
        in11 = *((z_word_t*) (input + i + (10 * sizeof(z_word_t)))) ^ next11 ^ in4;
        in12 = *((z_word_t*) (input + i + (11 * sizeof(z_word_t)))) ^ next12 ^ in1 ^ in5;
        in13 = *((z_word_t*) (input + i + (12 * sizeof(z_word_t)))) ^ next13 ^ in2 ^ in6;
        in14 = *((z_word_t*) (input + i + (13 * sizeof(z_word_t)))) ^ next14 ^ in3 ^ in7;
        in15 = *((z_word_t*) (input + i + (14 * sizeof(z_word_t)))) ^ next15 ^ in4 ^ in8;
        in16 = *((z_word_t*) (input + i + (15 * sizeof(z_word_t)))) ^ next16 ^ in5 ^ in9;
        in17 = *((z_word_t*) (input + i + (16 * sizeof(z_word_t)))) ^ next17 ^ in6 ^ in10;
        in18 = *((z_word_t*) (input + i + (17 * sizeof(z_word_t)))) ^ next18 ^ in7 ^ in11;
        in19 = *((z_word_t*) (input + i + (18 * sizeof(z_word_t)))) ^ next19 ^ in8 ^ in12;
        in20 = *((z_word_t*) (input + i + (19 * sizeof(z_word_t)))) ^ next20 ^ in9 ^ in13;
        in21 = *((z_word_t*) (input + i + (20 * sizeof(z_word_t)))) ^ next21 ^ in10 ^ in14;
        in22 = *((z_word_t*) (input + i + (21 * sizeof(z_word_t)))) ^ next22 ^ in11 ^ in15;
        in23 = *((z_word_t*) (input + i + (22 * sizeof(z_word_t)))) ^ in1 ^ in12 ^ in16 ^ bitbuffer[inoffset + 22];
        in24 = *((z_word_t*) (input + i + (23 * sizeof(z_word_t)))) ^ in2 ^ in13 ^ in17 ^ bitbuffer[inoffset + 23];
        in25 = *((z_word_t*) (input + i + (24 * sizeof(z_word_t)))) ^ in3 ^ in14 ^ in18 ^ bitbuffer[inoffset + 24];
        in26 = *((z_word_t*) (input + i + (25 * sizeof(z_word_t)))) ^ in4 ^ in15 ^ in19 ^ bitbuffer[inoffset + 25];
        in27 = *((z_word_t*) (input + i + (26 * sizeof(z_word_t)))) ^ in5 ^ in16 ^ in20 ^ bitbuffer[inoffset + 26];
        in28 = *((z_word_t*) (input + i + (27 * sizeof(z_word_t)))) ^ in6 ^ in17 ^ in21 ^ bitbuffer[inoffset + 27];
        in29 = *((z_word_t*) (input + i + (28 * sizeof(z_word_t)))) ^ in7 ^ in18 ^ in22 ^ bitbuffer[inoffset + 28];
        in30 = *((z_word_t*) (input + i + (29 * sizeof(z_word_t)))) ^ in8 ^ in19 ^ in23 ^ bitbuffer[inoffset + 29];
        in31 = *((z_word_t*) (input + i + (30 * sizeof(z_word_t)))) ^ in9 ^ in20 ^ in24 ^ bitbuffer[inoffset + 30];
        in32 = *((z_word_t*) (input + i + (31 * sizeof(z_word_t)))) ^ in10 ^ in21 ^ in25 ^ bitbuffer[inoffset + 31];

        next1 = in11 ^ in22 ^ in26;
        next2 = in12 ^ in23 ^ in27;
        next3 = in13 ^ in24 ^ in28;
        next4 = in14 ^ in25 ^ in29;
        next5 = in15 ^ in26 ^ in30;
        next6 = in16 ^ in27 ^ in31;
        next7 = in17 ^ in28 ^ in32;
        next8 = in18 ^ in29;
        next9 = in19 ^ in30;
        next10 = in20 ^ in31;
        next11 = in21 ^ in32;
        next12 = in22;
        next13 = in23;
        next14 = in24;
        next15 = in25;
        next16 = in26;
        next17 = in27;
        next18 = in28;
        next19 = in29;
        next20 = in30;
        next21 = in31;
        next22 = in32;

        bitbuffer[outoffset1 + 22] = in1;
        bitbuffer[outoffset1 + 23] = in2;
        bitbuffer[outoffset1 + 24] = in3;
        bitbuffer[outoffset1 + 25] = in4;
        bitbuffer[outoffset1 + 26] = in5;
        bitbuffer[outoffset1 + 27] = in6;
        bitbuffer[outoffset1 + 28] = in7;
        bitbuffer[outoffset1 + 29] = in8;
        bitbuffer[outoffset1 + 30] = in9;
        bitbuffer[outoffset1 + 31] = in10;
        bitbuffer[outoffset2 + 0] = in11;
        bitbuffer[outoffset2 + 1] = in12;
        bitbuffer[outoffset2 + 2] = in13;
        bitbuffer[outoffset2 + 3] = in14;
        bitbuffer[outoffset2 + 4] = in15;
        bitbuffer[outoffset2 + 5] = in16;
        bitbuffer[outoffset2 + 6] = in17;
        bitbuffer[outoffset2 + 7] = in18;
        bitbuffer[outoffset2 + 8] = in19;
        bitbuffer[outoffset2 + 9] = in20;
        bitbuffer[outoffset2 + 10] = in21;
        bitbuffer[outoffset2 + 11] = in22;
        bitbuffer[outoffset2 + 12] = in23;
        bitbuffer[outoffset2 + 13] = in24;
        bitbuffer[outoffset2 + 14] = in25;
        bitbuffer[outoffset2 + 15] = in26;
        bitbuffer[outoffset2 + 16] = in27;
        bitbuffer[outoffset2 + 17] = in28;
        bitbuffer[outoffset2 + 18] = in29;
        bitbuffer[outoffset2 + 19] = in30;
        bitbuffer[outoffset2 + 20] = in31;
        bitbuffer[outoffset2 + 21] = in32;
    }

    for(; i + (14870 + 64) * sizeof(z_word_t) < len; i += (32 * sizeof(z_word_t))) {
        z_word_t in1, in2, in3, in4, in5, in6, in7, in8;
        z_word_t in9, in10, in11, in12, in13, in14, in15, in16;
        z_word_t in17, in18, in19, in20, in21, in22, in23, in24;
        z_word_t in25, in26, in27, in28, in29, in30, in31, in32;
        int inoffset = (i / sizeof(z_word_t)) % bitbuffersizezwords;
        int outoffset1 = ((i / sizeof(z_word_t)) + 14848) % bitbuffersizezwords;
        int outoffset2 = ((i / sizeof(z_word_t)) + 14880) % bitbuffersizezwords;

        in1 = *((z_word_t*) (input + i + (0 * sizeof(z_word_t)))) ^ next1 ^ bitbuffer[inoffset + 0];
        in2 = *((z_word_t*) (input + i + (1 * sizeof(z_word_t)))) ^ next2 ^ bitbuffer[inoffset + 1];
        in3 = *((z_word_t*) (input + i + (2 * sizeof(z_word_t)))) ^ next3 ^ bitbuffer[inoffset + 2];
        in4 = *((z_word_t*) (input + i + (3 * sizeof(z_word_t)))) ^ next4 ^ bitbuffer[inoffset + 3];
        in5 = *((z_word_t*) (input + i + (4 * sizeof(z_word_t)))) ^ next5 ^ bitbuffer[inoffset + 4];
        in6 = *((z_word_t*) (input + i + (5 * sizeof(z_word_t)))) ^ next6 ^ bitbuffer[inoffset + 5];
        in7 = *((z_word_t*) (input + i + (6 * sizeof(z_word_t)))) ^ next7 ^ bitbuffer[inoffset + 6];
        in8 = *((z_word_t*) (input + i + (7 * sizeof(z_word_t)))) ^ next8 ^ in1 ^ bitbuffer[inoffset + 7];
        in9 = *((z_word_t*) (input + i + (8 * sizeof(z_word_t)))) ^ next9 ^ in2 ^ bitbuffer[inoffset + 8];
        in10 = *((z_word_t*) (input + i + (9 * sizeof(z_word_t)))) ^ next10 ^ in3 ^ bitbuffer[inoffset + 9];
        in11 = *((z_word_t*) (input + i + (10 * sizeof(z_word_t)))) ^ next11 ^ in4 ^ bitbuffer[inoffset + 10];
        in12 = *((z_word_t*) (input + i + (11 * sizeof(z_word_t)))) ^ next12 ^ in1 ^ in5 ^ bitbuffer[inoffset + 11];
        in13 = *((z_word_t*) (input + i + (12 * sizeof(z_word_t)))) ^ next13 ^ in2 ^ in6 ^ bitbuffer[inoffset + 12];
        in14 = *((z_word_t*) (input + i + (13 * sizeof(z_word_t)))) ^ next14 ^ in3 ^ in7 ^ bitbuffer[inoffset + 13];
        in15 = *((z_word_t*) (input + i + (14 * sizeof(z_word_t)))) ^ next15 ^ in4 ^ in8 ^ bitbuffer[inoffset + 14];
        in16 = *((z_word_t*) (input + i + (15 * sizeof(z_word_t)))) ^ next16 ^ in5 ^ in9 ^ bitbuffer[inoffset + 15];
        in17 = *((z_word_t*) (input + i + (16 * sizeof(z_word_t)))) ^ next17 ^ in6 ^ in10 ^ bitbuffer[inoffset + 16];
        in18 = *((z_word_t*) (input + i + (17 * sizeof(z_word_t)))) ^ next18 ^ in7 ^ in11 ^ bitbuffer[inoffset + 17];
        in19 = *((z_word_t*) (input + i + (18 * sizeof(z_word_t)))) ^ next19 ^ in8 ^ in12 ^ bitbuffer[inoffset + 18];
        in20 = *((z_word_t*) (input + i + (19 * sizeof(z_word_t)))) ^ next20 ^ in9 ^ in13 ^ bitbuffer[inoffset + 19];
        in21 = *((z_word_t*) (input + i + (20 * sizeof(z_word_t)))) ^ next21 ^ in10 ^ in14 ^ bitbuffer[inoffset + 20];
        in22 = *((z_word_t*) (input + i + (21 * sizeof(z_word_t)))) ^ next22 ^ in11 ^ in15 ^ bitbuffer[inoffset + 21];
        in23 = *((z_word_t*) (input + i + (22 * sizeof(z_word_t)))) ^ in1 ^ in12 ^ in16 ^ bitbuffer[inoffset + 22];
        in24 = *((z_word_t*) (input + i + (23 * sizeof(z_word_t)))) ^ in2 ^ in13 ^ in17 ^ bitbuffer[inoffset + 23];
        in25 = *((z_word_t*) (input + i + (24 * sizeof(z_word_t)))) ^ in3 ^ in14 ^ in18 ^ bitbuffer[inoffset + 24];
        in26 = *((z_word_t*) (input + i + (25 * sizeof(z_word_t)))) ^ in4 ^ in15 ^ in19 ^ bitbuffer[inoffset + 25];
        in27 = *((z_word_t*) (input + i + (26 * sizeof(z_word_t)))) ^ in5 ^ in16 ^ in20 ^ bitbuffer[inoffset + 26];
        in28 = *((z_word_t*) (input + i + (27 * sizeof(z_word_t)))) ^ in6 ^ in17 ^ in21 ^ bitbuffer[inoffset + 27];
        in29 = *((z_word_t*) (input + i + (28 * sizeof(z_word_t)))) ^ in7 ^ in18 ^ in22 ^ bitbuffer[inoffset + 28];
        in30 = *((z_word_t*) (input + i + (29 * sizeof(z_word_t)))) ^ in8 ^ in19 ^ in23 ^ bitbuffer[inoffset + 29];
        in31 = *((z_word_t*) (input + i + (30 * sizeof(z_word_t)))) ^ in9 ^ in20 ^ in24 ^ bitbuffer[inoffset + 30];
        in32 = *((z_word_t*) (input + i + (31 * sizeof(z_word_t)))) ^ in10 ^ in21 ^ in25 ^ bitbuffer[inoffset + 31];

        next1 = in11 ^ in22 ^ in26;
        next2 = in12 ^ in23 ^ in27;
        next3 = in13 ^ in24 ^ in28;
        next4 = in14 ^ in25 ^ in29;
        next5 = in15 ^ in26 ^ in30;
        next6 = in16 ^ in27 ^ in31;
        next7 = in17 ^ in28 ^ in32;
        next8 = in18 ^ in29;
        next9 = in19 ^ in30;
        next10 = in20 ^ in31;
        next11 = in21 ^ in32;
        next12 = in22;
        next13 = in23;
        next14 = in24;
        next15 = in25;
        next16 = in26;
        next17 = in27;
        next18 = in28;
        next19 = in29;
        next20 = in30;
        next21 = in31;
        next22 = in32;

        bitbuffer[outoffset1 + 22] = in1;
        bitbuffer[outoffset1 + 23] = in2;
        bitbuffer[outoffset1 + 24] = in3;
        bitbuffer[outoffset1 + 25] = in4;
        bitbuffer[outoffset1 + 26] = in5;
        bitbuffer[outoffset1 + 27] = in6;
        bitbuffer[outoffset1 + 28] = in7;
        bitbuffer[outoffset1 + 29] = in8;
        bitbuffer[outoffset1 + 30] = in9;
        bitbuffer[outoffset1 + 31] = in10;
        bitbuffer[outoffset2 + 0] = in11;
        bitbuffer[outoffset2 + 1] = in12;
        bitbuffer[outoffset2 + 2] = in13;
        bitbuffer[outoffset2 + 3] = in14;
        bitbuffer[outoffset2 + 4] = in15;
        bitbuffer[outoffset2 + 5] = in16;
        bitbuffer[outoffset2 + 6] = in17;
        bitbuffer[outoffset2 + 7] = in18;
        bitbuffer[outoffset2 + 8] = in19;
        bitbuffer[outoffset2 + 9] = in20;
        bitbuffer[outoffset2 + 10] = in21;
        bitbuffer[outoffset2 + 11] = in22;
        bitbuffer[outoffset2 + 12] = in23;
        bitbuffer[outoffset2 + 13] = in24;
        bitbuffer[outoffset2 + 14] = in25;
        bitbuffer[outoffset2 + 15] = in26;
        bitbuffer[outoffset2 + 16] = in27;
        bitbuffer[outoffset2 + 17] = in28;
        bitbuffer[outoffset2 + 18] = in29;
        bitbuffer[outoffset2 + 19] = in30;
        bitbuffer[outoffset2 + 20] = in31;
        bitbuffer[outoffset2 + 21] = in32;
    }

    bitbuffer[(i / sizeof(z_word_t) + 0) % bitbuffersizezwords] ^= next1;
    bitbuffer[(i / sizeof(z_word_t) + 1) % bitbuffersizezwords] ^= next2;
    bitbuffer[(i / sizeof(z_word_t) + 2) % bitbuffersizezwords] ^= next3;
    bitbuffer[(i / sizeof(z_word_t) + 3) % bitbuffersizezwords] ^= next4;
    bitbuffer[(i / sizeof(z_word_t) + 4) % bitbuffersizezwords] ^= next5;
    bitbuffer[(i / sizeof(z_word_t) + 5) % bitbuffersizezwords] ^= next6;
    bitbuffer[(i / sizeof(z_word_t) + 6) % bitbuffersizezwords] ^= next7;
    bitbuffer[(i / sizeof(z_word_t) + 7) % bitbuffersizezwords] ^= next8;
    bitbuffer[(i / sizeof(z_word_t) + 8) % bitbuffersizezwords] ^= next9;
    bitbuffer[(i / sizeof(z_word_t) + 9) % bitbuffersizezwords] ^= next10;
    bitbuffer[(i / sizeof(z_word_t) + 10) % bitbuffersizezwords] ^= next11;
    bitbuffer[(i / sizeof(z_word_t) + 11) % bitbuffersizezwords] ^= next12;
    bitbuffer[(i / sizeof(z_word_t) + 12) % bitbuffersizezwords] ^= next13;
    bitbuffer[(i / sizeof(z_word_t) + 13) % bitbuffersizezwords] ^= next14;
    bitbuffer[(i / sizeof(z_word_t) + 14) % bitbuffersizezwords] ^= next15;
    bitbuffer[(i / sizeof(z_word_t) + 15) % bitbuffersizezwords] ^= next16;
    bitbuffer[(i / sizeof(z_word_t) + 16) % bitbuffersizezwords] ^= next17;
    bitbuffer[(i / sizeof(z_word_t) + 17) % bitbuffersizezwords] ^= next18;
    bitbuffer[(i / sizeof(z_word_t) + 18) % bitbuffersizezwords] ^= next19;
    bitbuffer[(i / sizeof(z_word_t) + 19) % bitbuffersizezwords] ^= next20;
    bitbuffer[(i / sizeof(z_word_t) + 20) % bitbuffersizezwords] ^= next21;
    bitbuffer[(i / sizeof(z_word_t) + 21) % bitbuffersizezwords] ^= next22;

    for (int j = 14870; j < 14870 + 60; j++) {
        bitbuffer[(j + (i / sizeof(z_word_t))) % bitbuffersizezwords] = 0;
    }

    next1 = 0;
    next2 = 0;
    next3 = 0;
    next4 = 0;
    next5 = 0;
    uint8_t final[72] = {0};

    for(; i + 72 < len; i += 32) {
        uint64_t in1;
        uint64_t in2;
        uint64_t in3;
        uint64_t in4;
        uint64_t a1, a2, a3, a4;
        uint64_t b1, b2, b3, b4;
        uint64_t c1, c2, c3, c4;
        uint64_t d1, d2, d3, d4;

        uint64_t out1;
        uint64_t out2;
        uint64_t out3;
        uint64_t out4;
        uint64_t out5;

        in1 = *((uint64_t*) (input + i));
        in2 = *((uint64_t*) (input + i + (8*1)));
#if BYTE_ORDER == BIG_ENDIAN
        in1 = ZSWAP64(in1);
        in2 = ZSWAP64(in1);
#endif
        in1 ^= next1 ^ bitbuffer[(i / sizeof(uint64_t)) % bitbuffersizeqwords];
        in2 ^= next2 ^ bitbuffer[(i / sizeof(uint64_t) + 1) % bitbuffersizeqwords];

        a1 = (in1 << 17) ^ (in1 << 55);
        a2 = (in1 >> 47) ^ (in1 >> 9) ^ (in1 << 19);
        a3 = (in1 >> 45) ^ (in1 << 44);
        a4 = (in1 >> 20);
        
        b1 = (in2 << 17) ^ (in2 << 55);
        b2 = (in2 >> 47) ^ (in2 >> 9) ^ (in2 << 19);
        b3 = (in2 >> 45) ^ (in2 << 44);
        b4 = (in2 >> 20);

        in3 = *((uint64_t*) (input + i + (8*2)));
        in4 = *((uint64_t*) (input + i + (8*3)));
#if BYTE_ORDER == BIG_ENDIAN
        in1 = ZSWAP64(in1);
        in2 = ZSWAP64(in1);
#endif
        in3 ^= next3 ^ a1 ^ bitbuffer[(i / sizeof(uint64_t) + 2) % bitbuffersizeqwords];
        in4 ^= next4 ^ a2 ^ b1 ^ bitbuffer[(i / sizeof(uint64_t) + 3) % bitbuffersizeqwords];

        c1 = (in3 << 17) ^ (in3 << 55);
        c2 = (in3 >> 47) ^ (in3 >> 9) ^ (in3 << 19);
        c3 = (in3 >> 45) ^ (in3 << 44);
        c4 = (in3 >> 20);
        
        d1 = (in4 << 17) ^ (in4 << 55);
        d2 = (in4 >> 47) ^ (in4 >> 9) ^ (in4 << 19);
        d3 = (in4 >> 45) ^ (in4 << 44);
        d4 = (in4 >> 20);

        out1 = a3 ^ b2 ^ c1;
        out2 = a4 ^ b3 ^ c2 ^ d1;
        out3 = b4 ^ c3 ^ d2;
        out4 = c4 ^ d3;
        out5 = d4;

        next1 = next5 ^ out1;
        next2 = out2;
        next3 = out3;
        next4 = out4;
        next5 = out5;

    }

#if BYTE_ORDER == BIG_ENDIAN
    next1 = ZSWAP64(next1);
    next2 = ZSWAP64(next2);
    next3 = ZSWAP64(next3);
    next4 = ZSWAP64(next4);
    next5 = ZSWAP64(next5);
#endif

    memcpy(final, input+i, len-i);
    *((uint64_t*) (final + (0*8))) ^= next1;
    *((uint64_t*) (final + (1*8))) ^= next2;
    *((uint64_t*) (final + (2*8))) ^= next3;
    *((uint64_t*) (final + (3*8))) ^= next4;
    *((uint64_t*) (final + (4*8))) ^= next5;

    for(int j = 0; j<len-i; j++) {
        crc = crc_table[(crc ^ final[j] ^ bitbufferbytes[(j+i) % bitbuffersizebytes]) & 0xff] ^ (crc >> 8);
    }

    return crc;
}

/* Implement Chorba algorithm from https://arxiv.org/abs/2412.16398 */
uint32_t chorba_small_nondestructive (uint32_t crc, const uint8_t* buf, size_t len) {
    const uint8_t* input = buf;
    uint8_t final[72] = {0};
    uint64_t next1 = crc;
    crc = 0;
    uint64_t next2 = 0;
    uint64_t next3 = 0;
    uint64_t next4 = 0;
    uint64_t next5 = 0;

    int i = 0;
    for(; i + 72 < len; i += 32) {
        uint64_t in1;
        uint64_t in2;
        uint64_t in3;
        uint64_t in4;
        uint64_t a1, a2, a3, a4;
        uint64_t b1, b2, b3, b4;
        uint64_t c1, c2, c3, c4;
        uint64_t d1, d2, d3, d4;

        uint64_t out1;
        uint64_t out2;
        uint64_t out3;
        uint64_t out4;
        uint64_t out5;

        in1 = *((uint64_t*) (input + i)) ^ next1;
        in2 = *((uint64_t*) (input + i + (8*1))) ^ next2;

        a1 = (in1 << 17) ^ (in1 << 55);
        a2 = (in1 >> 47) ^ (in1 >> 9) ^ (in1 << 19);
        a3 = (in1 >> 45) ^ (in1 << 44);
        a4 = (in1 >> 20);
        
        b1 = (in2 << 17) ^ (in2 << 55);
        b2 = (in2 >> 47) ^ (in2 >> 9) ^ (in2 << 19);
        b3 = (in2 >> 45) ^ (in2 << 44);
        b4 = (in2 >> 20);

        in3 = *((uint64_t*) (input + i + (8*2))) ^ next3 ^ a1;
        in4 = *((uint64_t*) (input + i + (8*3))) ^ next4 ^ a2 ^ b1;

        c1 = (in3 << 17) ^ (in3 << 55);
        c2 = (in3 >> 47) ^ (in3 >> 9) ^ (in3 << 19);
        c3 = (in3 >> 45) ^ (in3 << 44);
        c4 = (in3 >> 20);
        
        d1 = (in4 << 17) ^ (in4 << 55);
        d2 = (in4 >> 47) ^ (in4 >> 9) ^ (in4 << 19);
        d3 = (in4 >> 45) ^ (in4 << 44);
        d4 = (in4 >> 20);

        out1 = a3 ^ b2 ^ c1;
        out2 = a4 ^ b3 ^ c2 ^ d1;
        out3 = b4 ^ c3 ^ d2;
        out4 = c4 ^ d3;
        out5 = d4;

        next1 = next5 ^ out1;
        next2 = out2;
        next3 = out3;
        next4 = out4;
        next5 = out5;

    }

    memcpy(final, input+i, len-i);
    *((uint64_t*) (final + (0*8))) ^= next1;
    *((uint64_t*) (final + (1*8))) ^= next2;
    *((uint64_t*) (final + (2*8))) ^= next3;
    *((uint64_t*) (final + (3*8))) ^= next4;
    *((uint64_t*) (final + (4*8))) ^= next5;

    for(int j = 0; j<len-i; j++) {
        crc = crc_table[(crc ^ final[j]) & 0xff] ^ (crc >> 8);
    }

    return crc;
}

Z_INTERNAL uint32_t PREFIX(crc32_braid)(uint32_t crc, const uint8_t *buf, size_t len) {
    uint32_t c;

    c = (~crc) & 0xffffffff;
    if(len > 512 * 1024)
        c = chorba_118960_nondestructive(c, buf, len);
    else if (len > 72)
        c = chorba_small_nondestructive(c, buf, len);
    else {
        c = crc32_braid_base(c, buf, len);
    }

    /* Return the CRC, post-conditioned. */
    return c ^ 0xffffffff;
}
