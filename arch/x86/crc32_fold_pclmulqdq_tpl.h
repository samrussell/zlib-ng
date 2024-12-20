/*
 * Compute the CRC32 using a parallelized folding approach with the PCLMULQDQ
 * instruction.
 *
 * A white paper describing this algorithm can be found at:
 *     doc/crc-pclmulqdq.pdf
 *
 * Copyright (C) 2013 Intel Corporation. All rights reserved.
 * Copyright (C) 2016 Marian Beermann (support for initial value)
 * Authors:
 *     Wajdi Feghali   <wajdi.k.feghali@intel.com>
 *     Jim Guilford    <james.guilford@intel.com>
 *     Vinodh Gopal    <vinodh.gopal@intel.com>
 *     Erdinc Ozturk   <erdinc.ozturk@intel.com>
 *     Jim Kukunas     <james.t.kukunas@linux.intel.com>
 *
 * For conditions of distribution and use, see copyright notice in zlib.h
 */

#ifdef COPY
Z_INTERNAL void CRC32_FOLD_COPY(crc32_fold *crc, uint8_t *dst, const uint8_t *src, size_t len) {
#else
Z_INTERNAL void CRC32_FOLD(crc32_fold *crc, const uint8_t *src, size_t len, uint32_t init_crc) {
#endif
    unsigned long algn_diff;
    __m128i xmm_t0, xmm_t1, xmm_t2, xmm_t3;
    __m128i xmm_crc0, xmm_crc1, xmm_crc2, xmm_crc3;
    __m128i xmm_crc_part = _mm_setzero_si128();
    char ALIGNED_(16) partial_buf[16] = { 0 };
#ifndef COPY
    __m128i xmm_initial = _mm_cvtsi32_si128(init_crc);
    int32_t first = init_crc != 0;

    /* The CRC functions don't call this for input < 16, as a minimum of 16 bytes of input is needed
     * for the aligning load that occurs.  If there's an initial CRC, to carry it forward through
     * the folded CRC there must be 16 - src % 16 + 16 bytes available, which by definition can be
     * up to 15 bytes + one full vector load. */
    assert(len >= 16 || first == 0);
#endif
    crc32_fold_load((__m128i *)crc->fold, &xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3);

    if (len < 16) {
        if (len == 0)
            return;

        memcpy(partial_buf, src, len);
        xmm_crc_part = _mm_load_si128((const __m128i *)partial_buf);
#ifdef COPY
        memcpy(dst, partial_buf, len);
#endif
        goto partial;
    }

    algn_diff = ((uintptr_t)16 - ((uintptr_t)src & 0xF)) & 0xF;
    if (algn_diff) {
        xmm_crc_part = _mm_loadu_si128((__m128i *)src);
#ifdef COPY
        _mm_storeu_si128((__m128i *)dst, xmm_crc_part);
        dst += algn_diff;
#else
        XOR_INITIAL128(xmm_crc_part);

        if (algn_diff < 4 && init_crc != 0) {
            xmm_t0 = xmm_crc_part;
            if (len >= 32) {
                xmm_crc_part = _mm_loadu_si128((__m128i*)src + 1);
                fold_1(&xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3);
                xmm_crc3 = _mm_xor_si128(xmm_crc3, xmm_t0);
            } else {
                memcpy(partial_buf, src + 16, len - 16);
                xmm_crc_part = _mm_load_si128((__m128i*)partial_buf);
                fold_1(&xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3);
                xmm_crc3 = _mm_xor_si128(xmm_crc3, xmm_t0);
                src += 16;
                len -= 16;
#ifdef COPY
                dst -= algn_diff;
#endif
                goto partial;
            }

            src += 16;
            len -= 16;
        }
#endif

        partial_fold(algn_diff, &xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3, &xmm_crc_part);

        src += algn_diff;
        len -= algn_diff;
    }

#ifdef X86_VPCLMULQDQ
    if (len >= 256) {
#ifdef COPY
        size_t n = fold_16_vpclmulqdq_copy(&xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3, dst, src, len);
        dst += n;
#else
        size_t n = fold_16_vpclmulqdq(&xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3, src, len,
            xmm_initial, first);
        first = 0;
#endif
        len -= n;
        src += n;
    }
#endif

        while (len >= 512 + 64 + 16*4) {
                __m128i shift544_shift480 = _mm_set_epi64x(0x1D9513D7, 0x8F352D95);
                // __m128i shift672_shift608 = _mm_set_epi64x(0xAE0B5394, 0x1C279815);
                // __m128i shift800_shift736 = _mm_set_epi64x(0x57C54819, 0xDF068DC2);
                __m128i shift1056_shift992 = _mm_set_epi64x(0x910EEEC1, 0x33FFF533);
            __m128i bonus4 = _mm_loadu_si128(src + (16 * 4));
            __m128i bonus3 = _mm_loadu_si128(src + (16 * 5));
            __m128i bonus2 = _mm_loadu_si128(src + (16 * 6));
            __m128i bonus1 = _mm_loadu_si128(src + (16 * 7));
            src += 16*4;
            len -= 16*4;
            xmm_t0 = _mm_loadu_si128(src + (16 * 4));
            xmm_t1 = _mm_loadu_si128(src + (16 * 5));
            xmm_t2 = _mm_loadu_si128(src + (16 * 6)) ^ bonus4;
            xmm_t3 = _mm_loadu_si128(src + (16 * 7)) ^ bonus3;

            // now we fold xmm_crc0 onto xmm_crc1
            // fold 8x because we stole a value
            __m128i fold_high1 = _mm_clmulepi64_si128(xmm_crc0, shift1056_shift992, 0x11);
            __m128i fold_low1 = _mm_clmulepi64_si128(xmm_crc0, shift1056_shift992, 0x00);
            xmm_crc0 = _mm_xor_si128(xmm_t0, fold_high1);
            xmm_crc0 = _mm_xor_si128(xmm_crc0, fold_low1);
            __m128i fold_high2 = _mm_clmulepi64_si128(xmm_crc1, shift1056_shift992, 0x11);
            __m128i fold_low2 = _mm_clmulepi64_si128(xmm_crc1, shift1056_shift992, 0x00);
            xmm_crc1 = _mm_xor_si128(xmm_t1, fold_high2);
            xmm_crc1 = _mm_xor_si128(xmm_crc1, fold_low2);
            __m128i fold_high3 = _mm_clmulepi64_si128(xmm_crc2, shift1056_shift992, 0x11);
            __m128i fold_low3 = _mm_clmulepi64_si128(xmm_crc2, shift1056_shift992, 0x00);
            xmm_crc2 = _mm_xor_si128(xmm_t2, fold_high3);
            xmm_crc2 = _mm_xor_si128(xmm_crc2, fold_low3);
            __m128i fold_high4 = _mm_clmulepi64_si128(xmm_crc3, shift1056_shift992, 0x11);
            __m128i fold_low4 = _mm_clmulepi64_si128(xmm_crc3, shift1056_shift992, 0x00);
            xmm_crc3 = _mm_xor_si128(xmm_t3, fold_high4);
            xmm_crc3 = _mm_xor_si128(xmm_crc3, fold_low4);
            xmm_t0 = _mm_loadu_si128(src + (16 * 8)) ^ bonus2;
            xmm_t1 = _mm_loadu_si128(src + (16 * 9)) ^ bonus1 ^ bonus4;
            xmm_t2 = _mm_loadu_si128(src + (16 * 10)) ^ bonus3 ^ bonus4;
            xmm_t3 = _mm_loadu_si128(src + (16 * 11)) ^ bonus2 ^ bonus3;

            // now we fold xmm_crc0 onto xmm_crc1
            fold_high1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x11);
            fold_low1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x00);
            xmm_crc0 = _mm_xor_si128(xmm_t0, fold_high1);
            xmm_crc0 = _mm_xor_si128(xmm_crc0, fold_low1);
            fold_high2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x11);
            fold_low2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x00);
            xmm_crc1 = _mm_xor_si128(xmm_t1, fold_high2);
            xmm_crc1 = _mm_xor_si128(xmm_crc1, fold_low2);
            fold_high3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x11);
            fold_low3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x00);
            xmm_crc2 = _mm_xor_si128(xmm_t2, fold_high3);
            xmm_crc2 = _mm_xor_si128(xmm_crc2, fold_low3);
            fold_high4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x11);
            fold_low4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x00);
            xmm_crc3 = _mm_xor_si128(xmm_t3, fold_high4);
            xmm_crc3 = _mm_xor_si128(xmm_crc3, fold_low4);
            xmm_t0 = _mm_loadu_si128(src + (16 * 12)) ^ bonus1 ^ bonus2;
            xmm_t1 = _mm_loadu_si128(src + (16 * 13)) ^ bonus1;
            xmm_t2 = _mm_loadu_si128(src + (16 * 14));
            xmm_t3 = _mm_loadu_si128(src + (16 * 15));

            // now we fold xmm_crc0 onto xmm_crc1
            fold_high1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x11);
            fold_low1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x00);
            xmm_crc0 = _mm_xor_si128(xmm_t0, fold_high1);
            xmm_crc0 = _mm_xor_si128(xmm_crc0, fold_low1);
            fold_high2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x11);
            fold_low2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x00);
            xmm_crc1 = _mm_xor_si128(xmm_t1, fold_high2);
            xmm_crc1 = _mm_xor_si128(xmm_crc1, fold_low2);
            fold_high3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x11);
            fold_low3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x00);
            xmm_crc2 = _mm_xor_si128(xmm_t2, fold_high3);
            xmm_crc2 = _mm_xor_si128(xmm_crc2, fold_low3);
            fold_high4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x11);
            fold_low4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x00);
            xmm_crc3 = _mm_xor_si128(xmm_t3, fold_high4);
            xmm_crc3 = _mm_xor_si128(xmm_crc3, fold_low4);
            xmm_t0 = _mm_loadu_si128(src + (16 * 16)) ^ bonus4;
            xmm_t1 = _mm_loadu_si128(src + (16 * 17)) ^ bonus3;
            xmm_t2 = _mm_loadu_si128(src + (16 * 18)) ^ bonus2;
            xmm_t3 = _mm_loadu_si128(src + (16 * 19)) ^ bonus1;

            // now we fold xmm_crc0 onto xmm_crc1
            fold_high1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x11);
            fold_low1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x00);
            xmm_crc0 = _mm_xor_si128(xmm_t0, fold_high1);
            xmm_crc0 = _mm_xor_si128(xmm_crc0, fold_low1);
            fold_high2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x11);
            fold_low2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x00);
            xmm_crc1 = _mm_xor_si128(xmm_t1, fold_high2);
            xmm_crc1 = _mm_xor_si128(xmm_crc1, fold_low2);
            fold_high3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x11);
            fold_low3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x00);
            xmm_crc2 = _mm_xor_si128(xmm_t2, fold_high3);
            xmm_crc2 = _mm_xor_si128(xmm_crc2, fold_low3);
            fold_high4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x11);
            fold_low4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x00);
            xmm_crc3 = _mm_xor_si128(xmm_t3, fold_high4);
            xmm_crc3 = _mm_xor_si128(xmm_crc3, fold_low4);
            xmm_t0 = _mm_loadu_si128(src + (16 * 20)) ^ bonus4;
            xmm_t1 = _mm_loadu_si128(src + (16 * 21)) ^ bonus3 ^ bonus4;
            xmm_t2 = _mm_loadu_si128(src + (16 * 22)) ^ bonus2 ^ bonus3 ^ bonus4;
            xmm_t3 = _mm_loadu_si128(src + (16 * 23)) ^ bonus1 ^ bonus2 ^ bonus3;

            // now we fold xmm_crc0 onto xmm_crc1
            fold_high1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x11);
            fold_low1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x00);
            xmm_crc0 = _mm_xor_si128(xmm_t0, fold_high1);
            xmm_crc0 = _mm_xor_si128(xmm_crc0, fold_low1);
            fold_high2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x11);
            fold_low2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x00);
            xmm_crc1 = _mm_xor_si128(xmm_t1, fold_high2);
            xmm_crc1 = _mm_xor_si128(xmm_crc1, fold_low2);
            fold_high3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x11);
            fold_low3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x00);
            xmm_crc2 = _mm_xor_si128(xmm_t2, fold_high3);
            xmm_crc2 = _mm_xor_si128(xmm_crc2, fold_low3);
            fold_high4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x11);
            fold_low4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x00);
            xmm_crc3 = _mm_xor_si128(xmm_t3, fold_high4);
            xmm_crc3 = _mm_xor_si128(xmm_crc3, fold_low4);
            xmm_t0 = _mm_loadu_si128(src + (16 * 24)) ^ bonus1 ^ bonus2 ^ bonus4;
            xmm_t1 = _mm_loadu_si128(src + (16 * 25)) ^ bonus1 ^ bonus3 ^ bonus4;
            xmm_t2 = _mm_loadu_si128(src + (16 * 26)) ^ bonus2 ^ bonus3;
            xmm_t3 = _mm_loadu_si128(src + (16 * 27)) ^ bonus1 ^ bonus2 ^ bonus4;

            // now we fold xmm_crc0 onto xmm_crc1
            fold_high1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x11);
            fold_low1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x00);
            xmm_crc0 = _mm_xor_si128(xmm_t0, fold_high1);
            xmm_crc0 = _mm_xor_si128(xmm_crc0, fold_low1);
            fold_high2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x11);
            fold_low2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x00);
            xmm_crc1 = _mm_xor_si128(xmm_t1, fold_high2);
            xmm_crc1 = _mm_xor_si128(xmm_crc1, fold_low2);
            fold_high3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x11);
            fold_low3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x00);
            xmm_crc2 = _mm_xor_si128(xmm_t2, fold_high3);
            xmm_crc2 = _mm_xor_si128(xmm_crc2, fold_low3);
            fold_high4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x11);
            fold_low4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x00);
            xmm_crc3 = _mm_xor_si128(xmm_t3, fold_high4);
            xmm_crc3 = _mm_xor_si128(xmm_crc3, fold_low4);
            xmm_t0 = _mm_loadu_si128(src + (16 * 28)) ^ bonus1 ^ bonus3 ^ bonus4;
            xmm_t1 = _mm_loadu_si128(src + (16 * 29)) ^ bonus2 ^ bonus3;
            xmm_t2 = _mm_loadu_si128(src + (16 * 30)) ^ bonus1 ^ bonus2 ^ bonus4;
            xmm_t3 = _mm_loadu_si128(src + (16 * 31)) ^ bonus1 ^ bonus3 ^ bonus4;

            // now we fold xmm_crc0 onto xmm_crc1
            fold_high1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x11);
            fold_low1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x00);
            xmm_crc0 = _mm_xor_si128(xmm_t0, fold_high1);
            xmm_crc0 = _mm_xor_si128(xmm_crc0, fold_low1);
            fold_high2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x11);
            fold_low2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x00);
            xmm_crc1 = _mm_xor_si128(xmm_t1, fold_high2);
            xmm_crc1 = _mm_xor_si128(xmm_crc1, fold_low2);
            fold_high3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x11);
            fold_low3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x00);
            xmm_crc2 = _mm_xor_si128(xmm_t2, fold_high3);
            xmm_crc2 = _mm_xor_si128(xmm_crc2, fold_low3);
            fold_high4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x11);
            fold_low4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x00);
            xmm_crc3 = _mm_xor_si128(xmm_t3, fold_high4);
            xmm_crc3 = _mm_xor_si128(xmm_crc3, fold_low4);
            xmm_t0 = _mm_loadu_si128(src + (16 * 32)) ^ bonus2 ^ bonus3 ^ bonus4;
            xmm_t1 = _mm_loadu_si128(src + (16 * 33)) ^ bonus1 ^ bonus2 ^ bonus3;
            xmm_t2 = _mm_loadu_si128(src + (16 * 34)) ^ bonus1 ^ bonus2;
            xmm_t3 = _mm_loadu_si128(src + (16 * 35)) ^ bonus1;

            // now we fold xmm_crc0 onto xmm_crc1
            fold_high1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x11);
            fold_low1 = _mm_clmulepi64_si128(xmm_crc0, shift544_shift480, 0x00);
            xmm_crc0 = _mm_xor_si128(xmm_t0, fold_high1);
            xmm_crc0 = _mm_xor_si128(xmm_crc0, fold_low1);
            fold_high2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x11);
            fold_low2 = _mm_clmulepi64_si128(xmm_crc1, shift544_shift480, 0x00);
            xmm_crc1 = _mm_xor_si128(xmm_t1, fold_high2);
            xmm_crc1 = _mm_xor_si128(xmm_crc1, fold_low2);
            fold_high3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x11);
            fold_low3 = _mm_clmulepi64_si128(xmm_crc2, shift544_shift480, 0x00);
            xmm_crc2 = _mm_xor_si128(xmm_t2, fold_high3);
            xmm_crc2 = _mm_xor_si128(xmm_crc2, fold_low3);
            fold_high4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x11);
            fold_low4 = _mm_clmulepi64_si128(xmm_crc3, shift544_shift480, 0x00);
            xmm_crc3 = _mm_xor_si128(xmm_t3, fold_high4);
            xmm_crc3 = _mm_xor_si128(xmm_crc3, fold_low4);

            len -= 512;
            src += 512;
        }

    while (len >= 64) {
        len -= 64;
        xmm_t0 = _mm_load_si128((__m128i *)src);
        xmm_t1 = _mm_load_si128((__m128i *)src + 1);
        xmm_t2 = _mm_load_si128((__m128i *)src + 2);
        xmm_t3 = _mm_load_si128((__m128i *)src + 3);
        src += 64;

        fold_4(&xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3);
#ifdef COPY
        _mm_storeu_si128((__m128i *)dst, xmm_t0);
        _mm_storeu_si128((__m128i *)dst + 1, xmm_t1);
        _mm_storeu_si128((__m128i *)dst + 2, xmm_t2);
        _mm_storeu_si128((__m128i *)dst + 3, xmm_t3);
        dst += 64;
#else
        XOR_INITIAL128(xmm_t0);
#endif

        xmm_crc0 = _mm_xor_si128(xmm_crc0, xmm_t0);
        xmm_crc1 = _mm_xor_si128(xmm_crc1, xmm_t1);
        xmm_crc2 = _mm_xor_si128(xmm_crc2, xmm_t2);
        xmm_crc3 = _mm_xor_si128(xmm_crc3, xmm_t3);
    }

    /*
     * len = num bytes left - 64
     */
    if (len >= 48) {
        len -= 48;

        xmm_t0 = _mm_load_si128((__m128i *)src);
        xmm_t1 = _mm_load_si128((__m128i *)src + 1);
        xmm_t2 = _mm_load_si128((__m128i *)src + 2);
        src += 48;
#ifdef COPY
        _mm_storeu_si128((__m128i *)dst, xmm_t0);
        _mm_storeu_si128((__m128i *)dst + 1, xmm_t1);
        _mm_storeu_si128((__m128i *)dst + 2, xmm_t2);
        dst += 48;
#else
        XOR_INITIAL128(xmm_t0);
#endif
        fold_3(&xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3);

        xmm_crc1 = _mm_xor_si128(xmm_crc1, xmm_t0);
        xmm_crc2 = _mm_xor_si128(xmm_crc2, xmm_t1);
        xmm_crc3 = _mm_xor_si128(xmm_crc3, xmm_t2);
    } else if (len >= 32) {
        len -= 32;

        xmm_t0 = _mm_load_si128((__m128i *)src);
        xmm_t1 = _mm_load_si128((__m128i *)src + 1);
        src += 32;
#ifdef COPY
        _mm_storeu_si128((__m128i *)dst, xmm_t0);
        _mm_storeu_si128((__m128i *)dst + 1, xmm_t1);
        dst += 32;
#else
        XOR_INITIAL128(xmm_t0);
#endif
        fold_2(&xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3);

        xmm_crc2 = _mm_xor_si128(xmm_crc2, xmm_t0);
        xmm_crc3 = _mm_xor_si128(xmm_crc3, xmm_t1);
    } else if (len >= 16) {
        len -= 16;
        xmm_t0 = _mm_load_si128((__m128i *)src);
        src += 16;
#ifdef COPY
        _mm_storeu_si128((__m128i *)dst, xmm_t0);
        dst += 16;
#else
        XOR_INITIAL128(xmm_t0);
#endif
        fold_1(&xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3);

        xmm_crc3 = _mm_xor_si128(xmm_crc3, xmm_t0);
    }

partial:
    if (len) {
        memcpy(&xmm_crc_part, src, len);
#ifdef COPY
        _mm_storeu_si128((__m128i *)partial_buf, xmm_crc_part);
        memcpy(dst, partial_buf, len);
#endif
        partial_fold(len, &xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3, &xmm_crc_part);
    }

    crc32_fold_save((__m128i *)crc->fold, &xmm_crc0, &xmm_crc1, &xmm_crc2, &xmm_crc3);
}
