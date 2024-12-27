/* crc32_acle.c -- compute the CRC-32 of a data stream
 * Copyright (C) 1995-2006, 2010, 2011, 2012 Mark Adler
 * Copyright (C) 2016 Yang Zhang
 * For conditions of distribution and use, see copyright notice in zlib.h
 *
*/

#ifdef ARM_ACLE
#include "acle_intrins.h"
#include "zbuild.h"

Z_INTERNAL Z_TARGET_CRC uint32_t crc32_acle(uint32_t crc, const uint8_t *buf, size_t len) {
    Z_REGISTER uint32_t c;
    Z_REGISTER const uint16_t *buf2;
    Z_REGISTER const uint32_t *buf4;
    Z_REGISTER const uint64_t *buf8;

    c = ~crc;

    if (UNLIKELY(len == 1)) {
        c = __crc32b(c, *buf);
        c = ~c;
        return c;
    }

    if ((ptrdiff_t)buf & (sizeof(uint64_t) - 1)) {
        if (len && ((ptrdiff_t)buf & 1)) {
            c = __crc32b(c, *buf++);
            len--;
        }

        if ((len >= sizeof(uint16_t)) && ((ptrdiff_t)buf & sizeof(uint16_t))) {
            buf2 = (const uint16_t *) buf;
            c = __crc32h(c, *buf2++);
            len -= sizeof(uint16_t);
            buf4 = (const uint32_t *) buf2;
        } else {
            buf4 = (const uint32_t *) buf;
        }

        if ((len >= sizeof(uint32_t)) && ((ptrdiff_t)buf & sizeof(uint32_t))) {
            c = __crc32w(c, *buf4++);
            len -= sizeof(uint32_t);
        }

        buf8 = (const uint64_t *) buf4;
    } else {
        buf8 = (const uint64_t *) buf;
    }

    while (len >= (sizeof(uint64_t) * (32 + 2))) {
        uint64_t chorba4 = *buf8++ ^ c;
        len -= sizeof(uint64_t);
        uint64_t chorba3 = *buf8++;
        len -= sizeof(uint64_t);
        uint64_t chorba2 = *buf8++;
        len -= sizeof(uint64_t);
        uint64_t chorba1 = *buf8++;
        len -= sizeof(uint64_t);
        c = __crc32d(0, *buf8++); // 1
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++); // 2
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba4); // 3
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba3); // 4
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba2); // 5
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba4); // 6
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba3 ^ chorba4); // 7
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba2 ^ chorba3); // 8
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba2); //9
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1); // 10
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++); // 11
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++); // 12
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba4); // 13
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba3); // 14
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba2); // 15
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1); // 16
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba4); // 17
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba3 ^ chorba4); // 18
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba2 ^ chorba3 ^ chorba4); // 19
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba2 ^ chorba3); // 20
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba2 ^ chorba4); // 21
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba3 ^ chorba4); // 22
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba2 ^ chorba3); // 23
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba2 ^ chorba4); // 24
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba3 ^ chorba4); // 25
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba2 ^ chorba3); // 26
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba2 ^ chorba4); // 27
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba3 ^ chorba4); // 28
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba2 ^ chorba3 ^ chorba4); // 29
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba2 ^ chorba3); // 30
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1 ^ chorba2); // 31
        len -= sizeof(uint64_t);
        c = __crc32d(c, *buf8++ ^ chorba1); // 32
        len -= sizeof(uint64_t);
    }

    while (len >= sizeof(uint64_t)) {
        c = __crc32d(c, *buf8++);
        len -= sizeof(uint64_t);
    }

    if (len >= sizeof(uint32_t)) {
        buf4 = (const uint32_t *) buf8;
        c = __crc32w(c, *buf4++);
        len -= sizeof(uint32_t);
        buf2 = (const uint16_t *) buf4;
    } else {
        buf2 = (const uint16_t *) buf8;
    }

    if (len >= sizeof(uint16_t)) {
        c = __crc32h(c, *buf2++);
        len -= sizeof(uint16_t);
    }

    buf = (const unsigned char *) buf2;
    if (len) {
        c = __crc32b(c, *buf);
    }

    c = ~c;
    return c;
}
#endif
