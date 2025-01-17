/* stb_image_write - v1.16 - public domain - http://nothings.org/stb
   Image writer for C/C++: PNG, BMP, TGA, JPG, HDR

   Standard library: <stdlib.h>, <stdio.h>, <string.h>, <math.h>, <assert.h>
*/

#ifndef INCLUDE_STB_IMAGE_WRITE_H
#define INCLUDE_STB_IMAGE_WRITE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef STB_IMAGE_WRITE_STATIC
#define STBIWDEF static
#else
#define STBIWDEF extern
#endif

#ifndef STBIW_UCHAR
typedef unsigned char stbiw_uint8;
#else
typedef STBIW_UCHAR stbiw_uint8;
#endif

#ifndef STBIW_UINT32
typedef unsigned int stbiw_uint32;
#endif

#ifndef STBIW_INT32
typedef int stbiw_int32;
#endif

typedef void stbi_write_func(void *context, void *data, int size);

STBIWDEF int stbi_write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes);
STBIWDEF int stbi_write_png_to_func(stbi_write_func *func, void *context, int w, int h, int comp, const void *data, int stride_in_bytes);
STBIWDEF void stbi_flip_vertically_on_write(int flip_boolean);

#ifdef __cplusplus
}
#endif

#endif // INCLUDE_STB_IMAGE_WRITE_H

#ifdef STB_IMAGE_WRITE_IMPLEMENTATION

#ifdef __cplusplus
extern "C" {
#endif

static int stbi__flip_vertically_on_write = 0;

STBIWDEF void stbi_flip_vertically_on_write(int flip_boolean)
{
   stbi__flip_vertically_on_write = flip_boolean;
}

typedef struct
{
   stbi_write_func *func;
   void *context;
   unsigned char buffer[64];
   int buf_used;
} stbi__write_context;

static void stbi__wbuf_flush(stbi__write_context *s)
{
   if (s->buf_used) {
      s->func(s->context, s->buffer, s->buf_used);
      s->buf_used = 0;
   }
}

static void stbi__wbuf(stbi__write_context *s, const void *d, int size)
{
   stbi__wbuf_flush(s);
   s->func(s->context, (void*)d, size);
}

static void stbi__write_flush(stbi__write_context *s)
{
   stbi__wbuf_flush(s);
}

static void stbi__stdio_write(void *context, void *data, int size)
{
   fwrite(data, 1, size, (FILE*) context);
}

// CRC32 and Adler32 calculation for PNG
static unsigned int stbi__crc32(const unsigned char *buffer, int len)
{
   static unsigned int crc_table[256];
   static int crc_table_computed = 0;
   unsigned int crc = 0xffffffff;
   int i, j;

   if (!crc_table_computed) {
      for (i = 0; i < 256; i++) {
         unsigned int c = i;
         for (j = 0; j < 8; j++) {
            if (c & 1)
               c = 0xedb88320L ^ (c >> 1);
            else
               c = c >> 1;
         }
         crc_table[i] = c;
      }
      crc_table_computed = 1;
   }

   for (i = 0; i < len; i++)
      crc = crc_table[(crc ^ buffer[i]) & 0xff] ^ (crc >> 8);

   return ~crc;
}

static unsigned int stbi__adler32(const unsigned char *buffer, int len)
{
   unsigned int s1 = 1, s2 = 0;
   int i;
   for (i = 0; i < len; ++i) {
      s1 = (s1 + buffer[i]) % 65521;
      s2 = (s2 + s1) % 65521;
   }
   return (s2 << 16) | s1;
}

// Minimal DEFLATE without compression (store blocks) for reliable PNG encoding
static unsigned char *stbi__deflate_store(const unsigned char *data, int len, int *out_len)
{
   int blocks = (len + 65534) / 65535;
   int total_len = 2 + blocks * 5 + len + 4; // zlib header (2) + blocks + adler32 (4)
   unsigned char *out = (unsigned char *) malloc(total_len);
   int p = 0;
   int remaining = len;
   int src_pos = 0;
   unsigned int adler = stbi__adler32(data, len);

   if (!out) return NULL;

   // zlib header
   out[p++] = 0x78; // CMF: Deflate, 32K window
   out[p++] = 0x01; // FLG: No compression

   while (remaining > 0) {
      int block_len = remaining > 65535 ? 65535 : remaining;
      int is_last = (block_len == remaining);
      out[p++] = (unsigned char)(is_last ? 0x01 : 0x00); // BFINAL and BTYPE=00 (uncompressed)
      out[p++] = (unsigned char)(block_len & 0xff);
      out[p++] = (unsigned char)((block_len >> 8) & 0xff);
      out[p++] = (unsigned char)((~block_len) & 0xff);
      out[p++] = (unsigned char)(((~block_len) >> 8) & 0xff);
      memcpy(out + p, data + src_pos, block_len);
      p += block_len;
      src_pos += block_len;
      remaining -= block_len;
   }

   // Adler32 (Big-Endian)
   out[p++] = (unsigned char)((adler >> 24) & 0xff);
   out[p++] = (unsigned char)((adler >> 16) & 0xff);
   out[p++] = (unsigned char)((adler >> 8) & 0xff);
   out[p++] = (unsigned char)(adler & 0xff);

   *out_len = p;
   return out;
}

static void stbi__write_chunk(stbi__write_context *s, const char *type, const unsigned char *data, int len)
{
   unsigned char head[8];
   unsigned int crc;
   head[0] = (unsigned char)((len >> 24) & 0xff);
   head[1] = (unsigned char)((len >> 16) & 0xff);
   head[2] = (unsigned char)((len >> 8) & 0xff);
   head[3] = (unsigned char)(len & 0xff);
   memcpy(head + 4, type, 4);
   stbi__wbuf(s, head, 8);

   if (data && len > 0) {
      stbi__wbuf(s, data, len);
      unsigned char *combined = (unsigned char *)malloc(4 + len);
      if (combined) {
         memcpy(combined, type, 4);
         memcpy(combined + 4, data, len);
         crc = stbi__crc32(combined, 4 + len);
         free(combined);
      } else {
         crc = 0;
      }
   } else {
      crc = stbi__crc32((const unsigned char *)type, 4);
   }
   unsigned char tail[4];
   tail[0] = (unsigned char)((crc >> 24) & 0xff);
   tail[1] = (unsigned char)((crc >> 16) & 0xff);
   tail[2] = (unsigned char)((crc >> 8) & 0xff);
   tail[3] = (unsigned char)(crc & 0xff);
   stbi__wbuf(s, tail, 4);
}

static int stbi_write_png_core(stbi__write_context *s, int x, int y, int comp, const void *data, int stride_bytes)
{
   if (stride_bytes == 0) stride_bytes = x * comp;

   // PNG signature
   unsigned char sig[8] = { 137, 80, 78, 71, 13, 10, 26, 10 };
   stbi__wbuf(s, sig, 8);

   // IHDR
   unsigned char ihdr[13];
   ihdr[0] = (unsigned char)((x >> 24) & 0xff);
   ihdr[1] = (unsigned char)((x >> 16) & 0xff);
   ihdr[2] = (unsigned char)((x >> 8) & 0xff);
   ihdr[3] = (unsigned char)(x & 0xff);
   ihdr[4] = (unsigned char)((y >> 24) & 0xff);
   ihdr[5] = (unsigned char)((y >> 16) & 0xff);
   ihdr[6] = (unsigned char)((y >> 8) & 0xff);
   ihdr[7] = (unsigned char)(y & 0xff);
   ihdr[8] = 8; // 8-bit depth
   ihdr[9] = (unsigned char)(comp == 1 ? 0 : (comp == 2 ? 4 : (comp == 3 ? 2 : 6))); // Color type
   ihdr[10] = 0; // Compression (deflate)
   ihdr[11] = 0; // Filter
   ihdr[12] = 0; // Interlace
   stbi__write_chunk(s, "IHDR", ihdr, 13);

   // Uncompressed raw scanlines with filter type 0 (None)
   int raw_stride = 1 + x * comp;
   int raw_len = y * raw_stride;
   unsigned char *raw_data = (unsigned char *)malloc(raw_len);
   if (!raw_data) return 0;

   const unsigned char *src = (const unsigned char *)data;
   for (int j = 0; j < y; ++j) {
      int row = stbi__flip_vertically_on_write ? (y - 1 - j) : j;
      unsigned char *dst_row = raw_data + j * raw_stride;
      dst_row[0] = 0; // Filter: None
      memcpy(dst_row + 1, src + row * stride_bytes, x * comp);
   }

   // Compress with zlib deflate
   int zlib_len = 0;
   unsigned char *zlib_data = stbi__deflate_store(raw_data, raw_len, &zlib_len);
   free(raw_data);

   if (!zlib_data) return 0;

   // IDAT
   stbi__write_chunk(s, "IDAT", zlib_data, zlib_len);
   free(zlib_data);

   // IEND
   stbi__write_chunk(s, "IEND", NULL, 0);

   stbi__write_flush(s);
   return 1;
}

STBIWDEF int stbi_write_png_to_func(stbi_write_func *func, void *context, int x, int y, int comp, const void *data, int stride_bytes)
{
   stbi__write_context s;
   s.func = func;
   s.context = context;
   s.buf_used = 0;
   return stbi_write_png_core(&s, x, y, comp, data, stride_bytes);
}

STBIWDEF int stbi_write_png(char const *filename, int x, int y, int comp, const void *data, int stride_bytes)
{
   FILE *f;
   if (!(f = fopen(filename, "wb"))) return 0;
   int r = stbi_write_png_to_func(stbi__stdio_write, (void *)f, x, y, comp, data, stride_bytes);
   fclose(f);
   return r;
}

#ifdef __cplusplus
}
#endif

#endif // STB_IMAGE_WRITE_IMPLEMENTATION
