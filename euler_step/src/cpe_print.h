#include <stdarg.h>
static void cpe_printf(const char *fmt, ...){
  volatile long vprintf_addr = (long)vprintf;
  int (*vprintf_ptr)(const char *, va_list) = (void*)vprintf_addr;
  va_list vlist;
  va_start(vlist, fmt);
  vprintf_ptr(fmt, vlist);
  va_end(vlist);
}
#define debug_int(fmt, ...) {if (DEBUG_COND) cpe_printf("%s:%d:" fmt, __VA_ARGS__)};
