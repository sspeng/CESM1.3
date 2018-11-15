#include <stdio.h>
#include <stdlib.h>
#define MAX(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define MIN(val1, val2) ((val1) < (val2) ? (val1) : (val2))
#define abs(value, ret) asm volatile ("fcpys $31, %1, %0" : "=r"(ret) : "r"(value))
#define sum(sum, array, len) {   \
  int i;   \
  for (i = 0; i < len; i++)     \
    sum = sum + array[i];  \
}

#define sum_array_multiply_1(sum, a, b, len, type) { \
  int i;            \
  type _array[len];    \
  for (i = 0; i < len; i++) { \
    _array[i] = a[i] * b[i];   \
    sum = sum + _array[i]; \
  }  \
}

#define sum_array_multiply(sum, a, b) { \
  int i;  \
  double _array[16];   \
  for (i = 0; i < 16; i++) {  \
    _array[i] = a[i] * b[i]; \
    sum += _array[i];  \
  } \
}

#define MAX4(a, a1, a2, a3, a4) { \
  double max = a; \
  max = MAX(max, a); \
  max = MAX(max, a1); \
  max = MAX(max, a2); \
  max = MAX(max, a3); \
  max = MAX(max, a4); \
  a = max; \
}

#define MIN4(a, a1, a2, a3, a4) { \
  double min = a; \
  min = MIN(min, a); \
  min = MIN(min, a1); \
  min = MIN(min, a2); \
  min = MIN(min, a3); \
  min = MIN(min, a4); \
  a = min; \
}
int main() {
  double a = 0.5, a1 = 10.1, a2 = 2.0, a3 = 8.3, a4 = 1.9;
  //MAX4(a, a1, a2, a3, a4);
  //printf("max:%lf\n",a);
  MIN4(a, a1, a2, a3, a4);
  printf("min:%lf\n",a);
  return 0;
}
