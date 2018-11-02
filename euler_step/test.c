#include <stdio.h>
#include <stdlib.h>
#define max(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define min(val1, val2) ((val1) < (val2) ? (val1) : (val2))
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
int main() {
  //double a, b, max, min;
  //a = 7.0f; b = 9.0f;
  //max = max(a, b);
  //min = min(a, b);
  //printf("max:%lf, min:%lf\n", max, min);
  //return 0;
  //int a[10][10];
  //int i, j, k;
  //for (i = 0; i < 10; i++) {
  //  for (j = 0; j < 10; j++) {
  //    a[i][j] = 10*i + j;
  //  }
  //}
  //for (i = 0; i < 10; i++) {
  //  for (j = 0; j < 10; j++) {
  //    if (j > 3 && j < 5)
  //      continue;
  //    printf("%d, ", a[i][j]);
  //    if (j > 6)
  //      break;
  //  }
  //  printf("\n");
  //}
  //return 0;
  //int sum_a = 0, sum_b = 0, sum_c = 0;
  //int a[16], b[16], c[16];
  //int test = -6;
  //int i = 0;
  //for (i = 0; i < 16; i++) {
  //  a[i] = 3;
  //  b[i] = 5;
  //}
  //for (i = 0; i < 16; i++) {
  //  sum_b = sum_b + b[i];
  //}
  ////test = abs(test);
  //sum(sum_a, a, 16);
  //for (i = 0; i < 16; i++) {
  //  c[i] = a[i] * b[i];
  //  sum_c = sum_c + c[i];
  //}
  //sum_c = 0;
  //int len = 16;
  //sum_array_multiply_1(sum_c, a, b, len, int);
  //double mass;
  //double cc[16], xx[16];
  //for (i = 0; i < 16; i++) {
  //  cc[i] = 6;
  //  xx[i] = 5;
  //}
  //sum_array_multiply(mass, cc, xx);
//
  //printf("sum_a:%d, sum_c:%d, mass:%lf\n", sum_a, sum_c, mass);
  //double a;
  //a = a + 1;
  //printf("a:%lf\n", a);
  //return 0;
  int a = 0;
  if(a) {
    printf("iam\n");
  }
}
