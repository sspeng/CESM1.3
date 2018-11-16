#include <slave.h>
#include <stdio.h>
#include <stdarg.h>

#include "cpe_print.h"

void test_(){
  if (_MYID == 0)
    cpe_printf("dxh cpe test!!!\n");
}
