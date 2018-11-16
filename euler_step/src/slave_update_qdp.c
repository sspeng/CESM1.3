#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define stripe_qdp        (NLEV*NP*NP)
#define qstep_Qten (NLEV*NP*NP)       // stripe in q axis of Qtens_biharmonic array
#define block  (12*NP*NP)


#define get_row_id(rid) asm volatile ("rcsr %0, 1" : "=r"(rid))
#define get_col_id(cid) asm volatile ("rcsr %0, 2" : "=r"(cid))
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0":"=r"(var))
#define MAX(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define MIN(val1, val2) ((val1) < (val2) ? (val1) : (val2))


typedef struct {
  double *gl_Qdp, *gl_Qtens_temp, *gl_spheremp;
  int nets, nete, qsize, qsize_d, step_elem;
} param_t;

void slave_update_qdp_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
  get_col_id(cid);
  get_row_id(rid);
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();
  double *gl_Qdp = param_d.gl_Qdp;
  double *gl_Qtens_temp = param_d.gl_Qtens_temp;
  double *gl_spheremp = param_d.gl_spheremp;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;
  int qsize_d = param_d.qsize_d;
  int step_elem = param_d.step_elem;

  double Qdp[block];
  double Qtens_temp[block];
  double spheremp[NP*NP];

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block + NLEV*NP*NP*3 + block + NLEV*2       \
        + 2*UC*NLEV);
    printf("slave_euler_step size_tol:%dk\n", size_tol/1024);
  }
#endif

  double *src_Qdp, *src_Qtens_temp, *src_spheremp;

}
