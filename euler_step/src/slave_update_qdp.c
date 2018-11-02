#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define UC 4         // a unit that divides nete by column direction
#define UR 3         // a unit that divides qsize by row direcion
#define NC 8
#define NR 8
#define stripe_qdp        (NLEV*NP*NP)
#define qstep_Qten (NLEV*NP*NP)       // stripe in q axis of Qtens_biharmonic array
#define block  (UC*NLEV*NP*NP)


#define get_row_id(rid) asm volatile ("rcsr %0, 1" : "=r"(rid))
#define get_col_id(cid) asm volatile ("rcsr %0, 2" : "=r"(cid))
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0":"=r"(var))
#define max(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define min(val1, val2) ((val1) < (val2) ? (val1) : (val2))


typedef struct {
  double *gl_Qdp, *gl_Qtens_temp, *gl_spheremp;
  int nets, nete, qsize, qsize_d, step_elem;
} param_t;

void slave_update_qdp_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
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

  rid = id / NC;
  cid = id % NC;
  int loop_r = ((nete - nets + 1) + UR*NR - 1)/(UR*NR);
  int loop_c = (qsize + UC*NC - 1)/(UC*NC);
  int istep_Qtens = qsize*NLEV*NP*NP;
  int c, r, i, j, k, q, ie, cbeg, cend, rbeg, rend, cn, rn, pos_dp, pos_qdp    \
      , pos_Qtens_bi, pos_qmax;

  // Divide ie-axis data on the row cpe with loop_r
  // Divide q-axis data on the cloumn cpe with loop_c
  for (r = 0; r < loop_r; r++) {
    rbeg = r*NR*UR + rid*UR;
    rend = r*NR*UR + (rid + 1)*UR;
    rend = rend < (nete - nets + 1) ? rend : (nete - nets + 1);
    rn = rend - rbeg;
    rn = rn < 0 ? 0 : rn;  // handling boundary issues, removing the case where rn < 0
    for (ie = 0; ie < rn; ie++) {
      src_spheremp = gl_spheremp + (rbeg + ie)*step_elem;
      pe_get(src_spheremp, spheremp, NP*NP*sizeof(double));
      dma_syn();
      for (c = 0; c < loop_c; c++) {
        cbeg = c*NC*UC + cid*UC;
        cend = c*NC*UC + (cid + 1)*UC;
        cend = cend < qsize ? cend : qsize;
        cn = cend - cbeg;
        if (cn > 0) {  // if cn < 0, the dma will get exceptional contribution
          src_Qdp = gl_Qdp + (rbeg + ie)*step_elem + cbeg*stripe_qdp;
          src_Qtens_temp = gl_Qtens_temp + (rbeg + ie)*istep_Qtens + cbeg*qstep_Qten;
          pe_get(src_Qdp, Qdp, block*sizeof(double));
          pe_get(src_Qtens_temp, Qtens_temp, (block*sizeof(double)));
          dma_syn();
          for (q = 0; q < cn; q++) {
            for (k = 0; k < NLEV; k++) {
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  int pos_spheremp = j*NP + i;
                  Qdp[pos] = spheremp[pos_spheremp]*Qtens_temp[pos];
                }
              }
            }
          }
          pe_put(src_Qdp, Qdp, (cn*NLEV*NP*NP*sizeof(double)));
          dma_syn();
        }
      }
    }
  }
#if 0
  if (id == 0) {
    printf("rhs_viss:%d, dt:%lf, nu_q:%lf\n", rhs_viss, dt, nu_q);
  }
#endif
}
