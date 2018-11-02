#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define UC 7         // a unit that divides nete by column direction
#define UR 3         // a unit that divides qsize by row direcion
#define NC 4
#define NR 16
#define block_Dinv        (2*2*NP*NP)
#define block_tensor      (2*2*NP*NP)
#define block_qtens       (UC*NLEV*NP*NP)
#define qstep_qtens       (NLEV*NP*NP)

#define abs(value, ret) asm volatile ("fcpys $31, %1, %0" : "=r"(ret) : "r"(value))

typedef struct {
  double *Dinv, *rspheremp, *spheremp, *variable_hyper, *tensorVisc, *qtens, *Dvv;
  double rrearth, hypervis_scaling, hypervis_power;
  int nets, nete, qsize, step_elem;
} param_t;

void slave_biharmonic_wk_scalar_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();

  double *gl_Dinv = param_d.Dinv;
  double *gl_rspheremp = param_d.rspheremp;
  double *gl_spheremp = param_d.spheremp;
  double *gl_variable_hyper = param_d.variable_hyper;
  double *gl_tensorVisc = param_d.tensorVisc;
  double *gl_qtens = param_d.qtens;
  double *gl_Dvv = param_d.Dvv;
  double rrearth = param_d.rrearth;
  double hypervis_scaling = param_d.hypervis_scaling;
  double hypervis_power = param_d.hypervis_power;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;
  int step_elem = param_d.step_elem;

  double Dinv[block_Dinv];
  double Dvv[NP*NP];
  double variable_hyper[NP*NP];
  double tensorVisc[block_tensor];
  double qtens[block_qtens];
  double rspheremp[NP*NP];
  double spheremp[NP*NP];
  double lap_p[NP*NP];
  double grads[2*NP*NP];
  double v1[NP*NP];
  double v2[NP*NP];
  double vtemp[NP*NP*2];

  double dsdx00, dsdy00; int var_coef = 1;
  double oldgrads_1, oldgrads_2, tensor_1, tensor_2;
  if (hypervis_scaling > 0) var_coef = 0;

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block_Dinv + block_tensor + block_qtens     \
        + 11*NP*NP);
    printf("slave_euler_step size_tol:%dk\n", size_tol/1024);
  }
#endif

  pe_get(gl_Dvv, Dvv, NP*NP*sizeof(double));
  dma_syn();

  double *src_Dinv, *src_rspheremp, *src_spheremp, *src_variable_hyper,        \
      *src_tensorVisc, *src_qtens;

  rid = id / NC;
  cid = id % NC;
  int loop_r = ((nete - nets + 1) + UR*NR - 1)/(UR*NR);
  int loop_c = (qsize + UC*NC - 1)/(UC*NC);
  int istep_qtens = qsize*NLEV*NP*NP;
  int c, r, i, j, d, l, m, n, k, q, ie, cbeg, cend, rbeg, rend, cn, rn, pos_dp \
      , pos_qdp, pos_Qtens_bi, pos_qmax;

  // Divide ie-axis data on the row cpe with loop_r
  // the value of loop_r rely on NR, UR; NP is the number of cloumn cpe,
  // UR is the unit that divides ie size by cloumn direcion,
  // Divide q-axis data on the cloumn cpe with loop_c
  // same as loop_r, the loop_c rely on NC, UC
  // UC is the unit that divides q size by row direction
  for (r = 0; r < loop_r; r++) {
    rbeg = r*NR*UR + rid*UR;
    rend = r*NR*UR + (rid + 1)*UR;
    rend = rend < (nete - nets + 1) ? rend : (nete - nets + 1);
    rn = rend - rbeg;
    rn = rn < 0 ? 0 : rn;  // handling boundary issues, removing the case where rn < 0
    for (ie = 0; ie < rn; ie++) {
      src_Dinv = gl_Dinv + (rbeg + ie)*step_elem;
      src_rspheremp = gl_rspheremp + (rbeg + ie)*step_elem;
      src_spheremp = gl_spheremp + (rbeg + ie)*step_elem;
      src_variable_hyper = gl_variable_hyper                 \
          + (rbeg + ie)*step_elem;
      src_tensorVisc = gl_tensorVisc + (rbeg + ie)*step_elem;
      pe_get(src_Dinv, Dinv, block_Dinv*sizeof(double));
      pe_get(src_rspheremp, rspheremp, NP*NP*sizeof(double));
      pe_get(src_spheremp, spheremp, NP*NP*sizeof(double));
      pe_get(src_tensorVisc, tensorVisc, block_tensor*sizeof(double));
      pe_get(src_variable_hyper, variable_hyper, NP*NP*sizeof(double));
      dma_syn();
      for (c = 0; c < loop_c; c++) {
        cbeg = c*NC*UC + cid*UC;
        cend = c*NC*UC + (cid + 1)*UC;
        cend = cend < qsize ? cend : qsize;
        cn = cend - cbeg;
        if (cn > 0) {  // if cn < 0, the dma will get exceptional contribution
          src_qtens = gl_qtens + (rbeg + ie)*istep_qtens + cbeg*qstep_qtens;
          pe_get(src_qtens, qtens, block_qtens*sizeof(double));
          dma_syn();
          for (q = 0; q < cn; q++) {
            for (k = 0; k < NLEV; k++) {
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_lap = j*NP + i;
                  int pos_qtens = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  lap_p[pos_lap] = qtens[pos_qtens];
                }
              }
              for (j = 0; j < NP; j++) {  // start gradient_sphere
                for (l = 0; l < NP; l++) {
                  dsdx00 = 0.0;
                  dsdy00 = 0.0;
                  for (i = 0; i < NP; i++) {
                    int pos_Dvv = l*NP + i;
                    int pos_lap_1 = j*NP + i;
                    int pos_lap_2 = i*NP + j;
                    dsdx00 = dsdx00 + Dvv[pos_Dvv]*lap_p[pos_lap_1];
                    dsdy00 = dsdy00 + Dvv[pos_Dvv]*lap_p[pos_lap_2];
                  }
                  int pos_v1 = j*NP + l;
                  int pos_v2 = l*NP + j;
                  v1[pos_v1] = dsdx00*rrearth;
                  v2[pos_v2] = dsdy00*rrearth;
                }
              }
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_grads_1 = j*NP + i;
                  int pos_grads_2 = NP*NP + j*NP + i;
                  int pos_v = j*NP + i;
                  int pos_Dinv_1 = j*NP + i;
                  int pos_Dinv_2 = NP*NP + j*NP + i;
                  int pos_Dinv_3 = 2*NP*NP + j*NP + i;
                  int pos_Dinv_4 = 3*NP*NP + j*NP + i;
                  grads[pos_grads_1] = Dinv[pos_Dinv_1]*v1[pos_v]              \
                      + Dinv[pos_Dinv_2]*v2[pos_v];
                  grads[pos_grads_2] = Dinv[pos_Dinv_3]*v1[pos_v]              \
                      + Dinv[pos_Dinv_4]*v2[pos_v];
                }
              }
              if (var_coef) {
                //if (id == 0)
                //  printf("var_coef_in:%d\n", var_coef);
                if (hypervis_power != 0) {
                  for (j = 0; j < NP; j++) {
                    for (i = 0; i < NP; i++) {
                      int pos_grads_1 = j*NP + i;
                      int pos_grads_2 = NP*NP + j*NP + i;
                      int pos = j*NP + i;
                      grads[pos_grads_1] = grads[pos_grads_1]*variable_hyper[pos];
                      grads[pos_grads_2] = grads[pos_grads_2]*variable_hyper[pos];
                    }
                  }
                } else if (hypervis_scaling != 0) {
                  for (j = 0; j < NP; j++) {
                    for (i = 0; i < NP; i++) {
                      int pos_grads_1 = j*NP + i;
                      int pos_grads_2 = NP*NP + j*NP + i;
                      int pos_tensor_1 = j*NP + i;
                      int pos_tensor_2 = NP*NP + j*NP + i;
                      int pos_tensor_3 = 2*NP*NP + j*NP + i;
                      int pos_tensor_4 = 3*NP*NP + j*NP + i;
                      oldgrads_1 = grads[pos_grads_1];
                      oldgrads_2 = grads[pos_grads_2];
                      grads[pos_grads_1] = oldgrads_1*tensorVisc[pos_tensor_1] \
                          + oldgrads_2*tensorVisc[pos_tensor_3];
                      grads[pos_grads_2] = oldgrads_1*tensorVisc[pos_tensor_2] \
                          + oldgrads_2*tensorVisc[pos_tensor_4];
                    }
                  }
                }
              }  // end of gradient_sphere
              // start divergence_sphere_wk
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_1 = j*NP + i;
                  int pos_2 = NP*NP + j*NP + i;
                  int pos_Dinv_1 = j*NP + i;
                  int pos_Dinv_2 = NP*NP + j*NP + i;
                  int pos_Dinv_3 = 2*NP*NP + j*NP + i;
                  int pos_Dinv_4 = 3*NP*NP + j*NP + i;
                  vtemp[pos_1] = Dinv[pos_Dinv_1]*grads[pos_1]                 \
                      + Dinv[pos_Dinv_3]*grads[pos_2];
                  vtemp[pos_2] = Dinv[pos_Dinv_2]*grads[pos_1]                 \
                      + Dinv[pos_Dinv_4]*grads[pos_2];
                }
              }
              for (n = 0; n < NP; n++) {
                for (m = 0; m < NP; m++) {
                  int pos_qtens = q*NLEV*NP*NP + k*NP*NP + n*NP + m;
                  qtens[pos_qtens] = 0;
                  for (j = 0; j < NP; j++) {
                    int pos_spheremp_1 = n*NP + j;
                    int pos_spheremp_2 = j*NP + m;
                    int pos_vtemp_1 = n*NP + j;
                    int pos_vtemp_2 = NP*NP + j*NP + m;
                    int pos_Dvv_1 = j*NP + m;
                    int pos_Dvv_2 = j*NP + n;
                    qtens[pos_qtens] = qtens[pos_qtens]
                        - (spheremp[pos_spheremp_1]*vtemp[pos_vtemp_1]*Dvv[pos_Dvv_1]         \
                        + spheremp[pos_spheremp_2]*vtemp[pos_vtemp_2]*Dvv[pos_Dvv_2])         \
                        *rrearth;
                  }
                }
              }  // end divergence_sphere_wk
            }  // end loop k
          }  // end loop q
          pe_put(src_qtens, qtens, cn*NLEV*NP*NP*sizeof(double));
          dma_syn();
        }  // end if cn > 0
      }  // end loop_c
    }  // end ie
  }  // end loop_r

#if 0
  if (id == 0)
    printf("qsize:%d, var_coef:%d\n", qsize, var_coef);
#endif
}
