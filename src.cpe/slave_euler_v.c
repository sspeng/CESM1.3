#include <slave.h>
#include "dma_macros.h"
//#include "ldm_alloc.h"

#define NP    4
#define NLEV  30
#define UC    1         // a unit that divides nete by column direction
#define UR    3         // a unit that divides qsize by row direcion
#define NC    4
#define NR    16
#define maxiter           (NP*NP - 1)
#define stripe_qdp        (NLEV*NP*NP)
#define istep_dp          (NLEV*NP*NP)
#define qstep_Qten        (NLEV*NP*NP)       // stripe in q axis of Qtens_biharmonic array
#define qstep_dp_star     (NLEV*NP*NP)
#define qstep_qmax        NLEV
#define qstep_Qtens       (NLEV*NP*NP)
#define block             (UC*NLEV*NP*NP)
#define block_dp          (NLEV*NP*NP)
#define block_vn0         (NLEV*2*NP*NP)
#define block_gradQ       (2*NP*NP)
#define block_dp_star     (NLEV*NP*NP)
#define block_Qtens     (NLEV*NP*NP)
#define block_Dinv        (2*2*NP*NP)
#define block_det         (NP*NP)

#define abs(value, ret) asm volatile ("fcpys $31, %1, %0" : "=r"(ret) : "r"(value))

#define sum_array(sum, array, len) {   \
  int i;   \
  sum = 0.0;   \
  for (i = 0; i < len; i++)     \
    sum = sum + array[i];  \
}

#define sum_array_multiply_1(sum, a, b, len, type) { \
  int i;            \
  sum = 0.0;         \
  type _array[len];    \
  for (i = 0; i < len; i++) { \
    _array[i] = a[i] * b[i];   \
    sum = sum + _array[i]; \
  }  \
}

#define sum_array_multiply(sum, a, b) { \
  int i;  \
  double _array[16];   \
  sum = 0;     \
  for (i = 0; i < 16; i++) {  \
    _array[i] = a[i] * b[i]; \
    sum = sum + _array[i];  \
  } \
}

typedef struct {
  double *gl_qdp, *gl_qdp_leap, *divdp_proj, *dp, *vn0, *Dvv, *Dinv, *metdet, \
      *rmetdet, *Qtens_biharmonic, *divdp, *dpdiss_biharmonic,   \
      *spheremp, *qmax, *qmin;
  double dt, rrearth, nu_p, nu_q;
  int nets, nete, rhs_multiplier, qsize, n0_qdp, np1_qdp, limiter_option       \
      , rhs_viss;
} param_t;

void slave_euler_v_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
	dma_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();

  double *gl_qdp = param_d.gl_qdp;
  double *gl_qdp_leap = param_d.gl_qdp_leap;
  double *gl_divdp_proj = param_d.divdp_proj;
  double *gl_dp = param_d.dp;
  double *gl_vn0 = param_d.vn0;
  double *gl_Dvv = param_d.Dvv;
  double *gl_Dinv = param_d.Dinv;
  double *gl_metdet = param_d.metdet;
  double *gl_rmetdet = param_d.rmetdet;
  double *gl_Qtens_biharmonic = param_d.Qtens_biharmonic;
  double *gl_divdp = param_d.divdp;
  double *gl_dpdiss_biharmonic = param_d.dpdiss_biharmonic;
  double *gl_spheremp = param_d.spheremp;
  double *gl_qmax = param_d.qmax;
  double *gl_qmin = param_d.qmin;
  double dt = param_d.dt;
  double rrearth = param_d.rrearth;
  double nu_p = param_d.nu_p;
  double nu_q = param_d.nu_q;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int rhs_multiplier = param_d.rhs_multiplier;
  int qsize = param_d.qsize;
  int n0_qdp = param_d.n0_qdp;
  int np1_qdp = param_d.np1_qdp;
  int limiter_option = param_d.limiter_option;
  int rhs_viss = param_d.rhs_viss;

  rid = id / NC;
  cid = id % NC;
  int loop_r = ((nete - nets + 1) + UR*NR - 1)/(UR*NR);
  int loop_c = (qsize + UC*NC - 1)/(UC*NC);
  int slice_qdp = (int)(gl_qdp_leap - gl_qdp);
  int istep_dp_star = qsize*NLEV*NP*NP;
  int istep_Qtens = qsize*NLEV*NP*NP;
  int istep_qmax = qsize*NLEV;
  int c, r, i, j, l, k, k1, iter, q, ie, cbeg, cend, rbeg, rend, cn, rn;

  double Qdp[block];
  double Qdp_np1[block];
  double dp[block_dp];
  double dp_temp[block_dp];
  double dp_star[block_dp_star];
  double vn0[block_vn0];
  double Vstar[block_vn0];                 // same size as vn0
  double divdp_proj[block_dp];             // same size as dp_tmp
  double gradQ[block_gradQ];
  double Dvv[NP*NP];
  double Dinv[2*2*NP*NP];
  double metdet[NP*NP];
  double rmetdet[NP*NP];
  double Qtens_temp[block_Qtens];
  double Qtens_biharmonic[block];
  double dpdiss_biharmonic[block_dp];     // same size as dp
  double divdp[block_dp];                 // same size as dp
  double dpdiss[NP*NP];
  double spheremp[NP*NP];
  double qmax[UC*NLEV];
  double qmin[UC*NLEV];

#if 0
  if (id == 0) {
    int size_tol = sizeof(double)*(block + block + block_dp + block_dp_star + block_vn0 \
        + block_vn0 + block_dp + block_gradQ + NP*NP + 4*NP*NP + NP*NP + NP*NP  \
        + block + block + block_dp + block_dp + NP*NP + NP*NP + NP*NP + NP*NP   \
        + UC*NLEV + UC*NLEV + 6*NP*NP);
    printf("size_tol:%dk\n", size_tol/1024);
  }
#endif

  pe_get(gl_Dvv, Dvv, NP*NP*sizeof(double));
  dma_syn();

  /* local Variables of deivergence_sphere */
  double dudx00, dvdy00;
  double gv[2*NP*NP], vvtemp[NP*NP];
  /* local Variables of deivergence_sphere */
  double cc[NP*NP], xx[NP*NP];
  double addmass, weightssum;
  double mass, sumc;
  double tol_limiter = 5e-14;

  double *src_qdp, *src_qdp_np1, *src_dp, *src_vn0, *src_divdp_proj, *src_Dinv,\
     *src_metdet, *src_rmetdet, *src_Qtens_biharmonic, *src_divdp,             \
     *src_dpdiss_biharmonic, *src_spheremp, *src_qmax, *src_qmin;              \
  double *gl_n0_qdp = gl_qdp + (n0_qdp - 1)*qsize*stripe_qdp;
  double *gl_np1_qdp = gl_qdp + (np1_qdp - 1)*qsize*stripe_qdp;

  // Divide ie-axis data on the row cpe with loop_r
  // Divide q-axis data on the cloumn cpe with loop_c
  for (r = 0; r < loop_r; r++) {
    rbeg = r*NR*UR + rid*UR;
    rend = r*NR*UR + (rid + 1)*UR;
    rend = rend < (nete - nets + 1) ? rend : (nete - nets + 1);
    rn = rend - rbeg;
    rn = rn < 0 ? 0 : rn;   // handling boundary issues, removing the case where rn < 0
    for (ie = 0; ie < rn; ie++) {
      src_dp = gl_dp + (rbeg + ie)*slice_qdp;
      src_vn0 = gl_vn0 + (rbeg + ie)*slice_qdp;
      src_divdp_proj = gl_divdp_proj + (rbeg + ie)*slice_qdp;
      src_divdp = gl_divdp + (rbeg + ie)*slice_qdp;
      src_dpdiss_biharmonic = gl_dpdiss_biharmonic + (rbeg + ie)*slice_qdp;
      src_Dinv = gl_Dinv + (rbeg + ie)*slice_qdp;
      src_metdet = gl_metdet + (rbeg + ie)*slice_qdp;
      src_rmetdet = gl_rmetdet + (rbeg + ie)*slice_qdp;
      src_spheremp = gl_spheremp + (rbeg + ie)*slice_qdp;
      pe_get(src_dp, dp, block_dp*sizeof(double));
      pe_get(src_vn0, vn0, block_vn0*sizeof(double));
      pe_get(src_divdp_proj, divdp_proj, block_dp*sizeof(double));
      pe_get(src_divdp, divdp, block_dp*sizeof(double));
      pe_get(src_dpdiss_biharmonic, dpdiss_biharmonic, block_dp*sizeof(double));
      pe_get(src_Dinv, Dinv, block_Dinv*sizeof(double));
      pe_get(src_metdet, metdet, block_det*sizeof(double));
      pe_get(src_rmetdet, rmetdet, block_det*sizeof(double));
      pe_get(src_spheremp, spheremp, NP*NP*sizeof(double));
      dma_syn();

      for (k = 0; k < NLEV; k++) {
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos_dp = k*NP*NP + j*NP + i;
            int pos_vn0_1 = k*2*NP*NP + j*NP + i;
            int pos_vn0_2 = k*2*NP*NP + NP*NP + j*NP + i;
            dp_temp[pos_dp] = dp[pos_dp] - rhs_multiplier*dt*divdp_proj[pos_dp];
            Vstar[pos_vn0_1] = vn0[pos_vn0_1]/dp_temp[pos_dp];
            Vstar[pos_vn0_2] = vn0[pos_vn0_2]/dp_temp[pos_dp];
          }
        }
      }
      for (c = 0; c < loop_c; c++) {
        cbeg = c*NC*UC + cid*UC;
        cend = c*NC*UC + (cid + 1)*UC;
        cend = cend < qsize ? cend : qsize;
        cn = cend - cbeg;
        if (cn > 0) {   // if cn < 0, the dma will get exceptional contribution
          src_qdp = gl_n0_qdp + (rbeg + ie)*slice_qdp + cbeg*stripe_qdp;
          src_qdp_np1 = gl_np1_qdp + (rbeg + ie)*slice_qdp + cbeg*stripe_qdp;
          src_Qtens_biharmonic = gl_Qtens_biharmonic + (rbeg + ie)*istep_Qtens \
              + cbeg*qstep_Qtens;
          src_qmax = gl_qmax + (rbeg + ie)*istep_qmax + cbeg*qstep_qmax;
          src_qmin = gl_qmin + (rbeg + ie)*istep_qmax + cbeg*qstep_qmax;
          pe_get(src_qdp, Qdp, block*sizeof(double));
          pe_get(src_Qtens_biharmonic, Qtens_biharmonic, block*sizeof(double));
          pe_get(src_qmax, qmax, UC*NLEV*sizeof(double));
          pe_get(src_qmin, qmin, UC*NLEV*sizeof(double));
          dma_syn();
          for (q = 0; q < cn; q++) {
            for (k = 0; k < NLEV; k++) {
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_qdp = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  int pos_vn0_1 = k*2*NP*NP + j*NP + i;
                  int pos_vn0_2 = k*2*NP*NP + NP*NP + j*NP + i;
                  int pos_gradQ_1 = j*NP + i;
                  int pos_gradQ_2 = NP*NP + j*NP + i;
                  int pos_Dinv_1 = j*NP + i;
                  int pos_Dinv_2 = NP*NP + j*NP + i;
                  int pos_Dinv_3 = 2*NP*NP + j*NP + i;
                  int pos_Dinv_4 = 3*NP*NP + j*NP + i;
                  int pos_gv_1 = j*NP + i;
                  int pos_gv_2 = NP*NP + j*NP + i;
                  int pos_det = j*NP + i;
                  gradQ[pos_gradQ_1] = Vstar[pos_vn0_1]*Qdp[pos_qdp];
                  gradQ[pos_gradQ_2] = Vstar[pos_vn0_2]*Qdp[pos_qdp];
                  gv[pos_gv_1] = metdet[pos_det]*(Dinv[pos_Dinv_1]*gradQ[pos_gradQ_1] \
                      + Dinv[pos_Dinv_3]*gradQ[pos_gradQ_2]);
                  gv[pos_gv_2] = metdet[pos_det]*(Dinv[pos_Dinv_2]*gradQ[pos_gradQ_1] \
                      + Dinv[pos_Dinv_4]*gradQ[pos_gradQ_2]);
                }
              }
              for (j = 0; j < NP; j++) {
                for (l = 0; l < NP; l++) {
                  dudx00 = 0.0;
                  dvdy00 = 0.0;
                  for (i = 0; i < NP; i++) {
                    int pos_Dvv = l*NP + i;
                    int pos_gv_1 = j*NP + i;
                    int pos_gv_2 = NP*NP + i*NP + j;
                    dudx00 = dudx00 + Dvv[pos_Dvv]*gv[pos_gv_1];
                    dvdy00 = dvdy00 + Dvv[pos_Dvv]*gv[pos_gv_2];
                  }
                  int pos_dp_star = k*NP*NP + j*NP + l;
                  int pos_vvtemp = l*NP + j;
                  dp_star[pos_dp_star] = dudx00;
                  vvtemp[pos_vvtemp] = dvdy00;
                }
              }
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_dp_star = k*NP*NP + j*NP + i;
                  int pos_vvtemp = j*NP + i;
                  int pos_det = j*NP + i;
                  dp_star[pos_dp_star] = (dp_star[pos_dp_star] + vvtemp[pos_vvtemp])   \
                      *(rmetdet[pos_det]*rrearth);
                }
              }   // end of divergence_sphere
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_Qtens = k*NP*NP + j*NP + i;
                  int pos_dp_star = k*NP*NP + j*NP + i;
                  int pos_qdp = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  Qtens_temp[pos_Qtens] = Qdp[pos_qdp] - dt*dp_star[pos_dp_star];
                }
              }
              if (rhs_viss != 0) {
                for (j = 0; j < NP; j++) {
                  for (i = 0; i < NP; i++) {
                    int pos_Qtens = k*NP*NP + j*NP + i;
                    int pos_bi = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                    Qtens_temp[pos_Qtens] = Qtens_temp[pos_Qtens] + Qtens_biharmonic[pos_bi];
                  }
                }
              }
            }   // end loop k
            if (limiter_option == 8) {
              for (k = 0; k < NLEV; k++) {
                for (j = 0; j < NP; j++) {
                  for (i = 0; i < NP; i++) {
                    int pos_dp_star = k*NP*NP + j*NP + i;
                    int pos_dp = k*NP*NP + j*NP + i;
                    dp_star[pos_dp_star] = dp_temp[pos_dp] - dt*divdp[pos_dp];
                  }
                }
                if (nu_p > 0 && rhs_viss != 0) {
                  for (j = 0; j < NP; j++) {
                    for (i = 0; i < NP; i++) {
                      int pos_dpdiss = j*NP + i;
                      int pos_spheremp = j*NP + i;
                      int pos_dpdiss_bi = k*NP*NP + j*NP + i;
                      int pos_dp_star = k*NP*NP + j*NP + i;
                      dpdiss[pos_dpdiss] = dpdiss_biharmonic[pos_dpdiss_bi];
                      dp_star[pos_dp_star] = dp_star[pos_dp_star]              \
                          - rhs_viss*dt*nu_q*dpdiss[pos_dpdiss]/spheremp[pos_spheremp];
                    }
                  }
                }
              }
              // start of function limiter_optim_iter_full
              for (k = 0; k < NLEV; k++) {
                for (k1 = 0; k1 < NP*NP; k1++) {
                  int pos_dp_star = k*NP*NP + k1;
                  int pos_Qtens = k*NP*NP + k1;
                  cc[k1] = spheremp[k1]*dp_star[pos_dp_star];
                  xx[k1] = Qtens_temp[pos_Qtens]/dp_star[pos_dp_star];
                }  // end loop k
                // implement the function sum which same as sum(c) in fortran
                // implement the function sum(c*x)
                // sum_array_multiply(mass, cc, xx, NP*NP, double)
                sum_array(sumc, cc, NP*NP);
                if (sumc <= 0) continue;
                sum_array_multiply(mass, cc, xx);
                int pos_qmax = q*NLEV + k;
                if (mass < qmin[pos_qmax]*sumc) qmin[pos_qmax] = mass/sumc;
                if (mass > qmax[pos_qmax]*sumc) qmax[pos_qmax] = mass/sumc;
                for (iter = 0; iter < maxiter; iter++) {
                  addmass = 0.0;
                  for (k1 = 0; k1 < NP*NP; k1++) {
                    if(xx[k1] > qmax[pos_qmax]) {
                      addmass = addmass + (xx[k1] - qmax[pos_qmax])*cc[k1];
                      xx[k1] = qmax[pos_qmax];
                    }
                    if(xx[k1] < qmin[pos_qmax]) {
                      addmass = addmass - (qmin[pos_qmax] - xx[k1])*cc[k1];
                      xx[k1] = qmin[pos_qmax];
                    }
                  }  // end loop k1
                  double addmass_abs;
                  double mass_abs;
                  abs(mass, mass_abs);
                  abs(addmass, addmass_abs);
                  if (addmass_abs <= tol_limiter*mass_abs) break;
                  weightssum = 0.0;
                  if (addmass > 0) {
                    for (k1 = 0; k1 < NP*NP; k1++) {
                      if (xx[k1] < qmax[pos_qmax])
                        weightssum = weightssum + cc[k1];
                    }

                    for (k1 = 0; k1 < NP*NP; k1++) {
                      if (xx[k1] < qmax[pos_qmax])
                        xx[k1] = xx[k1] + addmass/weightssum;
                    }

                  } else {
                    for (k1 = 0; k1 < NP*NP; k1++) {
                      if (xx[k1] > qmin[pos_qmax])
                        weightssum = weightssum + cc[k1];
                    }

                    for (k1 = 0; k1 < NP*NP; k1++) {
                      if (xx[k1] > qmin[pos_qmax])
                        xx[k1] = xx[k1] + addmass/weightssum;
                    }

                  }
                }  // end loop iter
                for (k1 = 0; k1 < NP*NP; k1++) {
                  int pos_Qtens = k*NP*NP + k1;
                  Qtens_temp[pos_Qtens] = xx[k1];
                }
              }   // end loop k
              for (k = 0; k < NLEV; k++) {
                for (k1 = 0; k1 < NP*NP; k1++) {
                  int pos_Qtens = k*NP*NP + k1;
                  int pos_dp_star = k*NP*NP + k1;
                  Qtens_temp[pos_Qtens] = Qtens_temp[pos_Qtens]*dp_star[pos_dp_star];
                }
              }  // end function limiter_optim_iter_full
            }   // end if limiter_option = 8
            for (k = 0; k < NLEV; k++) {
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_Qtens = k*NP*NP + j*NP + i;
                  int pos = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  int pos_spheremp = j*NP + i;
                  Qdp_np1[pos] = spheremp[pos_spheremp]*Qtens_temp[pos_Qtens];
                }
              }
            }
          }     // end loop q
          pe_put(src_qdp_np1, Qdp_np1, cn*NLEV*NP*NP*sizeof(double));
          pe_put(src_qmax, qmax, cn*NLEV*sizeof(double));
          pe_put(src_qmin, qmin, cn*NLEV*sizeof(double));
          dma_syn();
        }
      }
    }
  }

#if 0
  if (id == 0) {
    //double abs_org = -6.03504589485943894689494646;
    //double abs_test1 = 6.454934674584580238594;
    //double abs_aft;
    //abs(abs_org, abs_aft);
    //printf("abs_org:%.20lf, abs_aft:%.20lf\n", abs_org, abs_aft);
    double a;
    printf("a:%lf\n", a);
    a = a + 1;
    printf("a:%lf\n", a);
  }
#endif
}
