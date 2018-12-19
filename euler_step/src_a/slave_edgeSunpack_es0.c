#include <slave.h>
#include "dma_macros.h"
#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define UC 7         // a unit that divides nete by column direction
#define UR 3         // a unit that divides qsize by row direcion
#define NC 4
#define NR 16
#define west  1
#define east  2
#define south 3
#define north 4
#define swest 5
#define seast 6
#define nwest 7
#define neast 8
#define max_neigh_edges   8
#define max_corner_elem   1
#define block_Dinv        (2*2*NP*NP)
#define block_tensor      (2*2*NP*NP)
#define block_qtens       (UC*NLEV*NP*NP)
#define qstep_qtens       (NLEV*NP*NP)

#define MAX(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define MIN(val1, val2) ((val1) < (val2) ? (val1) : (val2))

#define MAX5(a, a1, a2, a3, a4) { \
  double max = a; \
  max = MAX(max, a); \
  max = MAX(max, a1); \
  max = MAX(max, a2); \
  max = MAX(max, a3); \
  max = MAX(max, a4); \
  a = max; \
}

#define MIN5(a, a1, a2, a3, a4) { \
  double min = a; \
  min = MIN(min, a); \
  min = MIN(min, a1); \
  min = MIN(min, a2); \
  min = MIN(min, a3); \
  min = MIN(min, a4); \
  a = min; \
}

typedef struct {
  double *qmax, *qmin, *receive;
  int *getmap;
  int nets, nete, qsize;
} param_t;

void slave_edgesunpack_es_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();
  double *gl_qmax = param_d.qmax;
  double *gl_qmin = param_d.qmin;
  double *gl_receive = param_d.receive;
  int *gl_getmap = param_d.getmap;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;

  double qmax[NLEV];
  double qmin[NLEV];
  double receive[NLEV];
  double receive_is[NLEV];
  double receive_ie[NLEV];
  double receive_in[NLEV];
  double receive_iw[NLEV];
  int getmap[max_neigh_edges];

  int k, n, ie, q, id_map, kptr, iptr, ll, is, iee, in, iw, ir;
  int N = ((nete - nets + 1)*qsize + 64 - 1)/64;

  double *src_receive, *src_qmax, *src_qmin;
  int *src_getmap;

  for (n = 0; n < N; n++) {
    id_map = n*64 + id;
    if (id_map < (nete - nets + 1)*qsize) {
      ie = id_map/qsize;
      q = id_map%qsize;
      //if (id == 25) printf("%d\n", q);
      //if (id == 1) printf("id:%d,id_map:%d,ie:%d,q:%d\n", id, id_map, ie, q);
      src_getmap = gl_getmap + ie*max_neigh_edges;
      src_qmax = gl_qmax + ie*qsize*NLEV + q*NLEV;
      src_qmin = gl_qmin + ie*qsize*NLEV + q*NLEV;
      pe_get(src_getmap, getmap, max_neigh_edges*sizeof(int));
      pe_get(src_qmax, qmax, NLEV*sizeof(double));
      pe_get(src_qmin, qmin, NLEV*sizeof(double));
      dma_syn();

      is  = getmap[south - 1];
      iee = getmap[east - 1];
      in  = getmap[north - 1];
      iw  = getmap[west - 1];

      src_receive = gl_receive + iee + q*NLEV;
      pe_get(src_receive, receive_ie, NLEV*sizeof(double));
      dma_syn();
      src_receive = gl_receive + is + q*NLEV;
      pe_get(src_receive, receive_is, NLEV*sizeof(double));
      dma_syn();
      src_receive = gl_receive + iw + q*NLEV;
      pe_get(src_receive, receive_iw, NLEV*sizeof(double));
      dma_syn();
      src_receive = gl_receive + in + q*NLEV;
      pe_get(src_receive, receive_in, NLEV*sizeof(double));
      dma_syn();
      for (k = 0; k < NLEV; k++) {
        MIN5(qmin[k], receive_ie[k], receive_is[k], receive_iw[k], receive_in[k]);
      }
      src_receive = gl_receive + iee + q*NLEV + qsize*NLEV;
      pe_get(src_receive, receive_ie, NLEV*sizeof(double));
      dma_syn();
      src_receive = gl_receive + is + q*NLEV + qsize*NLEV;
      pe_get(src_receive, receive_is, NLEV*sizeof(double));
      dma_syn();
      src_receive = gl_receive + iw + q*NLEV + qsize*NLEV;
      pe_get(src_receive, receive_iw, NLEV*sizeof(double));
      dma_syn();
      src_receive = gl_receive + in + q*NLEV + qsize*NLEV;
      pe_get(src_receive, receive_in, NLEV*sizeof(double));
      dma_syn();
      for (k = 0; k < NLEV; k++) {
        MAX5(qmax[k], receive_ie[k], receive_is[k], receive_iw[k], receive_in[k]);
      }
      // ---------------------------- SWEST ----------------------------------//
      for (ll = swest - 1; ll < (swest + max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q*NLEV;
          pe_get(src_receive, receive, NLEV*sizeof(double));
          dma_syn();
          for (k = 0; k < NLEV; k++) {
            qmin[k] = MIN(qmin[k], receive[k]);
          }
          src_receive = gl_receive + getmap[ll] + q*NLEV + qsize*NLEV;
          pe_get(src_receive, receive, NLEV*sizeof(double));
          dma_syn();
          for (k = 0; k < NLEV; k++) {
            qmax[k] = MAX(qmax[k], receive[k]);
          }
        }
      }  // end loop ll
      // ---------------------------- SEAST ----------------------------------//
      for (ll = swest + max_corner_elem - 1; ll < (swest + 2*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q*NLEV;
          pe_get(src_receive, receive, NLEV*sizeof(double));
          dma_syn();
          for (k = 0; k < NLEV; k++) {
            qmin[k] = MIN(qmin[k], receive[k]);
          }
          src_receive = gl_receive + getmap[ll] + q*NLEV + qsize*NLEV;
          pe_get(src_receive, receive, NLEV*sizeof(double));
          dma_syn();
          for (k = 0; k < NLEV; k++) {
            qmax[k] = MAX(qmax[k], receive[k]);
          }
        }
      } // end loop ll
      // ----------------------- NEAST -------------------------------- //
      for (ll = swest + 3*max_corner_elem - 1; ll < (swest + 4*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q*NLEV;
          pe_get(src_receive, receive, NLEV*sizeof(double));
          dma_syn();
          for (k = 0; k < NLEV; k++) {
            qmin[k] = MIN(qmin[k], receive[k]);
          }
          src_receive = gl_receive + getmap[ll] + q*NLEV + qsize*NLEV;
          pe_get(src_receive, receive, NLEV*sizeof(double));
          dma_syn();
          for (k = 0; k < NLEV; k++) {
            qmax[k] = MAX(qmax[k], receive[k]);
          }
        }
      } // end loop ll
      // ----------------------- NWEST -------------------------------- //
      for (ll = swest + 2*max_corner_elem - 1; ll < (swest + 3*max_corner_elem - 1); ll++) {
        if (getmap[ll] != -1) {
          src_receive = gl_receive + getmap[ll] + q*NLEV;
          pe_get(src_receive, receive, NLEV*sizeof(double));
          dma_syn();
          for (k = 0; k < NLEV; k++) {
            qmin[k] = MIN(qmin[k], receive[k]);
          }
          src_receive = gl_receive + getmap[ll] + q*NLEV + qsize*NLEV;
          pe_get(src_receive, receive, NLEV*sizeof(double));
          dma_syn();
          for (k = 0; k < NLEV; k++) {
            qmax[k] = MAX(qmax[k], receive[k]);
          }
        }
      } // end loop ll
      for (k = 0; k < NLEV; k++) {
        qmin[k] = MAX(qmin[k], 0.0);
      }
      pe_put(src_qmin, qmin, NLEV*sizeof(double));
      pe_put(src_qmax, qmax, NLEV*sizeof(double));
      dma_syn();
    } // end if
  } // end loop n
}
