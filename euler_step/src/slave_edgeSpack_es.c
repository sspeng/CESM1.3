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

typedef struct {
  double *qmax, *qmin, *buf;
  int *putmap;
  int nets, nete, qsize;
} param_t;

void slave_edgespack_es_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  dma_init();
  ldm_alloc_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();
  double *gl_qmax = param_d.qmax;
  double *gl_qmin = param_d.qmin;
  double *gl_buf = param_d.buf;
  int *gl_putmap = param_d.putmap;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int qsize = param_d.qsize;

  double qmax[NLEV];
  double qmin[NLEV];
  int putmap[max_neigh_edges];

  int n, ie, q, id_map, kptr, iptr, ll, is, iee, in, iw, ir;
  int N = ((nete - nets + 1)*qsize + 64 - 1)/64;

  double *src_buf, *src_qmax, *src_qmin;
  int *src_putmap;

  for (n = 0; n < N; n++) {
    id_map = n*64 + id;
    if (id_map < (nete - nets + 1)*qsize) {
      ie = id_map/qsize;
      q = id_map%qsize;
      //if (id == 1) printf("id:%d,id_map:%d,ie:%d,q:%d\n", id, id_map, ie, q);
      src_putmap = gl_putmap + ie*max_neigh_edges;
      src_qmax = gl_qmax + ie*qsize*NLEV + q*NLEV;
      src_qmin = gl_qmin + ie*qsize*NLEV + q*NLEV;
      pe_get(src_putmap, putmap, max_neigh_edges*sizeof(int));
      pe_get(src_qmax, qmax, NLEV*sizeof(double));
      pe_get(src_qmin, qmin, NLEV*sizeof(double));
      dma_syn();
      is = putmap[south - 1];
      iee = putmap[east - 1];
      in = putmap[north - 1];
      iw = putmap[west - 1];
      //for (k = 0; k < NLEV; k++) {
      //  buf[k] = qmin[k];
      //}
      src_buf = gl_buf + iee + q*NLEV;
      pe_put(src_buf, qmin, NLEV*sizeof(double)); // East for qmin
      dma_syn();
      src_buf = gl_buf + iee + q*NLEV + qsize*NLEV;
      pe_put(src_buf, qmax, NLEV*sizeof(double)); // East for qmax
      dma_syn();
      src_buf = gl_buf + is + q*NLEV;
      pe_put(src_buf, qmin, NLEV*sizeof(double)); // South for qmin
      dma_syn();
      src_buf = gl_buf + is + q*NLEV + qsize*NLEV;
      pe_put(src_buf, qmax, NLEV*sizeof(double)); // South for qmax
      dma_syn();
      src_buf = gl_buf + in + q*NLEV;
      pe_put(src_buf, qmin, NLEV*sizeof(double)); // North for qmin
      dma_syn();
      src_buf = gl_buf + in + q*NLEV + qsize*NLEV;
      pe_put(src_buf, qmax, NLEV*sizeof(double)); // North for qmax
      dma_syn();
      src_buf = gl_buf + iw + q*NLEV;
      pe_put(src_buf, qmin, NLEV*sizeof(double)); // West for qmin
      dma_syn();
      src_buf = gl_buf + iw + q*NLEV + qsize*NLEV;
      pe_put(src_buf, qmax, NLEV*sizeof(double)); // West for qmax
      dma_syn();

      // ---------------------------- SWEST ----------------------------------//
      for (ll = swest - 1; ll < (swest + max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          src_buf = gl_buf + putmap[ll] + q*NLEV;
          pe_put(src_buf, qmin, NLEV*sizeof(double));
          dma_syn();
          src_buf = gl_buf + putmap[ll] + q*NLEV + qsize*NLEV;
          pe_put(src_buf, qmax, NLEV*sizeof(double));
          dma_syn();
        }
      }  // end loop ll
      // ---------------------------- SEAST ----------------------------------//
      for (ll = swest + max_corner_elem - 1; ll < (swest + 2*max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          src_buf = gl_buf + putmap[ll] + q*NLEV;
          pe_put(src_buf, qmin, NLEV*sizeof(double));
          dma_syn();
          src_buf = gl_buf + putmap[ll] + q*NLEV + qsize*NLEV;
          pe_put(src_buf, qmax, NLEV*sizeof(double));
          dma_syn();
        }
      }
      // ----------------------- NEAST -------------------------------- //
      for (ll = swest + 3*max_corner_elem - 1; ll < (swest + 4*max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          src_buf = gl_buf + putmap[ll] + q*NLEV;
          pe_put(src_buf, qmin, NLEV*sizeof(double));
          dma_syn();
          src_buf = gl_buf + putmap[ll] + q*NLEV + qsize*NLEV;
          pe_put(src_buf, qmax, NLEV*sizeof(double));
          dma_syn();
        }
      }
      // ----------------------- NWEST -------------------------------- //
      for (ll = swest + 2*max_corner_elem - 1; ll < (swest + 3*max_corner_elem - 1); ll++) {
        if (putmap[ll] != -1) {
          src_buf = gl_buf + putmap[ll] + q*NLEV;
          pe_put(src_buf, qmin, NLEV*sizeof(double));
          dma_syn();
          src_buf = gl_buf + putmap[ll] + q*NLEV + qsize*NLEV;
          pe_put(src_buf, qmax, NLEV*sizeof(double));
          dma_syn();
        }
      }
    }
  }
}
