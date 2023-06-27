#ifndef _AUTOTUNING_
#define _AUTOTUNING_

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef C_S_A
#include "csa.h"
#define AT_wStart(at, po) \
  do {                    \
  AT_start(&at, &po)
#define AT_wEnd(at) \
  AT_end(&at);      \
  }                 \
  while (at.end == 0)
#endif

#ifndef OPT
#define OPT 5
// ERRO;
#endif
#ifndef ITE
#define ITE 300
// ERRO;
#endif

typedef struct {
  // Global Variables
  int min;        // Minimum value of search interval
  int max;        // Maximum value of search interval
  double *cost;   // Cost value [f(point)]
  int iteration;  // Iteration number
  int ignore;     // Numeber of ignore recive values
  int end;        // Flag for Auto-tuning end executation

  double t0, t1;  // Starting and Ending Times
  double T1;

  // CSA Variables
  CSA csa;             // Struct with valiables specific of the 'csa.h'
  int iOpt;            // Iterator for optimizers reciver to csa()
  double **optPoints;  // CSA optimizers points
  int *allPoints;

} Autotuning;

void AT_reset(Autotuning *at);
void AT_reset(Autotuning *at, int new_min, int new_max);
void AT_init(Autotuning *at, int _dim, int _min, int _max, int ignore,
             int _numOpt, int _numInt);  //, int **_points, int _num);
void AT_start(Autotuning *at, int *point);
void AT_end(Autotuning *at);
void AT_setPoint(Autotuning *at, int **points, int num);

void AT_finalize(Autotuning *at);
// void AT_bestReset(Autotuning *at);

#endif