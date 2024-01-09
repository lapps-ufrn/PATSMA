#include "autotuning.h"

#define reScaleToInt(in, out, min, max, size) \
  for (int ii = 0; ii < size; ii++)           \
    out[ii] = ((in[ii] + 1) / 2) * (max - min) + min;
#define reScaleToDouble(in, out, min, max, size) \
  for (int ii = 0; ii < size; ii++)              \
    out[ii] = ((1.0 * in[ii] - 1.0 * min) / (1.0 * max - 1.0 * min)) * 2 - 1;
#define isAlloc(v, s)                               \
  if (v == NULL) {                                  \
    printf("Failed allocating memory for " s "\n"); \
    exit(0);                                        \
  }
#define copy(a, b, s) \
  for (int ii = 0; ii < s; ii++) b[ii] = a[ii];

/*
 * at: in - Auto-tuning Global Variables
 */
void AT_reset(Autotuning *at) {
  csa_reset(at->csa);
  at->end = 0;
}

/*
 * at: in - Auto-tuning Global Variables
 * new_min: in - New Minimun Value
 * new_max: in - New Maximun Value
 */
void AT_reset(Autotuning *at, int new_min, int new_max) {
  csa_reset(at->csa);
  at->min = new_min != -1 ? new_min : at->min;
  at->max = new_max != -1 ? new_max : at->max;
  at->end = 0;
}

void AT_setPoint(Autotuning *at, int **points, int num) {
#ifdef _OPENMP
#pragma omp single nowait
#endif
  {
    for (int i = 0; i < num && i < at->csa.numOpt; i++) {
      reScaleToDouble(points[i], at->csa.opts[i].curSol, at->min, at->max,
                      at->csa.dim);
    }
  }
}

/*
 * at: in - Auto-tuning Global Variables
 * _dim: in - Cost Function Dimension
 * _min: in - Minimum value of search interval
 * _max: in - Maximum value of search interval
 * _ignore: in - Amount of ignore iterations, interlevead with one accpeted
 * _numOpt: in - Amount of optimizer
 */
void AT_init(Autotuning *at, int _dim, int _min, int _max, int _ignore,
             int _numOpt, int _numInt) {
#ifdef _OPENMP
#pragma omp single nowait
#endif
  {
    if (_dim < 1) {
      printf("Dimensional Value Invalid! Set _dim > 0.\n");
      exit(0);
    }
    if (_numOpt < 1) {
      printf("Optmizers Number Invalid! Set _numOpt > 0.\n");
      exit(0);
    }
    if (_ignore < 0) {
      printf("Ignore Value Invalid! Set _ignore > -1.\n");
      exit(0);
    }

    at->min = _min;
    at->max = _max;
    at->ignore = _ignore + 1;

    at->iOpt = 0;
    at->end = 0;
    at->iteration = 0;

    // at->point = (int *)malloc(sizeof(int) * _dim);
    // isAlloc(at->point,"at->point");

    at->cost = (double *)malloc(sizeof(double) * _numOpt);
    isAlloc(at->cost, "at->cost");

    at->optPoints = (double **)malloc(sizeof(double *) * _numOpt);
    isAlloc(at->optPoints, "at->optPoints");

    for (int i = 0; i < _numOpt; i++) {
      at->optPoints[i] = (double *)malloc(sizeof(double) * _dim);
      isAlloc(at->optPoints[i], "at->optPoints[i]");
    }

    at->optPoints = csa_init(&at->csa, _numOpt, _dim, (_numInt / at->ignore));
  }
}

/*
 * at: in - Auto-tuning Global Variables
 * _cost: in - Cost Value [f(point)]
 * return: out - Caculated point
 */
void AT_start(Autotuning *at, int *point) {
  // End execution
  if (at->end) return;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
  {
    // Itaration ignore test
    if (((++at->iteration) % at->ignore) == 0) {
      // Cost values save
      if (at->iOpt > 0 && at->iOpt <= at->csa.numOpt)
        at->cost[at->iOpt - 1] = at->T1;

      // Return to aplication the next solution
      if (at->iOpt < at->csa.numOpt) {
        reScaleToInt(at->optPoints[at->iOpt], point, at->min, at->max,
                     at->csa.dim);
        at->iOpt += 1;
      } else {
        at->optPoints = csa_exec(&at->csa, at->cost);

        if (at->csa.step == END)  // End execution, return best solution
        {
          at->end = 1;
          reScaleToInt(at->csa.bestSol, point, at->min, at->max, at->csa.dim);
        } else  // Return first optimazate
        {
          reScaleToInt(at->optPoints[0], point, at->min, at->max, at->csa.dim);
          at->iOpt = 1;
        }
      }
    }
  }

#ifdef _OPENMP
#pragma omp single nowait
  { at->t0 = omp_get_wtime(); }
#endif
}

/*
 * at: in - Auto-tuning Global Variables
 * _cost: in - Cost Value [f(point)]
 * return: out - Caculated point
 */
void AT_end(Autotuning *at) {
  if (at->end) return;

#ifdef _OPENMP
#pragma omp single nowait
  {
    at->t1 = omp_get_wtime();
    at->T1 = (at->t1 - at->t0);
  }
#endif
}

void AT_finalize(Autotuning *at) {
  // free(at->point);
  free(at->cost);

  for (int i = 0; i < at->csa.numOpt; i++) {
    if (at->optPoints[i]) {
      free(at->optPoints[i]);
    }
  }
  if (at->optPoints) {
    free(at->optPoints);
  }

  csa_finalize(&at->csa);
}