#include "Autotuning.hpp"

#include <omp.h>

#include <cmath>     // round
#include <iostream>  // cout, endl
#include <limits>    // numeric_limits

#define DBL_MAX std::numeric_limits<double>::max()
#define isAlloc(v, s)                                               \
  if ((v) == NULL) {                                                \
    std::cout << "Failed allocating memory for " s "" << std::endl; \
    exit(0);                                                        \
  }

inline void Autotuning::rescale(double *out, int *in) const {
  for (int i = 0; i < csa->dim; i++) {
    out[i] = (((double)in[i] + 1.0) / 2.0) * ((double)max - (double)min) +
             (double)min;
  }
}

inline void Autotuning::rescale(int *out, double *in) const {
  for (int i = 0; i < csa->dim; i++) {
    out[i] = (int)round(((in[i] + 1.0) / 2.0) * ((double)max - (double)min) +
                        (double)min);
  }
}

/*
 * at: in - Auto-tuning Global Variables
 * _dim: in - Cost Function Dimension
 * _min: in - Minimum value of search interval
 * _max: in - Maximum value of search interval
 * _ignore: in - Amount of ignore iterations, interlevead with one accpeted
 * _num_opt: in - Amount of optimizer
 */

Autotuning::Autotuning(int _dim, int _min, int _max, int _ignore, int _num_opt,
                       int _numInt)
    : finished_flag(false),
      min(_min),
      max(_max),
      iteration(0),
      t0(0),
      t1(0),
      runtime(0),
      iOpt(0) {
#ifdef _OPENMP
#pragma omp single nowait
#endif
  {
    if (_dim < 1) {
      std::cout << "Dimensional Value Invalid! Set _dim > 0." << std::endl;
      exit(0);
    }
    if (_num_opt < 1) {
      std::cout << "Optmizers Number Invalid! Set _num_opt > 0." << std::endl;
      exit(0);
    }
    if (_ignore < 0) {
      std::cout << "Ignore Value Invalid! Set _ignore > -1." << std::endl;
      exit(0);
    }

    ignore = _ignore + 1;

    cost = new double[_num_opt];
    isAlloc(cost, "cost");

    csa = new CSA(_num_opt, _dim, (_numInt / ignore));
#ifdef VERBOSE
    print();
#endif
  }
}

void Autotuning::reset(int level) {
#ifdef _OPENMP
#pragma omp single
#endif
  {
    iOpt = 0;
    finished_flag = false;
    iteration = 0;

    csa->reset(level);
  }
}

/*
 * at: in - Auto-tuning Global Variables
 * _cost: in - Cost Value [f(point)]
 * return: out - Caculated point
 */
void Autotuning::start(int *point) {
  // End execution
  if (!finished_flag) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
#ifdef _OPENMP
      t1 = omp_get_wtime();
#endif
      runtime = (t1 - t0);

      // Itaration ignore test
      if (((++iteration) % ignore) == 0) {
        // Cost values save
        if (iOpt > 0 && iOpt <= csa->num_opt) {
          cost[iOpt - 1] = runtime;
        }
        // Return to aplication the next solution
        if (iOpt < csa->num_opt) {
          // reScaleToInt(point, csa->solution[iOpt], min, max, csa->dim);
          rescale(point, csa->solution[iOpt]);
          iOpt += 1;
        } else {
          csa->partial_exec(cost);
          if (csa->step == END) {  // End execution, return best solution
            finished_flag = true;
            // reScaleToInt(point, csa->best_sol, min, max, csa->dim);
            rescale(point, csa->best_sol);
          } else {  // Return first optimazate
            // reScaleToInt(point, csa->solution[0], min, max, csa->dim);
            rescale(point, csa->solution[0]);
            iOpt = 1;
          }
        }
      }
#ifdef _OPENMP
      t0 = omp_get_wtime();
#endif
    }
  }
}

/*
 * at: in - Auto-tuning Global Variables
 * _cost: in - Cost Value [f(point)]
 * return: out - Caculated point
 */
void Autotuning::end() {
  if (!finished_flag) {
#ifdef _OPENMP
#pragma omp single
#endif
    {
#ifdef _OPENMP
      t1 = omp_get_wtime();
#endif
      runtime = (t1 - t0);
    }
  }
}

void Autotuning::print() {
  int *aux = new int[csa->dim];
  std::cout << "------------------- Autotuning Parameters -------------------"
            << std::endl;
  std::cout << "Dim: " << csa->dim << "\tNOpt: " << csa->num_opt
            << "\tNIgn: " << ignore
            << "\tNEval: " << (csa->max_iter * csa->num_opt * ignore)
            << std::endl;
  std::cout << "Min: " << min << "\tMax: " << max << std::endl;
  std::cout << "{";
  for (int i = 0; i < csa->num_opt; i++) {
    std::cout << "[" << i << ", ";
    std::cout << "(";
    // reScaleToInt(aux, csa->solution[i], min, max, csa->dim);
    rescale(aux, csa->solution[i]);
    for (int j = 0; j < csa->dim; j++) {
      std::cout << aux[j] << ((j < csa->dim - 1) ? "," : "");
    }
    std::cout << ")";
    std::cout << "]" << ((i < csa->num_opt - 1) ? "\t" : "");
  }
  std::cout << "}" << std::endl;
  std::cout << "-------------------------------------------------------------"
            << std::endl;
  delete[] aux;
}
