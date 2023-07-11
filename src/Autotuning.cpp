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

/// @brief Rescale point from int [min,max] to double [-1,1]
/// @param out Int point
/// @param in Double point
inline void Autotuning::rescale(double *out, int *in) const {
  for (int i = 0; i < csa->dim; i++) {
    out[i] = (((double)in[i] + 1.0) / 2.0) * ((double)max - (double)min) +
             (double)min;
  }
}

/// @brief Rescale point from double [-1,1] to int [min,max]
/// @param out Double point
/// @param in Int point
inline void Autotuning::rescale(int *out, double *in) const {
  for (int i = 0; i < csa->dim; i++) {
    out[i] = (int)round(((in[i] + 1.0) / 2.0) * ((double)max - (double)min) +
                        (double)min);
  }
}

/// @brief Auto-tuning constructor
/// @param _dim Cost Function Dimension
/// @param _min Minimum value of search interval
/// @param _max Maximum value of search interval
/// @param _ignore Amount of ignore iterations, interlevead with one accpeted
/// @param _num_opt Amount of optimizer
/// @param _numInt Number of iterations
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

/// @brief Reset the autotuning and CSA
/// @param level Select the level of reseting
///    level 2 - Reset the number of iteretions
///    level 1 - Reset the points and the temperatures (plus the previous ones)
///    level 0 - Remove the best solution (plus the previous ones)
void Autotuning::reset(int level) {
#ifdef _OPENMP
#pragma omp single
#endif
  {
    iOpt = 0;
    iteration = 0;
    finished_flag = false;

    csa->reset(level);
  }
}

/// @brief Goes before the cost function
/// @param point The point to be analyzed
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
            rescale(point, csa->best_sol);
          } else {  // Return first optimazate
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

/// @brief Goes after the cost function
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

/// @brief Printe class basic information
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
