#include "Autotuning.hpp"

#include <cmath>     // round
#include <iostream>  // cout, endl
#ifdef _OPENMP
#include <omp.h>
#else
#include <chrono>
#endif

#define AT_OPT_CSA 0x11
#define AT_OPT_NM 0x22

/// @brief Rescale point from int [min,max] to double [-1,1]
/// @param out Int point
/// @param in Double point
inline void Autotuning::rescale(double *out, int *in) const {
  for (int i = 0; i < optimizer->getDimension(); i++) {
    out[i] = (((double)in[i] + 1.0) / 2.0) * ((double)max - (double)min) +
             (double)min;
  }
}

/// @brief Rescale point from double [-1,1] to int [min,max]
/// @param out Double point
/// @param in Int point
inline void Autotuning::rescale(int *out, double *in) const {
  for (int i = 0; i < optimizer->getDimension(); i++) {
    out[i] = (int)round(((in[i] + 1.0) / 2.0) * ((double)max - (double)min) +
                        (double)min);
  }
}

/// @brief Auto-tuning constructor using CSA default optimizer
/// @param dim Cost Function Dimension
/// @param _min Minimum value of search interval
/// @param _max Maximum value of search interval
/// @param _ignore Amount of ignore iterations, interlevead with one accpeted
/// @param num_opt Amount of optimizer
/// @param max_iter Maximum number of iterations
Autotuning::Autotuning(int dim, int _min, int _max, int _ignore, int num_opt,
                       int max_iter)
    : Autotuning(_min, _max, _ignore, new CSA(num_opt, dim, max_iter)) {}

Autotuning::Autotuning(int _min, int _max, int _ignore,
                       NumericalOptimizer *_optimizer)
    : min(_min), max(_max), iteration(0), runtime(0) {
  if (_ignore < 0) {
    throw std::invalid_argument("Ignore Value Invalid! Set _ignore >= 0.");
  }

  ignore = _ignore + 1;
  optimizer = _optimizer;

#ifdef VERBOSE
  print();
#endif
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
    iteration = 0;
    optimizer->reset(level);
  }
}

/// @brief Goes before the cost function
/// @param point The point to be analyzed
void Autotuning::run(int *point, double _cost) {
  // End execution
  if (!optimizer->isEnd()) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
      // Itaration ignore test
      double *d_point;
      if (((++iteration) % ignore) == 0) {
        d_point = optimizer->run(_cost);
        rescale(point, d_point);
      }
    }
  }
}

/// @brief Goes before the cost function
/// @param point The point to be analyzed
void Autotuning::start(int *point) {
  // End execution
  if (!optimizer->isEnd()) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
      double *d_point;

      // Itaration ignore test
      if (((++iteration) % ignore) == 0) {
        d_point = optimizer->run(runtime);
        rescale(point, d_point);
      }
      // printf("%i ", point[0]);

      // t0 = std::chrono::high_resolution_clock::now();
      t0 = omp_get_wtime();
    }
  }
}

/// @brief Goes after the cost function
void Autotuning::end() {
  if (!optimizer->isEnd()) {
#ifdef _OPENMP
#pragma omp single
#endif
    {
      t1 = omp_get_wtime();
      // t1 = std::chrono::high_resolution_clock::now();
      // std::chrono::duration<double> elapsed = t1 - t0;
      // std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
      runtime = t0 - t1;
      printf("%lf %lf %lf\n", runtime, t0, t1);
    }
  }
}

/// @brief Printe class basic information
void Autotuning::print() {
  std::cout << "------------------- Autotuning Parameters -------------------"
            << std::endl;
  std::cout << "\tNIgn: " << ignore;
  std::cout << "Min: " << min << "\tMax: " << max << std::endl;
  optimizer->print();
  std::cout << "-------------------------------------------------------------"
            << std::endl;
}
