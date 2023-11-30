#include <omp.h>  // omp_get_wtime

#include <cmath>  // round

template <typename POINT>
void Autotuning::rescale(POINT *out, double *in) const {
  if (std::is_floating_point<POINT>::value) {
    for (int i = 0; i < p_optimizer->getDimension(); i++) {
      out[i] = static_cast<POINT>(
          round(((in[i] + 1.0) / 2.0) * (m_max - m_min) + m_min));
    }
  } else if (std::is_integral<POINT>::value) {
    for (int i = 0; i < p_optimizer->getDimension(); i++) {
      out[i] =
          static_cast<POINT>(((in[i] + 1.0) / 2.0) * (m_max - m_min) + m_min);
    }
  }
}

template <typename POINT, typename COST>
void Autotuning::run(POINT *point, COST _cost) {
  // End execution
  if (!p_optimizer->isEnd()) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
      // Itaration ignore test
      if (((++m_iter) % m_ignore) == 0) {
        double *d_point = p_optimizer->run(_cost);
        rescale<POINT>(point, d_point);
      }
    }
  }
}

template <typename POINT = int>
void Autotuning::start(POINT *point) {
  // End execution
  if (!p_optimizer->isEnd()) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
      // Itaration ignore test
      if (((++m_iter) % m_ignore) == 0) {
        double *d_point = p_optimizer->run(m_runtime);
        rescale<POINT>(point, d_point);
      }

      m_t0 = omp_get_wtime();
    }
  }
}

template <typename POINT = int, typename Func, typename... Args>
POINT *Autotuning::execOfflineRuntime(Func function, Args... args) {
  int *point = new int[p_optimizer->getDimension()];
  while (!isEnd()) {
    start(point);
    function(point, args...);
    end();
  }
  return point;
}

template <typename POINT = int, typename COST = double, typename Func,
          typename... Args>
POINT *Autotuning::execOffline(Func function, Args... args) {
  COST cost;
  POINT *point = new POINT[p_optimizer->getDimension()];
  while (!isEnd()) {
    run(point, cost);
    cost = function(point, args...);
  }
  return point;
}

template <typename POINT = int, typename Func, typename... Args>
void Autotuning::execOnlineRuntime(Func function, POINT *point, Args... args) {
  start(point);
  function(point, args...);
  end();
}

template <typename POINT = int, typename COST = double, typename Func,
          typename... Args>
void Autotuning::execOnline(Func function, POINT *point, Args... args) {
  COST cost;
  run(point, cost);
  cost = function(point, args...);
}