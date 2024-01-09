#include <omp.h>  // omp_get_wtime

#include <cmath>  // round

template <typename Point>
void Autotuning::rescale(Point *out, const double *in) const {
  if constexpr (std::is_floating_point<Point>::value) {
    for (int i = 0; i < p_optimizer->getDimension(); i++) {
      out[i] = static_cast<Point>(round(((in[i] + 1.0) / 2.0) * (m_max - m_min) + m_min));
    }
  } else if (std::is_integral<Point>::value) {
    for (int i = 0; i < p_optimizer->getDimension(); i++) {
      out[i] = static_cast<Point>(((in[i] + 1.0) / 2.0) * (m_max - m_min) + m_min);
    }
  }
}

template <typename Point>
void Autotuning::exec(Point *point, double cost) {
  static_assert(std::is_arithmetic<Point>::value,
                "Autotuning only supports integer and floating-point types");
  // End execution
  if (!isEnd()) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
      // Itaration ignore test
      if (((++m_iter) % m_ignore) == 0) {
        double *d_point = p_optimizer->run(cost);
        rescale(point, d_point);
      }
    }
  }
}

template <typename Point>
void Autotuning::start(Point *point) {
  static_assert(std::is_arithmetic<Point>::value,
                "Autotuning only supports integer and floating-point types");

  // End execution
  if (!isEnd()) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
      // Itaration ignore test
      if (((++m_iter) % m_ignore) == 0) {
        double *d_point = p_optimizer->run(m_runtime);
        rescale(point, d_point);
      }

      m_t0 = omp_get_wtime();
    }
  }
}

template <typename Point, typename Func, typename... Args>
void Autotuning::entireExecRuntime(Func function, Point *point, Args... args) {
  static_assert(std::is_arithmetic<Point>::value,
                "Autotuning only supports integer and floating-point types");
  while (!isEnd()) {
    start(point);
    function(args..., point);
    end();
  }
}

template <typename Point, typename Func, typename... Args>
void Autotuning::entireExec(Func function, Point *point, Args... args) {
  static_assert(std::is_arithmetic<Point>::value,
                "Autotuning only supports integer and floating-point types");
  while (!isEnd()) {
    exec(point, m_cost);
    m_cost = function(args..., point);
  }
}

template <typename Point, typename Func, typename... Args>
void Autotuning::singleExecRuntime(Func function, Point *point, Args... args) {
  static_assert(std::is_arithmetic<Point>::value,
                "Autotuning only supports integer and floating-point types");
  start(point);
  function(args..., point);
  end();
}

template <typename Point, typename Func, typename... Args>
void Autotuning::singleExec(Func function, Point *point, Args... args) {
  static_assert(std::is_arithmetic<Point>::value,
                "Autotuning only supports integer and floating-point types");
  exec(point, m_cost);
  m_cost = function(args..., point);
}