#include "Autotuning.hpp"

#include <omp.h>

#include <iostream>  // cout, endl

//
// inline  void
// Autotuning::rescale(double *out, Point *in) const {
//   for (int i = 0; i < p_optimizer->getDimension(); i++) {
//     out[i] = (((double)in[i] + 1.0) / 2.0) * ((double)m_max -
//     (double)m_min)
//     +
//              (double)m_min;
//   }
// }

Autotuning::Autotuning(int dim, double min, double max, int ignore, int num_opt, int max_iter)
    : Autotuning(min, max, ignore, new CSA(num_opt, dim, max_iter)) {}

Autotuning::Autotuning(double min, double max, int ignore, NumericalOptimizer *optimizer)
    : m_min(min),
      m_max(max),
      m_iter(0),
      m_ignore(ignore + 1),
      m_cost(0),
      m_t0(0),
      m_t1(0),
      m_runtime(0),
      p_optimizer(optimizer) {

  if (ignore < 0) {
    throw std::invalid_argument("Ignore Value Invalid! Set _ignore >= 0.");
  }

#ifdef VERBOSE
  print();
#endif
}

Autotuning::~Autotuning() {
  delete p_optimizer;
}

void Autotuning::reset(int level) {
#ifdef _OPENMP
#pragma omp single
#endif
  {
    m_iter = 0;
    p_optimizer->reset(level);
  }
}

void Autotuning::end() {
  if (!p_optimizer->isEnd()) {
#ifdef _OPENMP
#pragma omp single
#endif
    {
      m_t1 = omp_get_wtime();
      m_runtime = m_t0 - m_t1;
    }
  }
}

void Autotuning::print() const {
  std::cout << "------------------- Autotuning Parameters -------------------" << std::endl;
  std::cout << "\tNIgn: " << m_ignore;
  std::cout << "Min: " << m_min << "\tMax: " << m_max << std::endl;
  p_optimizer->print();
  std::cout << "-------------------------------------------------------------" << std::endl;
}
