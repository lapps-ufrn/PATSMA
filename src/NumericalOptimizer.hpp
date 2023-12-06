#pragma once

class NumericalOptimizer {
 public:
  NumericalOptimizer() = default;
  virtual ~NumericalOptimizer() = default;

  // Delete copy constructor and copy assignment operator
  NumericalOptimizer(const NumericalOptimizer& other) = delete;
  NumericalOptimizer& operator=(const NumericalOptimizer& other) = delete;

  // Delete move constructor and move assignment operator
  NumericalOptimizer(NumericalOptimizer&& other) noexcept = delete;
  NumericalOptimizer& operator=(NumericalOptimizer&& other) noexcept = delete;

  virtual double* run(double cost) = 0;
  virtual int getNumPoints() const = 0;
  virtual int getDimension() const = 0;
  virtual bool isEnd() const = 0;
  virtual void reset(int level){};
  virtual void print() const {}
};

#include "CSA.hpp"
#include "NelderMead.hpp"