#pragma once

class NumericalOptimizer {
 public:
  NumericalOptimizer(){};
  virtual ~NumericalOptimizer() {}

  virtual double* run(double cost) = 0;
  virtual int getNumPoints() const = 0;
  virtual int getDimension() const = 0;
  virtual bool isEnd() const = 0;
  virtual void reset(int level){};
  virtual void print() const {}
};

#include "CSA.hpp"
#include "NelderMead.hpp"