#ifndef AUTOTUNING_H
#define AUTOTUNING_H

#include "CSA.hpp"

#define wStart(AT, POINT)     \
  while (!(AT)->finished()) { \
    (AT)->start(&(POINT));
#define wEnd(AT) \
  (AT)->end();   \
  }
#define autotuning_run(AT, __FUNCTION__, VALUE) \
  while (!(AT)->finished()) {                   \
    (AT)->start(&(VALUE));                      \
    (__FUNCTION__);                             \
    (AT)->end();                                \
  }

class Autotuning {
 private:
  // Global Variables
  bool finished_flag;  // Flag for Auto-tuning end executation
  int min;             // Minimum value of search interval
  int max;             // Maximum value of search interval
  int iteration;       // Iteration number
  int ignore;          // Numeber of ignore recive values
  double *cost;        // Cost value [f(point)]

  double t0, t1;   // Starting and Ending Times
  double runtime;  // Total time of a task

  // CSA Variables
  CSA *csa;  // Class with valiables specific of the 'csa->h'
  int iOpt;  // Iterator for optimizers reciver to csa()
             // double **optPoints;  // CSA optimizers points
             // int *allPoints;

  void rescale(double *out, int *in) const;
  void rescale(int *out, double *in) const;

 public:
  auto finished() const -> bool { return finished_flag; }
  void start(int *point);
  void end();
  void print();
  void reset(int level);

  Autotuning operator=(Autotuning) = delete;
  auto operator=(Autotuning &&) -> Autotuning& = delete;	
  Autotuning(const Autotuning &) = delete;
  Autotuning(Autotuning &&) = delete;  

  Autotuning() = default;
  Autotuning(int _dim, int _min, int _max, int _ignore, int _num_opt,
             int _numInt);
  ~Autotuning() {
    delete[] cost;
    delete csa;
  }
};

#endif