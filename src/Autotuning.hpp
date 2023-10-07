#ifndef AUTOTUNING_H
#define AUTOTUNING_H

#include <chrono>  // time_point, now, duration

#include "NumericalOptimizer.hpp"

typedef void (*FunctionPointer)(int, ...);

class Autotuning {
  int min;        // Minimum value of search interval
  int max;        // Maximum value of search interval
  int iteration;  // Iteration number
  int ignore;     // Numeber of ignore recive values
  double cost;

  // std::chrono::high_resolution_clock::time_point t0,
  //     t1;          // Starting and Ending Times
  double t0, t1;
  double runtime;  // Total time of a task

  NumericalOptimizer *optimizer;

  void rescale(double *out, int *in) const;
  void rescale(int *out, double *in) const;

 public:
  bool isEnd() const { return optimizer->isEnd(); }
  void run(int *point, double _cost);
  void start(int *point);
  void end();
  void print();
  void reset(int level);

  template <typename Func, typename... Args>
  int *execOfflineRuntime(Func function, Args... args) {
    int *point = new int[optimizer->getDimension()];
    while (!isEnd()) {
      start(point);
      function(point, args...);
      end();
    }
    return point;
  }

  template <typename Func, typename... Args>
  int *execOffline(Func function, Args... args) {
    int *point = new int[optimizer->getDimension()];
    while (!isEnd()) {
      run(point, cost);
      cost = function(point, args...);
    }
    return point;
  }

  template <typename Func, typename... Args>
  void execOnlineRuntime(Func function, int *point, Args... args) {
    start(point);
    function(point, args...);
    end();
  }

  template <typename Func, typename... Args>
  void execOnline(Func function, int *point, Args... args) {
    run(point, cost);
    cost = function(point, args...);
  }

  auto operator=(Autotuning) -> Autotuning = delete;
  auto operator=(Autotuning &&) -> Autotuning & = delete;
  Autotuning(const Autotuning &) = delete;
  Autotuning(Autotuning &&) = delete;

  Autotuning() = default;
  Autotuning(int dim, int _min, int _max, int _ignore, int num_opt,
             int max_iter);
  Autotuning(int _min, int _max, int _ignore, NumericalOptimizer *_optimizer);
  ~Autotuning() { delete optimizer; }
};
#endif