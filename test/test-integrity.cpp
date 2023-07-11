#include <math.h>
#include <omp.h>
#include <unistd.h>

#include <algorithm>
#include <catch2/catch.hpp>
#include <iostream>

#include "Autotuning.hpp"

#define REQUIRE_MESSAGE(cond, msg) \
  do {                             \
    INFO(msg);                     \
    REQUIRE(cond);                 \
  } while ((void)0, 0)

class Test {
  int count;
  double time;
  int *values;
  double *times;

 public:
  Test(int n_iter, int n_opt) : count(0) {
    values = new int[3 * (n_iter + n_opt) + 5];
    times = new double[3 * (n_iter + n_opt) + 5];
  }
  ~Test() {
    delete values;
    delete times;
  }

  void test_function_01(double x) {
    time = omp_get_wtime();
    usleep(function_01(double(x)) * 1000);
    times[count] = omp_get_wtime() - time;
    values[count] = x;
    count++;
  }

  int get_min_value() {
    double element = *std::min_element(times, times + count - 1);
    for (int i = 0; i < count; i++) {
      if (element == times[i]) return values[i];
    }
    throw std::runtime_error("Element not found\n");
  }

  void zero_count() { count = 0; }

 private:
  double function_01(double x) {
    return 2 * sin(pow(x, 2)) + 5 * pow(sin(x), 3) + pow(cos(pow(x, 2)), 4) +
           pow(x, 2) / 1000 + 10;
  }
};

TEST_CASE("INTEGRITY") {
  const int dim = 1;
  const int min = -100, max = 100;
  const int ignore = 0;
  const int n_opt = 4;
  const int n_iter = 200;
  int value;

  Autotuning *at = new Autotuning(dim, min, max, ignore, n_opt, n_iter);
  Test *t = new Test(n_iter, n_opt);

  AUTOTUNING_RUN(at, value) { t->test_function_01(value); }
  REQUIRE_MESSAGE(value == t->get_min_value(), "Base Case");

  at->reset(2);

  AUTOTUNING_RUN(at, value) { t->test_function_01(value); }
  REQUIRE_MESSAGE(value == t->get_min_value(), "Reset 2 Case");

  at->reset(1);

  AUTOTUNING_RUN(at, value) { t->test_function_01(value); }
  REQUIRE_MESSAGE(value == t->get_min_value(), "Reset 1 Case");

  at->reset(0);
  t->zero_count();

  AUTOTUNING_RUN(at, value) { t->test_function_01(value); }
  REQUIRE_MESSAGE(value == t->get_min_value(), "Reset 0 Case");

  delete at;
  delete t;
}
