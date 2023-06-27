#include "Autotuning.hpp"
#include "catch.hpp"

#include <omp.h>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <algorithm>

double func01(double x) {
  return 2 * sin(pow(x, 2)) + 5 * pow(sin(x), 3) + pow(cos(pow(x, 2)), 4) +
         pow(x, 2) / 1000 + 10;
}

int find_pos(double element, double *vector, int size) {
  for (int i = 0; i < size; i++) {
    if (element == vector[i]) return i;
  }
  return -1;
}

TEST_CASE("INTEGRITY") {
  const int dim = 1;
  const int min = -100, max = 100;
  const int ignore = 0;
  const int n_opt = 4;
  const int n_iter = 100;
  int pos;
  int value;
  int count = 0;
  double time;
  int *values = new int[3 * (n_iter + n_opt) + 5];
  double *times = new double[3 * (n_iter + n_opt) + 5];

  Autotuning *at = new Autotuning(dim, min, max, ignore, n_opt, n_iter);

  wStart(at, value);
  time = omp_get_wtime();
  usleep(func01(double(value)) * 1000);
  values[count] = value;
  times[count++] = omp_get_wtime() - time;
  wEnd(at);

  pos = find_pos(*std::min_element(times, times + count - 1), times, count);  
  REQUIRE(value == values[pos]);

  at->reset(2);

  wStart(at, value);
  time = omp_get_wtime();
  usleep(func01(double(value)) * 1000);
  values[count] = value;
  times[count++] = omp_get_wtime() - time;
  wEnd(at);

  pos = find_pos(*std::min_element(times, times + count - 1), times, count);
  REQUIRE(value == values[pos]);

  at->reset(1);

  wStart(at, value);
  time = omp_get_wtime();
  usleep(func01(double(value)) * 1000);
  values[count] = value;
  times[count++] = omp_get_wtime() - time;
  wEnd(at);

  pos = find_pos(*std::min_element(times, times + count - 1), times, count);
  REQUIRE(value == values[pos]);

  at->reset(0);
  count = 0;

  wStart(at, value);
  time = omp_get_wtime();
  usleep(func01(double(value)) * 1000);
  values[count] = value;
  times[count++] = omp_get_wtime() - time;
  wEnd(at);

  pos = find_pos(*std::min_element(times, times + count - 1), times, count);
  REQUIRE(value == values[pos]);

  delete at;
}
