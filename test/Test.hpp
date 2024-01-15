#pragma once

#include <omp.h>

class Test {
  inline static std::vector<std::vector<int>> Points;
  inline static std::vector<double> Costs;

 public:
  static void sleep_function(int *point);

  static double function(int dim, int *p);

  static std::vector<int> get_min_point();

  static void reset();

 private:
  static double equation(int *p, int dim);
};

void Test::sleep_function(int *point) {
  double t1 = omp_get_wtime();
  std::vector<int> _p{*point};
  usleep((*point + 100) * 10000);
  Points.push_back(_p);
  Costs.push_back(omp_get_wtime() - t1);
}
double Test::function(int dim, int *p) {
  std::vector<int> _p;
  double cost = equation(p, dim);
  for (int i = 0; i < dim; i++) _p.push_back(p[i]);
  Points.push_back(_p);
  Costs.push_back(cost);
  return cost;
}

std::vector<int> Test::get_min_point() {
  double element = *std::min_element(Costs.begin(), Costs.end());
  for (int i = 0; i < (int)Points.size(); i++) {
    if (element == Costs[i]) {
      return Points[i];
    }
  }

  throw std::runtime_error("Element not found\n");
}

void Test::reset() {
  Costs.clear();
  Points.clear();
}

//! PRIVATE Methods !//
double Test::equation(int *p, int dim) {
  switch (dim) {
    case 1:
      return 2 * sin(pow(p[0], 2)) + 5 * pow(sin(p[0]), 3) + pow(cos(pow(p[0], 2)), 4) +
             pow(p[0], 2) / 1000 + 10;
    case 2:
      return pow(1 - p[0], 2) + 100 * pow(p[1] - p[0] * p[0], 2);
    case 3:
      return pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2);
    case 4:
      return pow(p[0] - 2, 2) + pow(p[1] + 3, 2) + pow(p[2] - 1, 2) + pow(p[3], 2);
    default:
      throw std::runtime_error("There is no function with dimension " + std::to_string(dim));
  }
}