#pragma once

class Test {
  inline static std::vector<std::vector<int>> Points;
  inline static std::vector<double> Costs;

 public:
  Test() {}
  ~Test() {}

  static double function(int dim, int *p) {
    std::vector<int> _p;
    double cost = equation(p, dim);
    for (int i = 0; i < dim; i++) _p.push_back(p[i]);
    Points.push_back(_p);
    Costs.push_back(cost);
    return cost;
  }

  static std::vector<int> get_min_point() {
    double element = *std::min_element(Costs.begin(), Costs.end());
    for (int i = 0; i < (int)Points.size(); i++) {
      if (element == Costs[i]) return Points[i];
    }
    throw std::runtime_error("Element not found\n");
  }

  static void reset() {
    Costs.clear();
    Points.clear();
  }

 private:
  static double equation(int *p, int dim) {
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
};