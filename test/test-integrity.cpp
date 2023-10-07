#include <math.h>
#include <omp.h>
#include <unistd.h>

#include <catch2/catch.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "Autotuning.hpp"

#define REQUIRE_MESSAGE(cond, msg) \
  do {                             \
    INFO(msg);                     \
    REQUIRE(cond);                 \
  } while ((void)0, 0)

class Test {
  inline static std::vector<std::vector<int>> Points;
  inline static std::vector<double> Costs;

 public:
  Test() {}
  ~Test() {}

  static void time_function(int *p) {
    double time = omp_get_wtime();
    usleep(equation(p, 1) * 1000);
    Costs.push_back(omp_get_wtime() - time);
    Points.push_back({p[0]});
  }

  static double function(int *p, int dim) {
    std::vector<int> _p;
    double cost = equation(p, dim);
    for (int i = 0; i < dim; i++) _p.push_back(p[i]);
    Points.push_back(_p);
    Costs.push_back(cost);
    return cost;
  }

  static std::vector<int> get_min_point() {
    double element = *std::min_element(Costs.begin(), Costs.end());
    for (int i = 0; i < Points.size(); i++) {
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
        return 2 * sin(pow(p[0], 2)) + 5 * pow(sin(p[0]), 3) +
               pow(cos(pow(p[0], 2)), 4) + pow(p[0], 2) / 1000 + 10;
      case 2:
        return pow(1 - p[0], 2) + 100 * pow(p[1] - p[0] * p[0], 2);
      case 3:
        return pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2);
      case 4:
        return pow(p[0] - 2, 2) + pow(p[1] + 3, 2) + pow(p[2] - 1, 2) +
               pow(p[3], 2);
      default:
        throw std::runtime_error("There is no function with dimension " +
                                 std::to_string(dim));
    }
  }
};

TEST_CASE("CSA") {
  int dim = 1;
  int min = -100, max = 100;
  int ignore = 0;
  int n_opt = 4;
  int n_iter = 20;

  SECTION("Test exception throwing") {
    Autotuning *at;
    REQUIRE_THROWS_AS(at = new Autotuning(min, max, -1, nullptr),
                      std::invalid_argument);

    CSA *csa;
    REQUIRE_THROWS_AS(csa = new CSA(0, dim, n_iter), std::invalid_argument);
    REQUIRE_THROWS_AS(csa = new CSA(n_opt, 0, n_iter), std::invalid_argument);
    REQUIRE_THROWS_AS(csa = new CSA(n_opt, dim, 0), std::invalid_argument);
  }

  SECTION("Reset Cases") {
    Test::reset();
    std::vector<int> reset_levels = {2, 1, 0};
    int point;

    Autotuning *at = new Autotuning(dim, min, max, ignore, n_opt, n_iter);

    point = at->execOffline(Test::function, dim)[0];
    REQUIRE_MESSAGE(point == Test::get_min_point()[0], "Base Case");

    for (auto &&level : reset_levels) {
      at->reset(level);
      if (level == 0) Test::reset();

      point = at->execOffline(Test::function, dim)[0];
      REQUIRE_MESSAGE(point == Test::get_min_point()[0],
                      "Reset level " + std::to_string(level));
    }

    delete at;
  }

  SECTION("MultiDimensional") {
    Test::reset();
    std::vector<int> num_iters = {20, 200, 2000};
    std::vector<int> dimensions = {1, 2, 3, 4};

    for (auto &&n_iter : num_iters) {
      for (auto &&dim : dimensions) {
        SECTION("Test with Dim = " + std::to_string(dim) +
                " and n_iter = " + std::to_string(n_iter)) {
          Test::reset();
          Autotuning *at = new Autotuning(dim, min, max, ignore, n_opt, n_iter);
          int *point = new int[dim];

          point = at->execOffline(Test::function, dim);

          std::vector<int> min_point = Test::get_min_point();
          for (int j = 0; j < dim; j++) {
            REQUIRE(point[j] == min_point[j]);
          }

          delete[] point;
          delete at;
        }
      }
    }
  }
}

TEST_CASE("Nelder-Mead") {
  int dim = 1;
  int min = -100, max = 100;
  int ignore = 0;
  double e = 10;

  SECTION("Test exception throwing") {
    Autotuning *at;
    REQUIRE_THROWS_AS(at = new Autotuning(min, max, -1, nullptr),
                      std::invalid_argument);

    NelderMead *nm;
    REQUIRE_THROWS_AS(nm = new NelderMead(0, e), std::invalid_argument);
    REQUIRE_THROWS_AS(nm = new NelderMead(dim, -.1), std::invalid_argument);
  }

  SECTION("Reset Cases") {
    Test::reset();

    int *point = new int[dim];
    std::vector<int> reset_levels = {1, 0};

    Autotuning *at = new Autotuning(min, max, ignore, new NelderMead(dim, e));

    point = at->execOffline(Test::function, dim);
    REQUIRE_MESSAGE(point[0] == Test::get_min_point()[0], "Base Case");

    for (auto &&level : reset_levels) {
      at->reset(level);
      if (level == 0) Test::reset();

      point = at->execOffline(Test::function, dim);
      REQUIRE_MESSAGE(point[0] == Test::get_min_point()[0],
                      "Reset level " + std::to_string(level));
    }
    delete[] point;
    delete at;
  }

  SECTION("MultiDimensional") {
    std::vector<double> erros = {10, 1, .1, .01, .001};
    std::vector<int> dimensions = {1, 2, 3, 4};

    for (auto &&e : erros) {
      for (auto &&dim : dimensions) {
        Test::reset();
        SECTION("Test with Dim = " + std::to_string(dim) +
                " and error = " + std::to_string(e)) {
          Autotuning *at =
              new Autotuning(min, max, ignore, new NelderMead(dim, e));
          int *point = new int[dim];

          point = at->execOffline(Test::function, dim);

          std::vector<int> min_point = Test::get_min_point();
          for (int j = 0; j < dim; j++) {
            REQUIRE(point[j] == min_point[j]);
          }

          delete at;
          delete[] point;
        }
      }
    }
  }
}