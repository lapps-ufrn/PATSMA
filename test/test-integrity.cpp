#include <math.h>
#include <omp.h>
#include <unistd.h>

#include <catch2/catch.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "Autotuning.hpp"
#include "Test.hpp"

#define CHECK_MESSAGE(cond, msg) \
  do {                           \
    INFO(msg);                   \
    CHECK(cond);                 \
  } while ((void)0, 0)

TEST_CASE("CSA") {
  int dim = 1;
  int min = -100, max = 100;
  int ignore = 0;
  int n_opt = 4;
  int n_iter = 20;

  SECTION("Test exception throwing") {
    CHECK_THROWS_AS(new Autotuning(min, max, -1, nullptr), std::invalid_argument);
    CHECK_THROWS_WITH(new Autotuning(min, max, -1, nullptr),
                      "Ignore Value Invalid! Set _ignore >= 0.");

    CHECK_THROWS_AS(new CSA(0, n_opt, n_iter), std::invalid_argument);
    CHECK_THROWS_WITH(new CSA(0, n_opt, n_iter), "Dimensional Value Invalid! Set dim > 0.");
    CHECK_THROWS_AS(new CSA(dim, 0, n_iter), std::invalid_argument);
    CHECK_THROWS_WITH(new CSA(dim, 0, n_iter), "Optmizers Number Invalid! Set num_opt > 0.");
    CHECK_THROWS_AS(new CSA(dim, n_opt, 0), std::invalid_argument);
    CHECK_THROWS_WITH(new CSA(dim, n_opt, 0),
                      "Max number of intereration Invalid! Set max_iter > 0.");
  }

  SECTION("Running methods") {
    int *point = new int[dim];
    auto *at = new Autotuning(min, max, ignore, dim, n_opt, n_iter);

    SECTION("singleExec") {
      at->singleExec(Test::function, point, dim);
      CHECK(point[0] == Test::get_min_point()[0]);
    }

    SECTION("entireExec") {
      at->entireExec(Test::function, point, dim);
      CHECK(point[0] == Test::get_min_point()[0]);
    }

    SECTION("singleExecRuntime") {
      at->singleExecRuntime(Test::sleep_function, point);
      CHECK(point[0] == Test::get_min_point()[0]);
    }

    SECTION("entireExecRuntime") {
      at->entireExecRuntime(Test::sleep_function, point);
      CHECK(point[0] == Test::get_min_point()[0]);
    }
  }

  SECTION("Reset Cases") {
    Test::reset();
    std::vector<int> reset_levels = {0, 1, 2};
    int *point = new int[dim];
    auto *at = new Autotuning(min, max, ignore, dim, n_opt, n_iter);

    for (auto &&level : reset_levels) {
      SECTION("singleExec") {
        at->singleExec(Test::function, point, dim);
        CHECK(point[0] == Test::get_min_point()[0]);
        at->reset(level);
        if (level == 2) Test::reset();
        at->singleExec(Test::function, point, dim);
        CHECK_MESSAGE(point[0] == Test::get_min_point()[0], "Reset level " + std::to_string(level));
      }

      SECTION("entireExec") {
        at->entireExec(Test::function, point, dim);
        CHECK(point[0] == Test::get_min_point()[0]);
        at->reset(level);
        if (level == 2) Test::reset();
        at->entireExec(Test::function, point, dim);
        CHECK_MESSAGE(point[0] == Test::get_min_point()[0], "Reset level " + std::to_string(level));
      }
    }
    delete at;
    delete point;
  }

  SECTION("MultiDimensional") {
    Test::reset();
    std::vector<int> num_iters = {20, 200, 2000};
    std::vector<int> dimensions = {1, 2, 3, 4};

    for (auto &&_n_iter : num_iters) {
      for (auto &&_dim : dimensions) {
        SECTION("Test with Dim = " + std::to_string(_dim) +
                " and n_iter = " + std::to_string(_n_iter)) {
          Test::reset();
          Autotuning *at = new Autotuning(min, max, ignore, _dim, n_opt, _n_iter);
          int *point = new int[_dim];

          at->singleExec(Test::function, point, _dim);

          std::vector<int> min_point = Test::get_min_point();
          for (int j = 0; j < _dim; j++) {
            CHECK(point[j] == min_point[j]);
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
    CHECK_THROWS_AS(new Autotuning(min, max, -1, nullptr), std::invalid_argument);
    CHECK_THROWS_WITH(new Autotuning(min, max, -1, nullptr),
                      "Ignore Value Invalid! Set _ignore >= 0.");

    CHECK_THROWS_AS(new NelderMead(0, e), std::invalid_argument);
    CHECK_THROWS_WITH(new NelderMead(0, e), "Dimensional Value Invalid! Set dim > 0.");
    CHECK_THROWS_AS(new NelderMead(dim, -.1), std::invalid_argument);
    CHECK_THROWS_WITH(new NelderMead(dim, -.1), "Invalid m_error! Set error >= 0.");
  }

  SECTION("Reset Cases") {
    Test::reset();

    int *point = new int[dim];
    std::vector<int> reset_levels = {0, 1};

    Autotuning *at = new Autotuning(min, max, ignore, new NelderMead(dim, e));

    at->singleExec(Test::function, point, dim);
    CHECK_MESSAGE(point[0] == Test::get_min_point()[0], "Base Case");

    for (auto &&level : reset_levels) {
      at->reset(level);
      if (level == 1) Test::reset();

      at->singleExec(Test::function, point, dim);
      CHECK_MESSAGE(point[0] == Test::get_min_point()[0], "Reset level " + std::to_string(level));
    }
    delete[] point;
    delete at;
  }

  SECTION("MultiDimensional") {
    std::vector<double> erros = {10, 1, .1, .01, .001};
    std::vector<int> dimensions = {1, 2, 3, 4};

    for (auto &&_e : erros) {
      for (auto &&_dim : dimensions) {
        Test::reset();
        SECTION("Test with Dim = " + std::to_string(_dim) + " and error = " + std::to_string(_e)) {
          Autotuning *at = new Autotuning(min, max, ignore, new NelderMead(_dim, _e));
          int *point = new int[_dim];

          at->singleExec(Test::function, point, _dim);

          std::vector<int> min_point = Test::get_min_point();
          for (int j = 0; j < _dim; j++) {
            CHECK(point[j] == min_point[j]);
          }

          delete at;
          delete[] point;
        }
      }
    }
  }
}