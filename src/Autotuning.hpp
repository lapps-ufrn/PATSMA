/**
 * @file Autotuning.hpp
 * @brief Header file for the Autotuning class
 */

#pragma once

#include "NumericalOptimizer.hpp"

/**
 * @brief Class for Autotuning
 */
class Autotuning {
  double m_min;  ///< Minimum value of the search interval
  double m_max;  ///< Maximum value of the search interval
  int m_iter;    ///< Iteration number
  int m_ignore;  ///< Number of iterations to ignore
  double m_cost;

  double m_t0, m_t1;  ///< Time starting and ending
  double m_runtime;   ///< Total time of a task

  NumericalOptimizer *p_optimizer;  ///< Numerical optimizer instance

  /**
   * @brief Rescale point from Type [-1,1] to Floating [min,max]
   * @tparam Point Type of optimization point
   * @param out Output int point
   * @param in Input double point
   */
  template <typename Point>
  void rescale(Point *out, const double *in) const;

  /**
   * @brief Check if the optimization has ended
   * @return True if optimization has ended, false otherwise
   */
  bool isEnd() const { return p_optimizer->isEnd(); }

 public:
  /**
   * @brief Run the autotuning algorithm
   * @tparam Point Type of optimization point
   * @param point Input/output array of tuning parameters
   * @param _cost Cost value for the current iteration
   */
  template <typename Point>
  void exec(Point *point, double cost);
  /**
   * @brief Start a new iteration of the autotuning algorithm
   * @tparam Point Type of optimization point
   * @param point Input/output array of tuning parameters
   */
  template <typename Point>
  void start(Point *point);

  /**
   * @brief End the current iteration of the autotuning algorithm
   */
  void end();

  /**
   * @brief Print basic information about the autotuning parameters
   */
  void print() const;

  /**
   * @brief Reset the autotuning and numerical optimizer
   * @param level Reset level:
   *   - level 2: Reset the number of iterations
   *   - level 1: Reset the points and the temperatures (plus the previous ones)
   *   - level 0: Remove the best solution (plus the previous ones)
   */
  void reset(int level);

  /**
   * @brief Execute the entire autotuning, recommended to be done off the
   * application loop. It uses "funciton" runtime as cost.
   * @tparam Point Type of optimization point
   * @tparam Func Type of the cost function
   * @tparam Args Types of additional arguments to the cost function
   * @param function Function to be optimized
   * @param point Input/output array of tuning parameters
   * @param args Additional arguments to the function, the first argument needs
   * to be the optimization point
   * @return Optimized point
   */
  template <typename Point = int, typename Func, typename... Args>
  void entireExecRuntime(Func function, Point *point, Args... args);

  /**
   * @brief Execute the entire autotuning, recommended to be done off the
   * application loop. It uses funciton return as cost.
   * @tparam Point Type of optimization point
   * @tparam Func Type of the cost function
   * @tparam Args Types of additional arguments to the cost function
   * @param function Cost function to be optimized. It needs to return the cost
   * and the first argument needs to be the optimization point
   * @param point Input/output array of tuning parameters
   * @param args Additional arguments to the cost function
   * @return Optimized point
   */
  template <typename Point = int, typename Func, typename... Args>
  void entireExec(Func function, Point *point, Args... args);

  /**
   * @brief Execute one iteration of the autotuning, appropriate to be done on
   * the application loop. It uses "funciton" runtime as cost.
   * @tparam Func Type of the function
   * @tparam Args Types of additional arguments to the function
   * @param function Function to be optimized. The first argument needs
   * to be the optimization point
   * @param point Input/output array of tuning parameters
   * @param args Additional arguments to the function
   */
  template <typename Point = int, typename Func, typename... Args>
  void singleExecRuntime(Func function, Point *point, Args... args);

  /**
   * @brief Execute one iteration of the autotuning, appropriate to be done on
   * the application loop. It uses funciton return as cost.
   * @tparam Point Type of optimization point
   * @tparam Func Type of the cost function
   * @tparam Args Types of additional arguments to the cost function
   * @param function Cost function to be optimized. It needs to return the cost
   * and the first argument needs to be the optimization point
   * @param point Input/output array of tuning parameters
   * @param args Additional arguments to the cost function
   */
  template <typename Point = int, typename Func, typename... Args>
  void singleExec(Func function, Point *point, Args... args);

  /**
   * @brief Deleted copy assignment operator
   * @param Autotuning Object to be assigned
   * @return Deleted assignment operator
   */
  auto operator=(Autotuning) -> Autotuning = delete;

  /**
   * @brief Deleted move assignment operator
   * @param Autotuning Object to be moved
   * @return Deleted move assignment operator
   */
  auto operator=(Autotuning &&) -> Autotuning & = delete;

  /**
   * @brief Deleted copy constructor
   * @param Autotuning Object to be copied
   */
  Autotuning(const Autotuning &) = delete;

  /**
   * @brief Deleted move constructor
   * @param Autotuning Object to be moved
   */
  Autotuning(Autotuning &&) = delete;

  /**
   * @brief Default constructor
   */
  Autotuning() = default;

  /**
   * @brief Parameterized constructor with CSA as the default optimizer
   * @param dim Cost Function Dimension
   * @param _min Minimum value of the search interval
   * @param _max Maximum value of the search interval
   * @param _ignore Number of iterations to ignore
   * @param num_opt Number of optimizers
   * @param max_iter Maximum number of iterations
   */
  Autotuning(int dim, double min, double max, int ignore, int num_opt, int max_iter);

  /**
   * @brief Parameterized constructor with a custom optimizer
   * @param _min Minimum value of the search interval
   * @param _max Maximum value of the search interval
   * @param _ignore Number of iterations to ignore
   * @param _optimizer Numerical optimizer instance
   */
  Autotuning(double min, double max, int ignore, NumericalOptimizer *optimizer);

  /**
   * @brief Destructor
   */
  ~Autotuning();
};

#include "Autotuning.tpp"
