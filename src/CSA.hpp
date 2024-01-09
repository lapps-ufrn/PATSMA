/**
 * @file CSA.hpp
 * @brief Definition of the CSA class, a Coupled Simulated Annealing
 * optimization algorithm.
 */

#pragma once

#include <cmath>
#include <ctime>  // drand48_data

#include "NumericalOptimizer.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327 /**< PI mathematical constant. */
#endif

// Constants for default temperatures
#ifndef TG
#define TG 0.1 /**< Default generation temperature. */
#endif
#ifndef TA
#define TA 0.9 /**< Default acceptance temperature. */
#endif

#define END 0x99 /**< Flag indicating the end of the optimization process. */

/**
 * @class CSA
 * @brief Coupled Simulated Annealing (CSA) optimization algorithm.
 *
 * This class implements the CSA algorithm, a heuristic optimization method that
 * simulates the annealing process coupled with multiple optimizers.
 */
class CSA : public NumericalOptimizer {
  struct Opt {
    int id = 0;                /**< Identifier for the optimizer. */
    double *curSol = nullptr;  /**< Current solution vector. */
    double *probSol = nullptr; /**< Probable solution vector. */
    double curCost = 0.0;      /**< Current cost associated with the current solution. */
    double probCost = 0.0;     /**< Cost associated with the probable solution. */
    // Auxiliaries
    double prob = 0.0;          /**< Probability value. */
    double result = 0.0;        /**< Result of a random number generation. */
    struct drand48_data buffer; /**< Buffer for random number generation. */
  };

  int m_step;    /**< Current step in the optimization process. */
  int m_iter;    /**< Iteration count. */
  int m_maxIter; /**< Maximum number of iterations. */
  int m_iOpt;    /**< Iterator for Optimizers. */
  int m_nOpt;    /**< Number of Optimizers. */
  int m_dim;     /**< Number of Dimensions. */

  double m_tGen;     /**< Generation Temperature. */
  double m_tAcc;     /**< Acceptance Temperature. */
  double m_gamma;    /**< Gamma value used in the algorithm. */
  double m_maxCost;  /**< Maximum cost value. */
  double m_tmp;      /**< Temporary variable. */
  double m_probVar;  /**< Probability variable. */
  double *m_bestSol; /**< Best Solution vector. */
  double m_bestCost; /**< Best Cost relative to the best solution. */

  struct Opt *m_opts; /**< Array of Optimizers. */
  double *m_point;    /**< Point to return. */

  /**
   * @brief Copy the solution vector.
   * @param out The output solution vector.
   * @param in The input solution vector.
   */
  void copy_solution(double *out, double *in) const;

  /**
   * @brief Make a round shift for values < -1 and > 1.
   * @param value The m_point.
   * @return The m_point between -1 and 1.
   */
  static double rotate(double value);

  /**
   * @brief Switch values in vector position [i] from current solution to
   * solution, same from current cost to cost, and check if this new cost is the
   * maximum.
   * @param i Switch position.
   */
  void swap_opt_info(int i);

  /**
   * @brief Perform partial execution of the CSA algorithm.
   */
  void partial_exec();

 public:
  /**
   * @brief Constructor for the CSA class.
   * @param _num_opt Number of optimizers.
   * @param _dim Dimension of the cost function.
   * @param _max_iter Maximum number of iterations.
   */
  CSA(int _num_opt, int _dim, int _max_iter);

  /**
   * @brief Destructor for the CSA class.
   */
  ~CSA() override;

  /**
   * @brief Get the number of points (optimizers) used in the CSA algorithm.
   * @return Number of optimizers.
   */
  int getNumPoints() const override { return m_nOpt; }

  /**
   * @brief Get the dimension of the cost function.
   * @return Dimension of the cost function.
   */
  int getDimension() const override { return m_dim; }

  /**
   * @brief Check if the optimization process has reached its end.
   * @return True if the optimization is complete, false otherwise.
   */
  bool isEnd() const override { return m_step == END; }

  /**
   * @brief Run the CSA algorithm to generate the next solution.
   * @param cost The cost of the current solution.
   * @return A pointer to the next solution generated by the algorithm.
   */
  double *run(double cost) override;

  /**
   * @brief Print the parameters and state of the CSA algorithm.
   */
  void print() const override;

  /**
   * @brief Reset the CSA algorithm to a specified level.
   * @param level The level of resetting:
   *    - level 2: Reset the number of iterations.
   *    - level 1: Reset the points and temperatures (plus the previous ones).
   *    - level 0: Remove the best solution (plus the previous ones).
   */
  void reset(int level) override;

  /**
   * @brief Default constructor.
   *
   * This constructor is explicitly deleted, meaning objects of this class
   * cannot be created using the default constructor.
   */
  CSA() = delete;

  /**
   * @brief Copy assignment operator.
   *
   * This copy assignment operator is explicitly deleted, indicating that
   * instances of this class cannot be assigned using the copy assignment syntax.
   *
   * @param other The object to be copied.
   * @return Deleted - no copy assignment allowed.
   */
  CSA operator=(CSA other) = delete;

  /**
   * @brief Move assignment operator.
   *
   * This move assignment operator is explicitly deleted, indicating that
   * instances of this class cannot be assigned using the move assignment syntax.
   *
   * @param other The object to be moved.
   * @return Deleted - no move assignment allowed.
   */
  CSA &operator=(CSA &&other) = delete;

  /**
   * @brief Copy constructor.
   *
   * This copy constructor is explicitly deleted, meaning instances of this
   * class cannot be created using the copy constructor.
   *
   * @param other The object to be copied.
   */
  CSA(const CSA &other) = delete;

  /**
   * @brief Move constructor.
   *
   * This move constructor is explicitly deleted, meaning instances of this
   * class cannot be created using the move constructor.
   *
   * @param other The object to be moved.
   */
  CSA(CSA &&other) = delete;
};