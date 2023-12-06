/**
 * @file NelderMead.hpp
 * @brief Declaration of the NelderMead class for numerical optimization using
 * the Nelder-Mead algorithm.
 */
#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include "NumericalOptimizer.hpp"

// Constants for algorithm parameters
#ifndef A_LFA
#define A_LFA 1
#endif
#ifndef G_AMA
#define G_AMA 2
#endif
#ifndef R_HO
#define R_HO 0.5
#endif
#ifndef S_GMA
#define S_GMA 0.5
#endif

#define END 0x99

/**
 * @class NelderMead
 * @brief Implementation of the Nelder-Mead algorithm for numerical
 * optimization.
 */
class NumericalOptimizer;
class NelderMead : public NumericalOptimizer {
  // Algorithm parameters
  double m_alpha;  ///< Used in reflection
  double m_gamma;  ///< Used in expansion
  double m_rho;    ///< Used in contraction
  double m_sigma;  ///< Used in reduction

  // Variables for algorithm state
  int m_step;      ///< Current step in the algorithm
  int m_worstID;   ///< Index of the worst solution
  int m_secondID;  ///< Index of the second worst solution
  int m_bestID;    ///< Index of the best solution

  // Problem-specific dimensions
  int m_nPoints;  ///< Number of solutions in the algorithm
  int m_dim;      ///< Dimensionality of the problem
  int m_iPoint;   ///< Index of the current point

  // Cost variables for different points
  double m_costReflection;   ///< Cost of reflected point
  double m_costExpansion;    ///< Cost of expanded point
  double m_costContraction;  ///< Cost of contracted point

  // Points generated by different steps
  double *m_pointReflection;   ///< Point generated by reflection
  double *m_pointExpansion;    ///< Point generated by expansion
  double *m_pointContraction;  ///< Point generated by contraction

  // Arrays to store points, costs, and centroid
  double **m_points;    ///< Vector where each line will be a solution
                        ///< 'dim'-dimensional
  double *m_bestPoint;  ///< Best solution found
  double *m_costs;      ///< Costs associated with each solution
  double *m_centroid;   ///< Centroid of the solutions

  double m_result;             ///< Temporary double variables
  struct drand48_data buffer;  ///< To use in random seed

  double m_error;  ///< Error threshold for terminating the algorithm

  // Private member functions
  void set_points();                           ///< Set initial points for the algorithm
  void sort_points();                          ///< Sort solutions based on their costs
  void calculate_centroid();                   ///< Calculate the centroid of the solutions
  double volume();                             ///< Calculate the volume of the solutions
  static void swap(double *&p1, double *&p2);  ///< Swap two pointers
  static void swap(double &p1, double &p2);    ///< Swap two values

 public:
  /**
   * @brief Get the number of points used in the algorithm.
   * @return Number of points.
   */
  int getNumPoints() const override { return m_nPoints; };

  /**
   * @brief Get the dimensionality of the problem.
   * @return Dimensionality.
   */
  int getDimension() const override { return m_dim; };

  /**
   * @brief Reset the state of the algorithm.
   * @param level Reset level, indicating how to reset the algorithm.
   */
  void reset(int level) override;

  /**
   * @brief Check if the optimization has reached the end.
   * @return True if the optimization has ended, false otherwise.
   */
  bool isEnd() const override { return m_step == END; }

  /**
   * @brief Run one iteration of the Nelder-Mead algorithm.
   * @param _cost The cost associated with the current solution.
   * @return Pointer to the next solution.
   */
  double *run(double _cost) override;

  /**
   * @brief Constructor for the NelderMead class.
   * @param dim Dimensionality of the problem.
   * @param error Error threshold for terminating the algorithm.
   */
  NelderMead(int dim, double error);

  /**
   * @brief Destructor for the NelderMead class.
   */
  ~NelderMead() override;

  /**
   * @brief Default constructor (deleted).
   */
  NelderMead() = delete;

  /**
   * @brief Copy constructor (deleted).
   */
  NelderMead(const NelderMead &) = delete;

  /**
   * @brief Move constructor (deleted).
   */
  NelderMead(NelderMead &&) = delete;

  /**
   * @brief Copy assignment operator (deleted).
   */
  NelderMead operator=(NelderMead) = delete;

  /**
   * @brief Move assignment operator (deleted).
   */
  NelderMead &operator=(NelderMead &&) = delete;
};
