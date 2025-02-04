#include "NelderMead.hpp"

#include <cstring>   // memcpy
#include <iostream>  // cout, endl, bad_alloc

#ifndef DBL_MAX
#include <limits>  // numeric_limits
#define DBL_MAX std::numeric_limits<double>::max()
#endif

NelderMead::NelderMead(int dim, double error) : m_dim(dim), m_error(error) {
  if (dim < 1) {
    throw std::invalid_argument("Dimensional Value Invalid! Set dim > 0.");
  }
  if (error < 0) {
    throw std::invalid_argument("Invalid m_error! Set error >= 0.");
  }
  m_nPoints = (dim > 2) ? (dim + 1) : 3;

  try {
    p_points = new double *[m_nPoints];
    p_costs = new double[m_nPoints];
    p_centroid = new double[m_dim];
    p_bestPoint = new double[m_dim];
    p_pointReflection = new double[m_dim];   // point generated by reflection
    p_pointExpansion = new double[m_dim];    // point generated by expansion
    p_pointContraction = new double[m_dim];  // point generated by contraction
    for (int i = 0; i < m_nPoints; i++) {
      p_points[i] = new double[m_dim];
    }
  } catch (const std::bad_alloc &e) {
    std::cout << "Memory Allocation"
              << " is failed: " << e.what() << std::endl;
  }

  m_worstID = m_nPoints - 1;   // index of the worst solution
  m_secondID = m_nPoints - 2;  // index of the second worst solution
  m_bestID = 0;                // index of the best solution

  m_costReflection = 0.0;   // Cost of reflected point
  m_costExpansion = 0.0;    // Cost of expanded point
  m_costContraction = 0.0;  // Cost of contracted point

  m_alpha = A_LFA;  // Used in reflection
  m_gamma = G_AMA;  // Used in expansion
  m_rho = R_HO;     // Used in contraction
  m_sigma = S_GMA;  // Used in reduction

  m_result = 0.0;
  m_iPoint = 0;
  m_step = 0;

  for (int i = 0; i < m_nPoints; i++) {
    p_costs[i] = DBL_MAX;
  }

  // Generate inital p_points
  srand48_r(time(NULL), &buffer);  // random seed
  for (int i = 0; i < m_nPoints; i++) {
    for (int j = 0; j < m_dim; j++) {
      drand48_r(&buffer, &m_result);            // random number in 'm_result'
      p_points[i][j] = (m_result * 2.0 - 1.0);  // numbers between -1 and 1.
    }
  }
}

double *NelderMead::run(double _cost) {
  do {
    switch (m_step) {
      case 0:  // Initialize Points
        if (m_iPoint > 0) {
          p_costs[m_iPoint - 1] = _cost;
        }
        if (m_iPoint < m_nPoints) {
          m_iPoint++;
          return p_points[m_iPoint - 1];
        }

      case 1:  // Initialize Variables
        sort_points();
        calculate_centroid();

      case 2:  // Reflection - Use m_alpha
        for (int j = 0; j < m_dim; j++) {
          p_pointReflection[j] =
              fmod(p_centroid[j] + m_alpha * (p_centroid[j] - p_points[m_worstID][j]), 1);
        }
        m_step = 3;
        return p_pointReflection;

      case 3:
        m_costReflection = _cost;
        if (p_costs[m_bestID] <= m_costReflection && m_costReflection < p_costs[m_secondID]) {
          swap(p_pointReflection, p_points[m_worstID]);
          p_costs[m_worstID] = m_costReflection;
          break;
        } else if (m_costReflection < p_costs[m_bestID]) {
          // Calculate Expansion - Use m_gamma
          for (int j = 0; j < m_dim; j++) {
            p_pointExpansion[j] =
                fmod(p_centroid[j] + m_gamma * (p_pointReflection[j] - p_centroid[j]), 1);
          }
          m_step = 4;
          return p_pointExpansion;
        } else if (p_costs[m_secondID] <= m_costReflection &&
                   m_costReflection < p_costs[m_worstID]) {
          // Calculate Outside Contraction - Use m_rho
          for (int j = 0; j < m_dim; j++) {
            p_pointContraction[j] =
                fmod(p_centroid[j] + m_rho * (p_pointReflection[j] - p_centroid[j]), 1);
          }
          m_step = 5;
          return p_pointContraction;
        } else if (m_costReflection >= p_costs[m_worstID]) {
          // Calculate Inside Contraction - Use m_rho
          for (int j = 0; j < m_dim; j++) {
            p_pointContraction[j] =
                fmod(p_centroid[j] - m_rho * (p_centroid[j]) - p_points[m_worstID][j], 1);
          }
          m_step = 6;
          return p_pointContraction;
        }

      case 4:
        m_costExpansion = _cost;
        if (m_costExpansion < m_costReflection) {
          swap(p_pointExpansion, p_points[m_worstID]);
          p_costs[m_worstID] = m_costExpansion;
        } else {
          swap(p_pointReflection, p_points[m_worstID]);
          p_costs[m_worstID] = m_costReflection;
        }
        break;

      case 5:
        m_costContraction = _cost;
        if (m_costContraction <= m_costReflection) {
          swap(p_pointContraction, p_points[m_worstID]);
          p_costs[m_worstID] = m_costContraction;
        } else {
          m_step = 7;
        }
        break;

      case 6:
        m_costContraction = _cost;
        if (m_costContraction < p_costs[m_worstID]) {
          swap(p_pointContraction, p_points[m_worstID]);
          p_costs[m_worstID] = m_costContraction;
          break;
        } else {
          m_step = 7;
          break;
        }
    }

    sort_points();
    if (m_step == 7) {
      // Replace all p_points, except the p_bestPoint
      for (int i = 0; i < m_nPoints; i++) {
        if (i != m_bestID) {
          for (int j = 0; j < m_dim; j++) {
            p_points[i][j] =
                fmod(p_points[m_bestID][j] + m_sigma * (p_points[i][j] - p_points[m_bestID][j]), 1);
          }
        }
      }
      m_step = 0;
      m_iPoint = 1;
    }

    calculate_centroid();
    m_step = 2;

  } while (volume() > m_error);

  m_step = END;
  memcpy(p_bestPoint, p_points[m_bestID], m_dim * sizeof(double));

  return p_bestPoint;
}

void NelderMead::sort_points() {
  // Put in order solutions, ie, p_points[0][x] is the best, ...,
  // p_points[n+1][x] is the worst Bubble sort
  for (int i = 0; i < m_nPoints - 1; i++) {
    for (int j = i + 1; j < m_nPoints; j++) {
      if (p_costs[i] > p_costs[j]) {
        swap(p_points[i], p_points[j]);
        swap(p_costs[i], p_costs[j]);
      }
    }
  }
}

void NelderMead::calculate_centroid() {
  // Calculate p_centroid
  for (int j = 0; j < m_dim; j++) {
    m_result = 0.0;
    p_centroid[j] = 0.0;  // Centroid = {0}
    for (int i = 0; i < m_nPoints; i++) {
      // The worst solution is not necessary
      if (i != m_worstID) {
        p_centroid[j] += p_points[i][j];
        m_result += 1.0;
      }
    }
    p_centroid[j] /= (m_result);
  }
}

void NelderMead::reset(int level) {
  m_costReflection = 0.0;   // Cost of reflected point
  m_costExpansion = 0.0;    // Cost of expanded point
  m_costContraction = 0.0;  // Cost of contracted point

  m_result = 0.0;
  m_iPoint = 0;
  m_step = 0;

  sort_points();
  switch (level) {
    case 1:  // Reset with random p_points removing the best point
      for (int j = 0; j < m_dim; j++) {
        drand48_r(&buffer, &m_result);                   // random number in 'm_result'
        p_points[m_bestID][j] = (m_result * 2.0 - 1.0);  // numbers between -1 and 1.
      }
    case 0:  // Reset with random p_points keeping the best point
      for (int i = 0; i < m_nPoints; i++) {
        if (i != m_bestID) {
          for (int j = 0; j < m_dim; j++) {
            drand48_r(&buffer, &m_result);            // random number in 'm_result'
            p_points[i][j] = (m_result * 2.0 - 1.0);  // numbers between -1 and 1.
          }
        }
      }
      break;

    default:
      throw std::runtime_error("There is not the Nelder-Mead reset option level " +
                               std::to_string(level));
      break;
  }
}

double NelderMead::volume() {
  double total = 0.0;
  for (int i = 0; i < m_nPoints; i++) {
    double value = 0.0;
    for (int j = 0; j < m_dim; j++) {
      value += pow(p_points[i][j] - p_centroid[j], 2.0);
    }
    value = sqrt(value);      // With this calculation, 'value' is the norm.
    value = pow(value, 2.0);  // value is equal to norm²
    total += value;
  }
  return sqrt(total / m_nPoints);
}

void NelderMead::swap(double *&p1, double *&p2) {
  double *temp = p1;
  p1 = p2;
  p2 = temp;
}

void NelderMead::swap(double &p1, double &p2) {
  double temp = p1;
  p1 = p2;
  p2 = temp;
}

NelderMead::~NelderMead() {
  for (int i = 0; i < m_nPoints; i++) {
    delete[] p_points[i];
  }
  delete[] p_points;
  delete[] p_costs;
  delete[] p_centroid;
  delete[] p_bestPoint;
  delete[] p_pointReflection;
  delete[] p_pointExpansion;
  delete[] p_pointContraction;
}