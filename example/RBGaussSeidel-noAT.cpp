#include <omp.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#define N 2001
#define MAX_ITER 1000

void initialize(double **A, int n) {
  int i, j;
  for (j = 0; j < n + 1; j++) {
    A[0][j] = 1.0;
  }
  for (i = 1; i < n + 1; i++) {
    A[i][0] = 1.0;
    for (j = 1; j < n + 1; j++) {
      A[i][j] = 0.0;
    }
  }
}

double matrix_calculation(double **A, int n) {
  double tmp;
  double diff = 0;
  int i, j;

#ifdef _OPENMP
#pragma omp parallel private(tmp, i, j)
  {
#pragma omp for reduction(+ : diff)
#endif
    for (i = 1; i <= n; ++i) {
      for (j = 1; j <= n; ++j) {
        if ((i + j) % 2 == 1) {
          tmp = A[i][j];
          A[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i - 1][j] + A[i][j + 1] + A[i + 1][j]);
          diff += fabs(A[i][j] - tmp);
        }
      }
    }
#ifdef _OPENMP
#pragma omp for reduction(+ : diff)
#endif
    for (i = 1; i <= n; ++i) {
      for (j = 1; j <= n; ++j) {
        if ((i + j) % 2 == 0) {
          tmp = A[i][j];
          A[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i - 1][j] + A[i][j + 1] + A[i + 1][j]);
          diff += fabs(A[i][j] - tmp);
        }
      }
    }
#ifdef _OPENMP
  }
#endif
  return diff;
}

void solve_parallel(double **A, int n) {
  double diff = 0;

  printf("\n\n-----------------------Parallel Red Black Solver-----------------------\n\n\n");

  for (int iters = 1; iters < MAX_ITER; ++iters) {
    diff = matrix_calculation(A, N - 1);
  }

  printf("Difference after %3d iterations: %f\n", iters, diff);
  printf("\n\nIteration LIMIT Reached...Exiting\n\n");
}

int main(int argc, char *argv[]) {
  int i;
  double t_start, t_end;
  double **A;
  A = new double *[N + 2];
  for (i = 0; i < N + 2; i++) {
    A[i] = new double[N + 2];
  }

  initialize(A, N);

  t_start = omp_get_wtime();
  solve_parallel(A, N - 1);
  t_end = omp_get_wtime();

  printf("Computation time for Parallel approach(secs) on %i threads: %f\n\n",
         omp_get_max_threads(), t_end - t_start);
}