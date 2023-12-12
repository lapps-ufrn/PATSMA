#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include "Autotuning.hpp"

#define Tolerance 0.000012
#define TRUE 1
#define FALSE 0

#define N 8001
#define THREAD_COUNT 6
#define PHI 1.61803398875

double **A, **B, diff;

void display(double **V, int n)
{
  for (int i = 0; i < n + 1; ++i)
  {
    for (int k = 0; k < n + 1; ++k)
    {
      printf("%f ", V[i][k]);
    }
    printf("\n");
  }
  printf("\n");
}

void initialize(double **A, int n)
{
  int i, j;

  for (j = 0; j < n + 1; j++)
  {
    A[0][j] = 1.0;
  }
  for (i = 1; i < n + 1; i++)
  {
    A[i][0] = 1.0;
    for (j = 1; j < n + 1; j++)
      A[i][j] = 0.0;
  }
}

long usecs(void)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec * 1000000 + t.tv_usec;
}

double matrix_calculation(double **A, int n, int *chunk)
{

  int iters = 0;
  int convergence = FALSE;
  double tmp;
  double diff;
  int i, j;

#pragma omp parallel num_threads(THREAD_COUNT) shared(chunk) private(tmp, i, j) reduction(+ : diff)
  {
#pragma omp for collapse(2) schedule(dynamic, *chunk)
    for (i = 1; i <= n; ++i)
    {
      for (j = 1; j <= n; ++j)
      {
        if ((i + j) % 2 == 1)
        {
          // printf("Thread: %d on (%d, %d)\n", omp_get_thread_num(), i, j);
          tmp = A[i][j];
          A[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i - 1][j] + A[i][j + 1] + A[i + 1][j]);
          diff += fabs(A[i][j] - tmp);
        }
      }
    }
#pragma omp barrier
  }

#pragma omp parallel num_threads(THREAD_COUNT) shared(chunk) private(tmp, i, j) reduction(+ : diff)
  {
#pragma omp for collapse(2) schedule(dynamic, *chunk)
    for (i = 1; i <= n; ++i)
    {
      for (j = 1; j <= n; ++j)
      {
        if ((i + j) % 2 == 0)
        {
          // printf("Thread: %d on (%d, %d)\n", omp_get_thread_num(), i, j);
          tmp = A[i][j];
          A[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i - 1][j] + A[i][j + 1] + A[i + 1][j]);
          diff += fabs(A[i][j] - tmp);
        }
      }
    }
#pragma omp barrier
  }
  return diff;
}

void solve_parallel(double **A, int n)
{
  /*PATSMA*/
  /* CSA parameters*/
  int dim = 1;
  int min = 1, max = N / (omp_get_num_threads() * 2);
  int ignore = 0;
  int n_opt = 4;
  int n_iter = 20;

  /*PATSMA INSTANCE*/
  Autotuning *at = new Autotuning(min,max,ignore,dim,n_opt,n_iter);

  /*Formula chunkExpert*/
  int number_iter = N;
  int number_threads = omp_get_num_threads();
  double f = std::floor(std::log2(static_cast<double>(number_iter) / number_threads)) * (1.0 / PHI);
  int chunk_size = static_cast<int>(std::floor(number_iter / (std::pow(2.0, f) * 2 * number_threads)));

  int *chunk = new int[dim];
  *chunk = chunk_size;
  // chunk = at->execOffline(chunk,matrix_calculation, A, N - 1);
  at->entireExecRuntime(matrix_calculation, chunk, A, N - 1);
  // at->singleExecRuntime(matrix_calculation,chunk, A, N - 1);

  printf("\n\n-----------------------Parallel Red Black Solver-----------------------\n\n\n");
  int for_iters;
  int iters = 0;
  int convergence = FALSE;
  double tmp;
  double diff;
  int i, j;
  for (for_iters = 1; for_iters < 21; ++for_iters)
  {
    diff = 0;

    diff = matrix_calculation(A, N - 1, chunk);
    iters++;

    printf("Difference after %3d iterations: %f\n", iters, diff);
    if (diff / ((double)N * (double)N) < Tolerance)
    {
      printf("\nConvergence achieved after %d iterations....Now exiting\n\n", iters);
      return;
    }
  }
  printf("\n\nIteration LIMIT Reached...Exiting\n\n");
}

int main(int argc, char *argv[])
{

  int i;
  long t_start, t_end;
  double time;
  A = new double *[N + 2];
  B = new double *[N + 2];
  for (i = 0; i < N + 2; i++)
  {
    A[i] = new double[N + 2];
    B[i] = new double[N + 2];
  }

  initialize(A, N);

  t_start = usecs();
  solve_parallel(A, N - 1);
  t_end = usecs();

  time = ((double)(t_end - t_start)) / 1000000;
  printf("Computation time for Parallel approach(secs): %f\n\n", time);
}