#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>

using namespace std;

// Simple helper class to measure time.
struct Timer {
  using clock = std::chrono::steady_clock;

  clock::time_point startTime;

  Timer() : startTime(clock::now()) {}

  double getSeconds() {
    return std::chrono::duration_cast<std::chrono::duration<float>>(clock::now() - startTime).count();
  }
};

// Verify result
bool verify(double* result, double* ref, long size) {
  for (int i = 0; i < size; i++) {
    if (fabs(result[i] - ref[i]) >= 1e-3)
        return false;
  }
  return true;
}


// Write reference file
bool writeToFile(std::string file, double* vals, long size) {
  std::fstream fout(file, std::ios_base::out);
  if (!fout.is_open()) {
    std::cout << "Unable to open file" << std::endl;
    return false;
  }
  for (int i = 0; i < size; i++) {
    fout << vals[i] << " ";
  }
  return true;
}

// Read reference file
bool readFromFile(std::string file, double* vals, long size) {
  std::fstream fin(file, std::ios_base::in);
  long idx = 0;
  double val;
  while (fin >> val) {
    if (idx >= size) {
        std::cout << "Too many values in reference file" << std::endl;
        return false;
    }
    vals[idx++] = val;
  }
  if (size != idx) {
    std::cout << "Too few values in reference file" << std::endl;
    return false;
  }
  return true;
}

int main(int argc, char *argv[])

//
//  Purpose:
//
//    MAIN is the main program for HEATED_PLATE_MPI.
//
//  Discussion:
//
//    This code solves the steady state heat equation on a rectangular region.
//
//    The sequential version of this program needs approximately
//    18/epsilon iterations to complete.
//
//
//    The physical region, and the boundary conditions, are suggested
//    by this diagram;
//
//                   W = 0
//             +------------------+
//             |                  |
//    W = 100  |                  | W = 100
//             |                  |
//             +------------------+
//                   W = 100
//
//    The region is covered with a grid of M by N nodes, and an M by N
//    array W is used to record the temperature.  The correspondence between
//    array indices and locations in the region is suggested by giving the
//    indices of the four corners:
//
//                  I = 0
//          [0][0]-------------[0][N-1]
//             |                  |
//      J = 0  |                  |  J = N-1
//             |                  |
//        [M-1][0]-----------[M-1][N-1]
//                  I = M-1
//
//    The steady state solution to the discrete heat equation satisfies the
//    following condition at an interior grid point:
//
//      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )
//
//    where "Central" is the index of the grid point, "North" is the index
//    of its immediate neighbor to the "north", and so on.
//
//    Given an approximate solution of the steady state heat equation, a
//    "better" solution is given by replacing each interior point by the
//    average of its 4 neighbors - in other words, by using the condition
//    as an ASSIGNMENT statement:
//
//      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )
//
//    If this process is repeated often enough, the difference between
//    successive estimates of the solution will go to zero.
//
//    This program carries out such an iteration, using a tolerance specified by
//    the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2011
//
//  Author:
//
//    Original C version by Michael Quinn.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Michael Quinn,
//    Parallel Programming in C with MPI and OpenMP,
//    McGraw-Hill, 2004,
//    ISBN13: 978-0071232654,
//    LC: QA76.73.C15.Q55.
//
//  Local parameters:
//
//    Local, double DIFF, the norm of the change in the solution from one
//    iteration to the next.
//
//    Local, double MEAN, the average of the boundary values, used to initialize
//    the values of the solution in the interior.
//
//    Local, double U[M][N], the solution at the previous iteration.
//
//    Local, double W[M][N], the solution computed at the latest iteration.
//
{
#define M 500
#define N 500

  double diff;
  double epsilon = 0.001;
  int i;
  int iterations;
  int iterations_print;
  int j;
  double mean;

  // TODO: Each process may only allocate a part of the matrix
  // you may want to allocate the matrix wit malloc depening on the number of
  // lines each process may get
  double u[M][N];
  double w[M][N];
  double wtime;


  // TODO: Make sure that only one MPI process prints the output
  cout << "\n";
  cout << "HEATED_PLATE\n";
  cout << "  C++ version\n";
  cout
      << "  A program to solve for the steady state temperature distribution\n";
  cout << "  over a rectangular plate.\n";
  cout << "\n";
  cout << "  Spatial grid of " << M << " by " << N << " points.\n";
  cout << "  The iteration will be repeated until the change is <= " << epsilon
       << "\n";

  // TODO: Get number of MPI processes
  long num_procs = 1;
  cout << "  Number of processes available = " << num_procs << "\n";

  std::string refFile;

  if (argc >= 2) {
    refFile = argv[1];
    std::cout << "  Reading reference results from " << refFile << std::endl;
  } else {
    std::cout << "  No reference file provided - skipping verification." << std::endl;
  }

  //
  //  Set the boundary values, which don't change.
  // TODO: Be aware that each process hold a different part of the matrix and may only
  // initialize his own part
  mean = 0.0;

  for (i = 1; i < M - 1; i++) {
    w[i][0] = 100.0;
  }
  for (i = 1; i < M - 1; i++) {
    w[i][N - 1] = 100.0;
  }
  for (j = 0; j < N; j++) {
    w[M - 1][j] = 100.0;
  }
  for (j = 0; j < N; j++) {
    w[0][j] = 0.0;
  }
  //
  //  Average the boundary values, to come up with a reasonable
  //  initial value for the interior.
  //  TODO: Keep in mind that you need to average over all processes
  for (i = 1; i < M - 1; i++) {
    mean = mean + w[i][0] + w[i][N - 1];
  }
  for (j = 0; j < N; j++) {
    mean = mean + w[M - 1][j] + w[0][j];
  }

  mean = mean / (double)(2 * M + 2 * N - 4);
  cout << "\n";
  cout << "  MEAN = " << mean << "\n";
  //
  //  Initialize the interior solution to the mean value.
  //
  for (i = 1; i < M - 1; i++) {
    for (j = 1; j < N - 1; j++) {
      w[i][j] = mean;
    }
  }
  //
  //  Iterate until the  new solution W differs from the old solution U
  //  by no more than EPSILON.
  //
  iterations = 0;
  iterations_print = 1;
  cout << "\n";
  cout << " Iteration  Change\n";
  cout << "\n";

  Timer timer;

  diff = epsilon;

cout << "------------------" << std::endl;

  // TODO: This statement should yield the same result on all processes, so that all
  // processes stop at the same iteration
  while (epsilon <= diff) {
    //
    //  Save the old solution in U.
    //
    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        u[i][j] = w[i][j];
      }
    }
    //
    //  Determine the new estimate of the solution at the interior points.
    //  The new solution W is the average of north, south, east and west
    //  neighbors.
    //  TODO: Here you may need parts of the matrix that are part of other processes
    for (i = 1; i < M - 1; i++) {
      for (j = 1; j < N - 1; j++) {
        w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;
      }
    }

    diff = 0.0;

    // TODO: Be aware that each process only computes its local diff here. You may
    // need to combine the results from all processes
    for (i = 1; i < M - 1; i++) {
      for (j = 1; j < N - 1; j++) {
        if (diff < fabs(w[i][j] - u[i][j])) {
          diff = fabs(w[i][j] - u[i][j]);
        }
      }
    }

    iterations++;
    if (iterations == iterations_print) {
      cout << "  " << setw(8) << iterations << "  " << diff << "\n";
      iterations_print = 2 * iterations_print;
    }
  }
  
  // TODO: Insert a Barrier before time measurement
  wtime = timer.getSeconds();

  // TODO: Make sure that only one MPI process prints the output
  cout << "\n";
  cout << "  " << setw(8) << iterations << "  " << diff << "\n";
  cout << "\n";
  cout << "  Error tolerance achieved.\n";
  cout << "  Wallclock time = " << setprecision(3) << wtime << "s\n";
  //
  //  Terminate.
  //
  cout << "\n";
  cout << "HEATED_PLATE:\n";
  cout << "  Normal end of execution.\n";


  if (!refFile.empty()) {
    double ref[M][N];
    if (!readFromFile(refFile, &ref[0][0], M * N)) {
      std::cout << "  Verification failed - could not load reference file." << std::endl;
      return 0;
    } 
    if (verify(&w[0][0], &ref[0][0], M * N)) {
      std::cout << "  Result is valid." << std::endl;
    } else {
      std::cout << "  Verification failed!" << std::endl;
    }
  }


  return 0;

#undef M
#undef N
}
