#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>
#include <mpi.h>

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
bool verify(double* result, double* ref, long size, int rank) {
    for (int i = 0; i < size; i++) {
        if (fabs(result[i] - ref[i]) >= 1e-3)
        {
            std::cout<<"DEBUG R" << rank << ":\tRes: " << result[i] << ", Ref: " << ref[i] << ", Fabs: " << fabs(result[i] - ref[i]) << std::endl;
            return false;
        }
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

    MPI_Init(&argc, &argv);
    int mpi_init = 0;
    MPI_Initialized(&mpi_init);
    double diff;
    double epsilon = 0.001;
    int i;
    int iterations;
    int iterations_print;
    int j;
    double mean;
    int rank;
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Request request0;
    MPI_Request requestN;
    MPI_Request requestK;
    MPI_Request requestK2;

    // TODO: Each process may only allocate a part of the matrix
    // you may want to allocate the matrix wit malloc depening on the number of
    // lines each process may get
    if(M < num_procs)
    {
        MPI_Finalize();
        return 1;
    }
    int rows = M / num_procs;
    if(rank < M % num_procs) //Distribute the remaining rows in numerical order
        rows += 1;

    double u[rows][N];
    double w[rows][N];
    double wtime;
    std::string refFile;

    // TODO: Make sure that only one MPI process prints the output
    if(rank == 0) {
        cout << "\n";
        cout << "HEATED_PLATE Asynchronous\n";
        cout << "  C++ version\n";
        cout
                << "  A program to solve for the steady state temperature distribution\n";
        cout << "  over a rectangular plate.\n";
        cout << "\n";
        cout << "  Spatial grid of " << M << " by " << N << " points.\n";
        cout << "  The iteration will be repeated until the change is <= " << epsilon
             << "\n";

        // TODO: Get number of MPI processes
        cout << " Number of processes available = " << num_procs << "\n";
    }
    if (argc >= 2) {
        refFile = argv[1];
        std::cout << "Rank: \t" << rank << "  Reading reference results from " << refFile << std::endl;
    } else {
        std::cout << "Rank: \t" << rank << "  No reference file provided - skipping verification." << std::endl;
    }


    //
    //  Set the boundary values, which don't change.
    // TODO: Be aware that each process hold a different part of the matrix and may only
    // initialize his own part

    mean = 0.0;

    if(rank == 0)
    {
        for (i = 1; i < rows; i++) {
            w[i][0] = 100.0;
            w[i][N - 1] = 100.0;
        }
        for (j = 0; j < N; j++) {
            w[0][j] = 0.0;
        }
    }
    else if(rank == num_procs-1)
    {
        for (i = 0; i < rows - 1; i++) {
            w[i][N - 1] = 100.0;
            w[i][0] = 100.0;
        }
        for (j = 0; j < N; j++) {
            w[rows - 1][j] = 100.0;
        }
    }
    else
    {
        for (i = 0; i < rows; i++) {
            w[i][0] = 100.0;
            w[i][N - 1] = 100.0;
        }
    }

    //
    //  Average the boundary values, to come up with a reasonable
    //  initial value for the interior.
    //  TODO: Keep in mind that you need to average over all processes
    if(rank == 0)
    {

        for (j = 0; j < N; j++) {
            mean = mean + w[0][j];
        }
        for (i = 1; i < rows; i++) {
            mean = mean + w[i][0] + w[i][N - 1];
        }
    }
    else if(rank == num_procs-1)
    {
        for (j = 0; j < N; j++) {
            mean = mean + w[rows - 1][j];
        }
        for (i = 0; i < rows - 1; i++) {
            mean = mean + w[i][0] + w[i][N - 1];
        }
    }
    else
    {
        for (i = 0; i < rows; i++) {
            mean = mean + w[i][0] + w[i][N - 1];
        }
    }

    double out;
    MPI_Allreduce(&mean, &out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    mean = out;
    mean = mean / (double) (2 * M + 2 * N - 4);

    if(rank == 0)
    {
        cout << "\n";
        cout << "  MEAN = " << mean << "\n";
    }
    //
    //  Initialize the interior solution to the mean value.
    //
    if(rank == 0)
    {
        for (i = 1; i < rows; i++) {
            for (j = 1; j < N - 1; j++) {
                w[i][j] = mean;
            }
        }
    }
    else if(rank == num_procs-1)
    {
        for (i = 0; i < rows - 1; i++) {
            for (j = 1; j < N - 1; j++) {
                w[i][j] = mean;
            }
        }
    }
    else
    {
        for (i = 0; i < rows; i++) {
            for (j = 1; j < N - 1; j++) {
                w[i][j] = mean;
            }
        }
    }


    //
    //  Iterate until the  new solution W differs from the old solution U
    //  by no more than EPSILON.
    //
    iterations = 0;
    iterations_print = 1;
    if(rank == 0) {
        cout << "\n";
        cout << " Iteration  Change\n";
        cout << "\n";
        cout << "------------------" << std::endl;
    }
    Timer timer;

    diff = epsilon;

    // TODO: This statement should yield the same result on all processes, so that all
    // processes stop at the same iteration
    while (epsilon <= diff) {
        //
        //  Save the old solution in U.
        //
        for (i = 0; i < rows; i++) {
            for (j = 0; j < N; j++) {
                u[i][j] = w[i][j];
            }
        }
        //
        //  Determine the new estimate of the solution at the interior points.
        //  The new solution W is the average of north, south, east and west
        //  neighbors.
        //  TODO: Here you may need parts of the matrix that are part of other processes

        if(rank == 0)
            MPI_Isend(u[rows-1], N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request0);
        else if(rank == num_procs-1)
            MPI_Isend(u[0], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requestN);
        else
        {
            MPI_Isend(u[rows-1], N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &requestK2);
            MPI_Isend(u[0], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requestK);
        }


        //Do all interior nodes, which can be done without communication
        for(i = 1; i < rows-1; i++)
        {
            for (j = 1; j < N - 1; j++)
            {
                w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;
            }
        }

        if(rank == 0) //does only need row below and must send the lowest row only
        {
            double below[N];
            MPI_Irecv(below, N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request0);
            MPI_Wait(&request0, MPI_STATUS_IGNORE);
            //calc remaining row
            for (j = 1; j < N - 1; j++)
            {
                w[rows-1][j] = (u[rows-2][j] + below[j] + u[rows-1][j - 1] + u[rows-1][j + 1]) / 4.0;
            }

        }
        else if(rank == num_procs-1) //does only need row above and must send the highest row only
        {
            double above[N];
            MPI_Irecv(above, N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requestN);
            MPI_Wait(&requestN, MPI_STATUS_IGNORE);
            //calc remaining row
            for (j = 1; j < N - 1; j++)
            {
                w[0][j] = (above[j] + u[1][j] + +u[0][j - 1] + u[0][j + 1]) / 4.0;
            }
        }
        else //does need the row above and below and must send both highest and lowest row
        {
            double above[N];
            double below[N];
            MPI_Irecv(above, N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requestK);
            MPI_Irecv(below, N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &requestK2);
            MPI_Wait(&requestK2, MPI_STATUS_IGNORE);
            MPI_Wait(&requestK, MPI_STATUS_IGNORE);
            //calc remaining rows
            if(rows == 1)
            {
                for (j = 1; j < N - 1; j++)
                {
                    w[0][j] = ( above[j] + below[j] + + u[0][j - 1] + u[0][j + 1]) / 4.0;
                }
            }
            else
            {
                for (j = 1; j < N - 1; j++)
                {
                    w[0][j] = ( above[j] + u[1][j] + u[0][j - 1] + u[0][j + 1]) / 4.0;
                    w[rows-1][j] = (u[rows-2][j] + below[j] + u[rows-1][j - 1] + u[rows-1][j + 1]) / 4.0;
                }
            }

        }

        diff = 0.0;

        // TODO: Be aware that each process only computes its local diff here. You may
        // need to combine the results from all processes

        if(rank == 0)
        {
            for (i = 1; i < rows; i++) {
                for (j = 1; j < N - 1; j++) {
                    if (diff < fabs(w[i][j] - u[i][j])) {
                        diff = fabs(w[i][j] - u[i][j]);
                    }
                }
            }
        }
        else if(rank == num_procs-1)
        {
            for (i = 0; i < rows - 1; i++) {
                for (j = 1; j < N - 1; j++) {
                    if (diff < fabs(w[i][j] - u[i][j])) {
                        diff = fabs(w[i][j] - u[i][j]);
                    }
                }
            }
        }
        else
        {
            for (i = 0; i < rows; i++) {
                for (j = 1; j < N - 1; j++) {
                    if (diff < fabs(w[i][j] - u[i][j])) {
                        diff = fabs(w[i][j] - u[i][j]);
                    }
                }
            }
        }

        double rec_diff;
        MPI_Allreduce(&diff, &rec_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        diff = rec_diff;
        iterations++;
        if(rank == 0){
            if (iterations == iterations_print) {
                cout << "  " << setw(8) << iterations << "  " << diff << "\n";
                iterations_print = 2 * iterations_print;
            }
        }

    }
    // TODO: Insert a Barrier before time measurement
    MPI_Barrier(MPI_COMM_WORLD);
    wtime = timer.getSeconds();
    // TODO: Make sure that only one MPI process prints the output
    if(rank == 0)
    {
        cout << "\n";
        cout << "  " << setw(8) << iterations << "  " << diff << "\n";
        cout << "\n";
        cout << "  Error tolerance achieved.\n";
        cout << "  Wallclock time = " << setprecision(3) << wtime << "s\n";
        //
        //  Terminate.
        //
        cout << "\n";
        cout << "HEATED_PLATE Asynchronous:\n";
        cout << "  Normal end of execution.\n";
    }
    if (!refFile.empty()) {
        double ref[M][N];
        if (!readFromFile(refFile, &ref[0][0], M * N)) {
            std::cout << "Rank: \t" << rank << "  Verification failed - could not load reference file." << std::endl;
            MPI_Finalize();
            return 0;
        }
        int index;
        //because of adding the remaining rows prior,
        // the index hast to be calculated by differentiate between processes which have been given an extra row and which haven't
        if(rank < M % num_procs)
        {
            index = rank * rows;
        }
        else
        {
            // M%num_procs          Equals the number of processes which got an extra row
            // rows+1               Equals the row size of processes which got an extra row
            // rank - (M%num_procs) Number of processes until this process without additional row
            index = (M % num_procs) * (rows+1) + rows * (rank - (M % num_procs) );
        }
        if (verify(&w[0][0], &ref[index][0], rows * N, rank)) {
            std::cout << "Rank: \t" << rank << "  Result is valid." << std::endl;
        } else {
            std::cout << "Rank: \t" << rank << "  Verification failed!" << std::endl;
        }
    }
    MPI_Finalize();
    return 0;

#undef M
#undef N
}