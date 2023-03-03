//  Author:
//
//    Yannic Fischler
//    Modfied for WS21 by Sebastian Kreutzer
//    Completed by Group 3

#include <stdlib.h>
#include <stdio.h>
#include <time.h>  // für die Zufallszahlen-Seeds
#include <mpi.h>


int world_size;
int world_rank;
double sumSoFar;
double mySum = 0;  // mySum = Summe aller y-Koordinaten in diesem Prozess 

// Die folgende Funktion berechnet Zufallszahlen zwischen low und high.
double getRandomXValue( double low, double high )
{
  return low + ((double)rand() / (double)RAND_MAX) * ( high - low );
}

// Die Funktion f sowie die untere und obere Grenze des Integrals werden spezifiert.
double low = -1;
double high = 1;
double f(double x) {
  return 2 / (1 + x*x);
}

// Unsere Alternative zu MPI_Reduce für die gegebene Aufgabe ohne kollektive Operationen
double myReduceFunction() {
  if (world_rank != world_size-1) {
    MPI_Recv(&sumSoFar, 1, MPI_DOUBLE, world_rank+1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    sumSoFar += mySum;
  } else {
    sumSoFar = mySum;
  }

  if (world_rank != 0) {
    MPI_Send(&sumSoFar, 1, MPI_DOUBLE, world_rank-1, 42, MPI_COMM_WORLD);
  }
}

// Beginn der main-Methode
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  // DONE: Store number of processes in world_size
  // DONE: Store current rank in world_rank
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

  if (argc != 2) {
    printf("Add number of sampling points as parameter\n");
    return 1;
  }

  int numberOfPoints = atoi(argv[1]);
  
  // DONE: Make sure that numberOfPoints is valid
  if (numberOfPoints <= 0 && world_rank == 0) {
    printf ("Number of points must be positive!");
    return 1;
  }
  // Falls weniger Punkte angefordert werden als Prozesse gibt, werden stattdessen genau so viele Punkte simuliert wie es Prozesse 
  // gibt, da andernfalls einige Prozesse überhaupt keine Arbeit hätten, vgl. Fußnote 2 in Aufgabenstellung
  if (numberOfPoints < world_size) {
    numberOfPoints = world_size;
  }

  if (world_rank == 0) {
    printf("Running with %d processes\n", world_size);
  }

  double result = 0;

  // TODO: Implement the Monte-Carlo simulation as described in the task.
  //       Store the solution in "result".

  // x ist das ganzzahlige Ergebnis der Division numberOfPoints/world_size, y der Rest, d.h. es gilt x * world_size + y = numberOfPoints.
  // Die ersten y Prozesse kümmern sich bei uns um x+1 Punkte, die anderen worldsize-y Prozesse um x Punkte.
  // Insgesamt simulieren die Prozesse also y*(x+1)+(worldsize-y)*x = y + world_size * x = numberOfPoints Punkte, wie gefordert. Außerdem simulieren sie
  // allesamt ähnlich viele Punkte (x bzw. x+1 Punkte) und brauchen somit vermutlich auch ähnlich lang, was der Perfomance zugute kommt.
  int x = numberOfPoints / world_size;
  int y = numberOfPoints % world_size;
  int myPoints;    // myPoints = Anzahl der von diesem Prozess simulierten Punkte
  if (world_rank <= y-1) {
    myPoints = x+1;
  } else {
    myPoints = x;
  }
  
  srand(world_rank+time(NULL));  // Für jeden Prozess wird ein separater Zufallszahlen-Seed erstellt.
  double randomXValue, yValue;
  for (int i=0; i < myPoints; i++) {
    randomXValue = getRandomXValue(low, high);
    //printf ("%f -> ", randomXValue);
    yValue = f(randomXValue);
    //printf ("%f | ", yValue);
    mySum += yValue;
  }

  //MPI_Reduce(&mySum, &completeSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  myReduceFunction();

  if (world_rank == 0) {
    double completeSum = sumSoFar;  // completeSum = Summe aller y-Koordinaten in allen Prozessen
    result = ((high - low) / numberOfPoints) * completeSum; 
    printf("%f\n", result);
  }

  MPI_Finalize();

  return 0;
}






