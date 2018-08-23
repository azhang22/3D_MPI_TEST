	//+++++++++++++++++++++++++++filename: main.cpp +++++++++++++++++++++++++++++++//
#include <stdlib.h>
#include <string.h>
#include "airpore.h"
#include "math.h"
#include "mpi.h"

int main(int argc,char *argv[])
{
	int mpi_rank,mpi_size;
	char infile[100]=" ";
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	strcpy_s(infile,"input.txt");
	airpore *AirPore;
	AirPore=new airpore(infile);
//	goto loop;
	AirPore->get_output(mpi_rank,mpi_size);
//	return 0;
//	loop:
	if (mpi_rank==0) AirPore->get_FFT_y1(2);// for the long-distance sound propagation
	delete AirPore;
	MPI_Finalize();
	return 0;
}
