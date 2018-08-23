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
	MPI::Init(argc,argv);
	mpi_rank=MPI::COMM_WORLD.Get_rank();
	mpi_size=MPI::COMM_WORLD.Get_size();
	strcpy(infile,"input.txt");
	airpore *AirPore;
	AirPore=new airpore(infile);
//	goto loop;
	AirPore->get_output(mpi_rank,mpi_size);
//	return 0;
//	loop:
	if (mpi_rank==0) AirPore->get_FFT_y1(2);// for the long-distance sound propagation
	delete AirPore;
	MPI::Finalize();
	return 0;
}
