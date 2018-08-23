//++++++++++++++++++++++++++++++filename: airpore.h +++++++++++++++++++++++++++++++++++++++//

//----------------------------description of variable----------------------------------//
/***************************************************************************************/
/* the number of points FFT_N be a power of 2 (2^FFT_m)
/* p_t_filename and p_f_filename is filename of pressure output within time-domain and
/* frequency domain respectively.
/* AirStruct is struct data containing air structure parameter such as density.
/* it is the same for PoreStruct for porous media
/* pr,pi record real and imaginary pressure at receiver
/*************************************************************************************/
#ifndef _AIRPORE_H
#define _AIRPORE_H
#include "mygeneral.h"

class calculation_2D;

class airpore
{
private:
	int time_domain,FFT_m,FFT_N,out_difft;
	int CaseNo;
	bool MovingFrame;
	double gauss_width,frequency;
	scheme_list scheme;
	struct MovingFrame move_frame;
	double velocity_coef;
	int velocity_method;
	int restart,restart_out,Rec_count,restart_timer;
	char restart_infile[20],restart_infile1[200];

	calculation_2D *AirMedia;
	calculation_2D *west_pore,*east_pore,*south_pore,*north_pore,*front_pore,*back_pore;

	boundary air_boundary;
	AirStruct AirPara;
	PoreStruct PorePara;
	PML_boundary PML_xbd,PML_ybd,PML_zbd;
	Position source,receiver;

	DifferenceStep AirStep,SouthStep;
	hillpore hill,canopy,trunk; // add the porous media for hill;

	// define mpi_rank,mpi_size at the airpore.cpp
	MPI_info mpi_var;
	int mpi_rank,mpi_size;
public:
	airpore();
	airpore(char *input_file);
	~airpore();
	void get_output(int mpi_rank1,int mpi_size1);
	int Read_restart();
	void Save_restart(int index);
	void save_restartfile(char *restartfile);
	void input_restartfile(int *point_pois);
	void MPI_Initialize();
	void SetInitialCond();
	void SetMovingFrame(int MF_count);
	void Cal_pressure();
	void Cal_velocity();
	void UpdateInitialCond(int MF_count);
	void UpdateSource(double f0,int N, Position);
	void mpisend_data_contour(int MF_count);
	void get_data_contour(int n);
	void get_FFT_y1(int out_type);
	void copyFile(const char* src_loc, const char* dest_loc);
};

#endif

