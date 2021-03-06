//++++++++++++++++++++++++++++filename: mygeneral.h +++++++++++++++++++++++++++++++++++//

/***************************************************************************************/
/* the head file provide basic variables type, which is convenient for whole calculation
/* all variables with respect to the length, their unit are m
/**************************************************************************************/
#ifndef _MYGENERAL_H
#define _MYGENERAL_H

#include <complex>
#include "mpi.h"
using namespace std;
typedef complex<double> Complex;

//public variables
typedef double (*pf)(double,double,double);
#define PI 3.1415926

struct Position{
	double x;
	double y;
	double z;
};
struct AirStruct{
	double Pav;//pa
	double adiabatic_coef;
	double Aav;//m3/kg
	double sound_speed;//m/s
};
struct PoreStruct{
	double Cs;//structure factor of porous medium
	double eff_density;//kg/m3;
	double porosity;
	double resistivity;
	double Kp;//m.s2/kg
	double q;
	double Z_inf;
	double tau;
};
struct DifferenceStep{
	double diff_t;//second
	double diff_x;//m
	double diff_y;//m
	double diff_z;//m
};

//enum medium_list {AirMedium,PoreMedium,Coupling,WestMedium,NorthMedium,SouthMedium,
//EastMedium,WSMedium,WNMedium};//coupling is the interface between air and porous
enum scheme_list {FB_p_v,FB_v_p,FB_vp,LeapTrap,ABM,RK4,RK2};
enum boundary_type{rigid,absorbing,porous_media,radiation,TDBC};
enum boundary_location{NorthBC,SouthBC,WestBC,EastBC,FrontBC,BackBC};

struct boundary{
	boundary_type air_west;
	boundary_type air_east;
	boundary_type air_south;
	boundary_type air_north;
	boundary_type air_front;
	boundary_type air_back;
};

struct MovingFrame{
	int Judge;
	int JMAX;
	int lead_DJ;
	int trail_DJ;
	double ratio;
};

struct hillpore{
	int curve_judge;
	int curve_points;
	double curve_coefa;
	double curve_coefb;
	double curve_coefc;
	double curve_coefd;
	double curve_coefe;
	double curve_coeff;
	double curve_coefg;
	double curve_coefh;
	double curve_coefi;
	double curve_coefj;
	double curve_coefk;
	double curve_coefl;
	double Cs;
	double eff_density;
	double porosity;
	double resistivity;
	double Kp;
	double q;
	double curve_xstar;
	double curve_xstop;
	double curve_ystar;
	double curve_ystop;
	double curve_zstar;
	double curve_zstop;
};
// this is material known for density and speed of sound
struct material{
	double density;
	double speed;
};
struct PML_boundary{
	int judge;
	double alpha;
	double resistivity;
	int time_step;
};
struct MPI_info{
	int mpi_rank;
	int mpi_size;
	int mpi_xarea;
	int mpi_yarea;
	int mpi_zarea;
	//MPI::Cartcomm mpi_comm3d;
	MPI_Comm mpi_comm3d;
	int mpi_nleft;
	int mpi_nright;
	int mpi_nbottom;
	int mpi_ntop;
	int mpi_nfront;
	int mpi_nback;
	int *mpi_coords;
};
//public function
void FFT(short int dir,long m,double x[],double y[]);
void FFT_output(long m,double x[],double y[],char* filename,double diff_t);
double Geo_over_var(double geo_var1,double geo_var2,double geo_var3,
					double geo_var4,double geo_jaco);
double Velo_over_var(double geo_var1,double geo_var2,double geo_var3,
			double Uav_up, double Uav_down,double Vav_up,double Vav_down,
			double Wav_up, double Wav_down,double geo_x,double geo_y,double geo_z);
double Coef_velo(double geo_var1,double geo_var2,double geo_var3,
						double Uav_vel,double Vav_vel,double Wav_vel,
						double geo_xt,double geo_yt,double geo_zt);
void mpe_decomp1d(int n,int numprocs,int myid,int &s,int &e);
int qu_min(int a,int b);


/***************************************************************************************/
/*
	0.6 //gaussian_width
	501 2 2 //Moving domain:MF_IMAX; Moving distance: MF_I
	2 //scheme:1.FB_v_p 2.FB_p_v 3.FB_vp 4.LeapTrap 5.RK4 6.RK2 7.ABM ;
	2 3 2 1 //NorthBC SouthBC WestBC EastBC and
			//B.C type: 1.rigid 2.absorbing layer 3.porous media 4.radiation boundary condition;
	3950 50 17 300 //time_domain;out_difft(output each this value);FFT_m(2^FFT_m=N);frequency;

	1.0 5 //b(Vav=b*z);
	1.4 1.0e5 340 //AirPara.adiabatic_coef;AirPara.Pav;AirPara.sound_speed;
	1.0 1.0 500 //AbsorbPara.Cs;AbsorbPara.porosity;AbsorbPara.resistivity;
	3.0 0.3 2.0e5 //PorePara.Cs;PorePara.porosity;PorePara.resistivity;
	0 2 0 2 //airpore: source.y;source.z;receiver.y;receiver.z;

	0.5 0.1 0.1 //AirStep.diff_t;AirStep.diff_y;AirStep.diff_z;
	0.5 0.1 0.1 //TransientStep2.diff_t;TransientStep2.diff_y;TransientStep2.diff_z;
	0.5 0.1 0.025 //TransientStep2.diff_t;TransientStep2.diff_y;TransientStep2.diff_z;
	0 50 30 2 //AirFrame.left;AirFrame.right;AirFrame.upper;AirFrame.lower;
	0 50 2 2 //TransientFrame2.left;TransientFrame2.right;TransientFrame2.upper;TransientFrame2.lower;
	0 50 2 0 //TransientFrame1.left;TransientFrame1.right;TransientFrame1.upper;TransientFrame1.lower;

	0 50 40 30 //NorthFrame.left;NorthFrame.right;NorthFrame.upper;NorthFrame.lower;
	0 50 0 -10 //SouthFrame.left;SouthFrame.right;SouthFrame.upper;SouthFrame.lower;
	-10 0 30 2 //WestFrame.left;WestFrame.right;WestFrame.upper;WestFrame.lower;
	0 0 0 0 //WestSouthFrame.left;WestSouthFrame.right;WestSouthFrame.upper;WestSouthFrame.lower;
	-10 0 40 30 //WestNorthFrame.left;WestNorthFrame.right;WestNorthFrame.upper;WestNorthFrame.lower;

	porous_media;
	0 50 -2 -10 //SouthFrame.left;SouthFrame.right;SouthFrame.upper;SouthFrame.lower;
	0.5 0.1 0.1 //SouthStep.diff_t;SouthStep.diff_y;SouthStep.diff_z;
	7 8 4.0 0.5	  //EastFrame.left;EastFrame.right;EastFrame.upper;EastFrame.lower;
	0.15e-3 0.1 0.1 //EastFrame.difft;EastFrame.diff_y;EastFrame.diff_z;

*/
/*
description of input.txt
Moving frame=0 denoting calculation will be performed only right boundary moving
Moving frame=1 denoting the calculation will be performed with moving frame (both left and right
boundary move). Besides, under the condition, if MF_IMAX=whole_IMAX, no boundary moves.
if the input AirStep.diff_t < 0.01, the input represents n for the formula: diff_t=n*diff_z/c.
n may be 0.5, 0.25 etc. Thus(move_frame.trail_DI-1)*AirStep.diff_y/AirPara.sound_speed/AirStep.diff_t
is competely integer without any approximation;
*/
/***************************************************************************************/

#endif
