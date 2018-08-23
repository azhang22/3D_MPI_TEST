//++++++++++++++++++++++++++++++filename: calculation_2D.h +++++++++++++++++++++++++++++++++++++++//

//----------------------------description of variable----------------------------------//
/***************************************************************************************/
/* v_n,w_n,p_n will store velocity of y direction, velocity of z direction and pressure
/* at old time v_nn,w_nn,p_nn will store velocity of y direction, velocity of z direction
/* and pressure at new time diff_t_air,diff_y_air,diff_z_air is time step, space step
/* of y direction, space step of z direction Vav_over_y,Vav_over_z,Wav_over_y,Wav_over_z
/* is differential of average wind speed to axial y and z
/* IMAX,JMAX the max number of points for y and z direction within moving frame
/* IMAX=frame_JMAX reqiure the number is even number so as to be divided as two part when computing
/* whole_IMAX and whole_JMAX are max number of points for whole domain
/* Y_over_y is partial differential of Y to y, where Y is coordinate within physical plane
/* y is coordinate within computational plane.
/* Y_over_z, Z_over_y and Z_over_z have the same meaning as Y_over_y
/*************************************************************************************/
#ifndef _CALCULATION_2D_H
#define _CALCULATION_2D_H
#include "mygeneral.h"

#include <iostream>
#include <fstream>
using namespace std;

class calculation_2D
{
	protected:
		int IMAX,IMAX1,IMAX2,JMAX,JMAX1,JMAX2,JMAX3,KMAX,KMAX1,KMAX2,IJKMAX;
		double diff_t,diff_x,diff_y,diff_z;//time step and uniform spacing steps
		double diff_x1,diff_y1,diff_z1;
		scheme_list scheme;
		//enum Col_Method {FB1,FB2,LeapTrap,ABM,RK4,RK2} MyMethod;
		double ***u_nn,***v_nn,***w_nn,***p_nn;
		double ***u_n,***v_n,***w_n,***p_n;
		double *X,*Y,*Z;
		double *X_over_x,*Y_over_y,*Z_over_z;

		double X_over_y,X_over_z;
		double Y_over_x,Y_over_z;
		double Z_over_x,Z_over_y;
		double ***Jacobian;
		double x_over_X,x_over_Y,x_over_Z;
		double y_over_X,y_over_Y,y_over_Z;
		double z_over_X, z_over_Y,z_over_Z;
		double gaussian_coef;
		//definition for boundary values
		double *axisbound_X,*axisbound_Y,*axisbound_Z;
		int judge_porous1;
		// output data file for pressure contour
		int IMAX_out,JMAX_out,KMAX_out,out_n;
		double *X_out,*Y_out,*Z_out;
		double *p_in,*p_out;
		//this is added for moving frame
		struct MovingFrame move_frame;

		//define variables for MPI
		int mpi_rank,mpi_size;
		int mpi_xarea,mpi_yarea,mpi_zarea;//the number of subdomains in each direction
		int *mpi_coords; // the coordinates of each domain
		int mpi_nleft,mpi_nright,mpi_nbottom,mpi_ntop,mpi_nfront,mpi_nback;
		int mpi_i1,mpi_i2,mpi_j1,mpi_j2,mpi_k1,mpi_k2;
		int mpi_IMAX,mpi_IMAX1,mpi_IMAX2,mpi_JMAX,mpi_JMAX1,mpi_JMAX2,mpi_KMAX,mpi_KMAX1,mpi_KMAX2;
		//MPI::Cartcomm mpi_comm3d;
		//MPI::Status status;
		MPI_Comm mpi_comm3d;
		MPI_Status status;
		// save the data of sub domain at the boundary for moving frame
		double *uer,*uws,*ver,*vws,*wer,*wws,*per,*pws;
		double *ps;

	public:
		calculation_2D();
		calculation_2D(scheme_list,DifferenceStep,char *coordi,MovingFrame MF,
				const double,MPI_info mpi_var,int judge_porous);
		~calculation_2D();
		double cal_InitialPressure(double x,double y,double z);
		void Set_InitialCond(Position source,int judge);
		double get_pressure(int ii,int jj,int kk,int N);
		int get_whole_IMAX(){return IMAX;}
		int get_whole_JMAX(){return IJKMAX;}
		int get_whole_KMAX(){return KMAX;}
		void get_position(int &ii,int &jj,int &kk, Position receiver);
		virtual void cal_PML_var(PML_boundary,PML_boundary,PML_boundary){};
		virtual void cal_hill_var(hillpore){};
		virtual void Update_PML_Qvw(){};
		virtual void Update_PML_Qp(){};
		void UpdateInitialCond(int N);
		void UpdateSource(double f0, int N, Position source);  
		void save_restart_cal(char *restartfile);
		void input_restart_cal(char *restart_file,int*point_pois);
		void SetMovingFrame(int N);
		void mpi_send_data(int N);

		virtual void SetWindProfile(int N){};
		virtual void SetQMove(int N){};
		virtual void Set_Timer(int N){};
		virtual void save_restart_air(char *restartfile){};
		virtual void save_restart_pml(char *restartfile){};
		virtual void input_restart_air(char *restart_infile,int *point_pois){};
		virtual void input_restart_pml(char *restart_infile,int *point_pois){};

		void Cal_velocity(void);
		void Cal_pressure(void);
		void MPIexchange_pressure();
		void MPIexchange_velocity();
		void UpdateBC_pressure(boundary_location BC);
		void UpdateBC_velocity(boundary_location BC);
		void UpdateBC_velocity(boundary_location BC,calculation_2D& porous);
		void UpdateBC_pressure(boundary_location BC,calculation_2D& porous);

		virtual void cal_fuvw(double ***temp_u,double ***temp_v,double ***temp_w,
			double ***temp_p){};
		virtual void cal_fp(double ***temp_u,double ***temp_v,double ***temp_w,
			double ***temp_p){};
		
		friend ostream& operator<<(ostream& output_stream, calculation_2D& my_object);
		friend void getOstream(ostream& output_stream, calculation_2D& my_object);
		friend ostream& coordi(ostream& output_stream,calculation_2D& my_object);
};

#endif
