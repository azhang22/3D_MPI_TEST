//++++++++++++++++++++++++++++++filename: air.h +++++++++++++++++++++++++++++++++++++++//

//----------------------------description of variable----------------------------------//
/***************************************************************************************/ 
/* Vav,Wav is average speed. Vav_over_y means patial (Vav)/patial (y)
/* Pav is average atmospheric pressure
/* Aav is a reversal of average density
/* adiabatic_coef is adiabatic coefficient
/* Profile_Vav,Profile_Wav point to the average wind speed profile
/**************************************************************************************/ 
#ifndef _AIR_H
#define _AIR_H
#include "calculation_2D.h"

//**************************wind speed profile*****************************//
//z is with respect to local ground level
class SpeedProfile{
private:
	int velocity_method;
	double y0,z0,L,alpha;//m
	double b,c,circulation,vorticity;//b is Mach number for vortex back flow
public:
	SpeedProfile(){};
	SpeedProfile(double,double,int);
	~SpeedProfile();
	double Uav(double x,double y,double z);
	double Vav(double x,double y,double z);
	double Wav(double x,double y,double z);
};

class air:public calculation_2D
{
	protected:
		double ***Uav,***Vav,***Wav;
		double Uav_over_X,Uav_over_Y,Uav_over_Z,Vav_over_X,Vav_over_Y;
		double Vav_over_Z,Wav_over_X,Wav_over_Y,Wav_over_Z;
		double coef1,coef2,coef3;
		double Pav,Aav,adiabatic_coef,sound_speed;
		int velocity_method;
		SpeedProfile *SP;
		
		double ***Qux1,***Qux2,***Qvx1,***Qvx2,***Qwx1,***Qwx2,***Qpx1,***Qpx2;
		double ***Quy1,***Quy2,***Qvy1,***Qvy2,***Qwy1,***Qwy2,***Qpy1,***Qpy2;
		double ***Quz1,***Quz2,***Qvz1,***Qvz2,***Qwz1,***Qwz2,***Qpz1,***Qpz2;
		int PML_judgez,PML_judgey,PML_judgex;
		// PML boundary prameters for the z direction
		double PML_AbsorbZmax,PML_alphaz;
		double PML_Height1,PML_z1,PML_Height2,PML_z2;
		// PML boundary prameters for the y direction
		double PML_AbsorbYmax,PML_alphay;
		double PML_Length1,PML_y1,PML_Length2,PML_y2;
		// PML boundary prameters for the x direction
		double PML_AbsorbXmax,PML_alphax;
		double PML_Width1,PML_x1,PML_Width2,PML_x2;
		int PML_IMAX1,PML_IMAX2,PML_JMAX1,PML_JMAX2,PML_JMAX1_time,PML_KMAX1,PML_KMAX2;

		// for the porous for the hill
		double *cr_y,*cr_z;
		int ***cr_uc;
		int cr_judge,cr_points;
		double cr_a,cr_b,cr_c,cr_d,cr_e,cr_f,cr_g,cr_h,cr_i,cr_j,cr_k,cr_l;
		double crpo_eff_density,crpo_porosity,crpo_resistivity;
		double crpo_eff_density1,crpo_porosity1,crpo_resistivity1;
		double crpo_eff_density2,crpo_porosity2,crpo_resistivity2;
		double crpo_Kp,crpo_Beta,crpo_Gama,crpo_invBeta;
		double crpo_Kp1,crpo_Beta1,crpo_Gama1,crpo_invBeta1;
		double crpo_Kp2,crpo_Beta2,crpo_Gama2,crpo_invBeta2;
		int istar,istop,jstar,jstop,kstar,kstop;

	public:
		air(){};
		air(scheme_list,DifferenceStep,char *coordi,MovingFrame MF,const double,
		double velocity_coef,int velocity_method1,AirStruct AirPara,PML_boundary PML_xbd,
			PML_boundary PML_ybd,PML_boundary PML_zbd,hillpore hill,hillpore canopy,hillpore trunk,
			MPI_info mpi_var,int judge_porous, int restart_timer);
		~air();
		virtual void cal_fuvw(double ***temp_u,double ***temp_v,double ***temp_w,
			double ***temp_p);
		virtual void cal_fp(double ***temp_u,double ***temp_v,double ***temp_w,
			double ***temp_p);
		virtual void save_restart_air(char *restartfile);
		virtual void input_restart_air(char *restartfile,int *point_pois);
		virtual void SetWindProfile(int N);
		
		virtual void cal_hill_var(hillpore,hillpore,hillpore);
		virtual void cal_PML_var(PML_boundary,PML_boundary,PML_boundary);
		virtual void Update_PML_Qvw();
		virtual void Update_PML_Qp();
		virtual void SetQMove(int N);
		virtual void Set_Timer(int N);
			
};

#endif
