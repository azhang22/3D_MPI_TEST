//+++++++++++++++++++++++++++filename: porous.cpp ++++++++++++++++++++++++++++++//

//-----------------------porous media acoustic wave equation---------------------//
/*******************************************************************************/ 
/* the program just apply Navier-Stokes equation to calculate acoustic pressure 
/* when sound propagation in the ground (porous media) with Eulerian time-domain model 
/*******************************************************************************/ 

//-----------------------scheme from Victor W.Sparrow--------------------------//
/*******************************************************************************/ 
/* Victor W.Sparrow
/* calculate pressure first, and then calculate velocity
/* apply non-staggered grid, velocity and pressure are at the same grid point
/* point pieces at y direction(air): 0,1,...IMAX-1;
/* point pieces at z direction(air): 0,1,...JMAX-1;
/********************************************************************************/

//---------------------------scheme from Erik M.Salomons-----------------------//
/*******************************************************************************/ 
/* Eulerian Time-Domain Model for sound propagation over 
/* a finite impedance ground surface. 
/* Comparison with frequency-domain models. 
/* Author: Erik M.Salomons.
/* ACTA Vol.88(2002) 483-492
/* calculate velocity first, and then calculate pressure
/* apply staggered grid, velocity and pressure are at the different grid point
/* point pieces at y direction(air): 0,1,...IMAX-1;
/* point pieces at z direction(air): 0,1,...JMAX-1;
/* east and west boundary point of v is v[1][j] and 
/* north and south boundary point of w is w[i][1] and w[i][JMAX-1]
/* the four boundary point v[1][j], v[IMAX-1][j], w[i][1], w[i][JMAX-1]
/* are calculated through equation, not from interpolation. This is different from
/* above colocated scheme
/*******************************************************************************/ 
#include <stdio.h>
#include <math.h>
#include "porous.h"


porous::porous(scheme_list scheme1,DifferenceStep Step,char *coordi,MovingFrame MF,
			   const double gauss_width,PoreStruct PorePara,MPI_info mpi_var,int judge_porous)
			   :calculation_2D(scheme1,Step,coordi,MF,gauss_width,mpi_var,judge_porous)
{
	eff_density=1.0/PorePara.eff_density;
	porosity=PorePara.porosity;
	resistivity=PorePara.resistivity;
	Kp=PorePara.Kp;
	Beta=resistivity*diff_t*eff_density;
	Gama=diff_t/Kp/porosity;
	invBeta=1.0/(Beta+1.0);
}

porous::~porous()
{
}

//----------------------------calculate velocity and pressure------------------//
void porous::cal_fuvw(double ***temp_u,double ***temp_v,double ***temp_w,double ***temp_p)
{
	int i,j,k;
	for(k=1;k<mpi_KMAX-1;k++){
		for(j=1;j<mpi_JMAX-1;j++){
			for(i=1;i<mpi_IMAX-1;i++){
				x_over_X=Geo_over_var(Y_over_y[j],Z_over_z[k],Y_over_z,
									Z_over_y,Jacobian[i][j][k]);
				x_over_Y=Geo_over_var(X_over_z,Z_over_y,X_over_y,
									Z_over_z[k],Jacobian[i][j][k]);
				x_over_Z=Geo_over_var(X_over_y,Y_over_z,X_over_z,
									Y_over_y[j],Jacobian[i][j][k]);
				y_over_X=Geo_over_var(Y_over_z,Z_over_x,Y_over_x,
									Z_over_z[k],Jacobian[i][j][k]);
				y_over_Y=Geo_over_var(X_over_x[i],Z_over_z[k],X_over_z,
									Z_over_x,Jacobian[i][j][k]);
				y_over_Z=Geo_over_var(X_over_z,Y_over_x,X_over_x[i],
									Y_over_z,Jacobian[i][j][k]);
				z_over_X=Geo_over_var(Y_over_x,Z_over_y,Y_over_y[j],
									Z_over_x,Jacobian[i][j][k]);
				z_over_Y=Geo_over_var(X_over_y,Z_over_x,X_over_x[i],
									Z_over_y,Jacobian[i][j][k]);
				z_over_Z=Geo_over_var(X_over_x[i],Y_over_y[j],X_over_y,
									Y_over_x,Jacobian[i][j][k]);

				u_nn[i][j][k]=-1.0*invBeta*diff_t*eff_density*
										(x_over_X*(temp_p[i][j][k]-temp_p[i-1][j][k])*diff_x
										+y_over_X*(temp_p[i][j][k]-temp_p[i][j-1][k])*diff_y
										+z_over_X*(temp_p[i][j][k]-temp_p[i][j][k-1])*diff_z)+
										temp_u[i][j][k]*invBeta-temp_u[i][j][k];
				v_nn[i][j][k]=-1.0*invBeta*diff_t*eff_density*
										(x_over_Y*(temp_p[i][j][k]-temp_p[i-1][j][k])*diff_x
										+y_over_Y*(temp_p[i][j][k]-temp_p[i][j-1][k])*diff_y
										+z_over_Y*(temp_p[i][j][k]-temp_p[i][j][k-1])*diff_z)+
										temp_v[i][j][k]*invBeta-temp_v[i][j][k];

				w_nn[i][j][k]=-1.0*invBeta*diff_t*eff_density*
										(x_over_Z*(temp_p[i][j][k]-temp_p[i-1][j][k])*diff_x
										+y_over_Z*(temp_p[i][j][k]-temp_p[i][j-1][k])*diff_y
										+z_over_Z*(temp_p[i][j][k]-temp_p[i][j][k-1])*diff_z)+
										temp_w[i][j][k]*invBeta-temp_w[i][j][k];


			}
		} 
	}
}


void porous::cal_fp(double ***temp_u,double ***temp_v,double ***temp_w,double ***temp_p)
{
	int i,j,k;
	for(k=1;k<mpi_KMAX-1;k++){
		for(j=1;j<mpi_JMAX-1;j++){
			for(i=1;i<mpi_IMAX-1;i++){
				x_over_X=Geo_over_var(Y_over_y[j],Z_over_z[k],Y_over_z,
									Z_over_y,Jacobian[i][j][k]);
				x_over_Y=Geo_over_var(X_over_z,Z_over_y,X_over_y,
									Z_over_z[k],Jacobian[i][j][k]);
				x_over_Z=Geo_over_var(X_over_y,Y_over_z,X_over_z,
									Y_over_y[j],Jacobian[i][j][k]);
				y_over_X=Geo_over_var(Y_over_z,Z_over_x,Y_over_x,
									Z_over_z[k],Jacobian[i][j][k]);
				y_over_Y=Geo_over_var(X_over_x[i],Z_over_z[k],X_over_z,
									Z_over_x,Jacobian[i][j][k]);
				y_over_Z=Geo_over_var(X_over_z,Y_over_x,X_over_x[i],
									Y_over_z,Jacobian[i][j][k]);
				z_over_X=Geo_over_var(Y_over_x,Z_over_y,Y_over_y[j],
									Z_over_x,Jacobian[i][j][k]);
				z_over_Y=Geo_over_var(X_over_y,Z_over_x,X_over_x[i],
									Z_over_y,Jacobian[i][j][k]);
				z_over_Z=Geo_over_var(X_over_x[i],Y_over_y[j],X_over_y,
									Y_over_x,Jacobian[i][j][k]);
				p_nn[i][j][k]=-Gama*(x_over_X*diff_x*(temp_u[i+1][j][k]-temp_u[i][j][k])+
									y_over_X*diff_y*(temp_u[i][j+1][k]-temp_u[i][j][k])+
									z_over_X*diff_z*(temp_u[i][j][k+1]-temp_u[i][j][k])+
									x_over_Y*diff_x*(temp_v[i+1][j][k]-temp_v[i][j][k])+
									y_over_Y*diff_y*(temp_v[i][j+1][k]-temp_v[i][j][k])+
									z_over_Y*diff_z*(temp_v[i][j][k+1]-temp_v[i][j][k])+
									x_over_Z*diff_x*(temp_w[i+1][j][k]-temp_w[i][j][k])+
									y_over_Z*diff_y*(temp_w[i][j+1][k]-temp_w[i][j][k])+
									z_over_Z*diff_z*(temp_w[i][j][k+1]-temp_w[i][j][k]));
			}
		}
	}
}

