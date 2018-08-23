//++++++++++++++++++++++++++filename: air.cpp +++++++++++++++++++++++++++++++++//

//-----------------------air acoustic wave equation--------------------------//
/*******************************************************************************/
/* the program just apply Navier-Stokes equation to calculate acoustic pressure
/* when sound propagation with Eulerian time-domain model
/*******************************************************************************/

//-----------------------scheme from Victor W.Sparrow--------------------------//
/*******************************************************************************/
/* Victor W.Sparrow
/* calculate pressure first, and then calculate velocity
/* apply non-staggered grid(colocated), velocity and pressure are at the same grid point
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
/* east and west boundary point of v is v[1][j] and v[IMAX-1][j]
/* north and south boundary point of w is w[i][1] and w[i][JMAX-1]
/* the four boundary point v[1][j], v[IMAX-1][j], w[i][1], w[i][JMAX-1]
/* are calculated through equation, not from interpolation. This is different from
/* above colocated scheme
/*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "air.h"
#include "mygeneral.h"

SpeedProfile::SpeedProfile(double velocity_coef, double sound_speed, int VM)
{
	alpha = 1.256431;
	b = velocity_coef; c = sound_speed; velocity_method = VM;
	/*
	switch(velocity_method){
	case 2:
	{
	y0=0;z0=0;L=0.03;
	vorticity=b*(1+2*alpha)*c/L;
	circulation=vorticity*L*L*PI/alpha;
	}
	break;
	case 3:
	{
	y0=0.66;z0=0;L=0.03;
	circulation=b;
	}
	break;
	case 6:
	{
	y0=3.5;z0=4;L=1;
	vorticity=b*(1+2*alpha)*c/L;
	circulation=vorticity*L*L*PI/alpha;
	}
	break;
	default:
	;
	}
	*/
}
double SpeedProfile::Uav(double x, double y, double z)
{
	switch (velocity_method){
	case 1:
	{
		return 0;
	}
	break;
	case 2:
	{
		return 0;
		/*
					double r,Vr;
					r=sqrt(pow((y-y0),2)+pow((z-z0),2));
					if((r/L)<=0.000001) r=1e-10;
					Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
					return -(y-y0)/r*Vr;
					*/
	}
	break;
	case 3:
	{
		double r, Vr;
		r = sqrt(pow((y - y0), 2) + pow((z - z0), 2));
		if (r <= L){//Vr=circulation/(2*PI*L)*r/L;Wav=-(y-y0)/r*Vr;
			return circulation / (2 * PI*L) / L*(-(y - y0));
		}
		Vr = circulation / (2 * PI) / r;
		return -(y - y0) / r*Vr;
	}
	break;
	case 4:
	{
		return 0;
	}
	break;
	case 5:
	{
		return b*y;
	}
	break;
	case 6:
	{
		double r, Vr;
		r = sqrt(pow((y - y0), 2) + pow((z - z0), 2));
		if ((r / L) <= 0.00001) return 0;
		Vr = circulation / (2 * PI*r)*(1 - exp(-alpha*r*r / L / L));
		return -(y - y0) / r*Vr;
	}
	break;
	default:
		return 0;//b*log(fabs(27-y)/0.1+1)
	}
}
double SpeedProfile::Vav(double x, double y, double z)
{
	switch (velocity_method){
	case 1://"Eulerian time-domain ...",Acta Acustica united with Acustica,Vol.88(2002) 483-492
	{
		//////only for H=4////
		if (z >= 2.8){
			return b*2.8;
		}
		else if (z >= 0.0){
			return b*z;
		}
		else{
			return 0.0;
		}
	}
	break;
	case 2://"study of the sound-vortex interaction ...", Eur.Phys.J.B 32,237-242(2003)
	{
		double Rr;
		Rr = sqrt(pow(y - 5.0, 2) + pow(z - 2.0, 2));
		if (abs(y - 5.0) <= 1.0e-6 && abs(z - 2.0) <= 1.0e-6) alpha = 0.0;
		else if (y - 5.0 >= 0.0) alpha = asin((y - 5.0) / Rr);
		else if (y - 5.0 <= 0.0) alpha = asin((y - 5.0) / Rr) + 3.1415926;
		return -1.0*2.0 / (Rr + 2.0)*sin(alpha);
		/*
					double r,Vr;
					r=sqrt(pow((y-y0),2)+pow((z-z0),2));
					if((r/L)<=0.000001) r=1e-20;
					Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
					return (z-z0)/r*Vr;
					*/
	}
	break;
	case 3://"transmission of sound through a single vortex", Eur.Phys.J.B 37,229-239(2004)
	{
		double r, Vr;
		r = sqrt(pow((y - y0), 2) + pow((z - z0), 2));
		if (r <= L){//Vr=circulation/(2*PI*L)*r/L;Vav=(z-z0)/r*Vr
			return circulation / (2 * PI*L) / L*(z - z0);
		}
		Vr = circulation / (2 * PI) / r;
		return (z - z0) / r*Vr;
	}
	break;
	case 4://reverse coordinate from case 1
	{
		return 0.0;
	}
	break;
	case 5://book of Solomans: terrain
	{
		return b*log(fabs(z) / 0.1 + 1);
	}
	break;
	case 6:
	{
		double r, Vr;
		r = sqrt(pow((y - y0), 2) + pow((z - z0), 2));
		if ((r / L) <= 0.00001) return 0;
		Vr = circulation / (2 * PI*r)*(1 - exp(-alpha*r*r / L / L));
		return (z - z0) / r*Vr;
	}
	break;
	default:
		return 0;
	}
}

double SpeedProfile::Wav(double x, double y, double z)
{
	switch (velocity_method){
	case 1:
	{
		return 0;
	}
	break;
	case 2:
	{
		double Rr;
		Rr = sqrt(pow(y - 5.0, 2) + pow(z - 2.0, 2));
		if (abs(y - 5.0) <= 1.0e-6 && abs(z - 2.0) <= 1.0e-6) alpha = 0.0;
		else if (y - 5.0 >= 0.0) alpha = asin((y - 5.0) / Rr);
		else if (y - 5.0 <= 0.0) alpha = asin((y - 5.0) / Rr) + 3.1415926;
		return 1.0*2.0 / (Rr + 2.0)*cos(alpha);
		/*
					double r,Vr;
					r=sqrt(pow((y-y0),2)+pow((z-z0),2));
					if((r/L)<=0.000001) r=1e-10;
					Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
					return -(y-y0)/r*Vr;
					*/
	}

	break;
	case 3:
	{
		double r, Vr;
		r = sqrt(pow((y - y0), 2) + pow((z - z0), 2));
		if (r <= L){//Vr=circulation/(2*PI*L)*r/L;Wav=-(y-y0)/r*Vr;
			return circulation / (2 * PI*L) / L*(-(y - y0));
		}
		Vr = circulation / (2 * PI) / r;
		return -(y - y0) / r*Vr;
	}
	break;
	case 4:
	{
		return 0;
	}
	break;
	case 5:
	{
		return b*y;
	}
	break;
	case 6:
	{
		double r, Vr;
		r = sqrt(pow((y - y0), 2) + pow((z - z0), 2));
		if ((r / L) <= 0.00001) return 0;
		Vr = circulation / (2 * PI*r)*(1 - exp(-alpha*r*r / L / L));
		return -(y - y0) / r*Vr;
	}
	break;
	default:
		return 0;//b*log(fabs(27-y)/0.1+1)
	}
}
SpeedProfile::~SpeedProfile()
{
}
//---------------------------------air member function---------------------------------------//

air::air(scheme_list scheme1, DifferenceStep Step, char *coordi, MovingFrame MF, const double gauss_width,
	double velocity_coef, int velocity_method1, AirStruct AirPara, PML_boundary PML_xbd,
	PML_boundary PML_ybd, PML_boundary PML_zbd, hillpore hill, hillpore canopy, hillpore trunk, MPI_info mpi_var, int judge_porous, int restart_timer)
	:calculation_2D(scheme1, Step, coordi, MF, gauss_width, mpi_var, judge_porous)
{
	int i, j, k;
	adiabatic_coef = AirPara.adiabatic_coef;
	Pav = AirPara.Pav; Aav = AirPara.Aav; sound_speed = AirPara.sound_speed;
	velocity_method = velocity_method1;
	SP = new SpeedProfile(velocity_coef, AirPara.sound_speed, velocity_method);

	//Allocate the meomerty to velocity and Q velocity for PML;
	Uav = new double**[mpi_IMAX]; Vav = new double**[mpi_IMAX]; Wav = new double**[mpi_IMAX];
	for (i = 0; i < mpi_IMAX; i++){
		Uav[i] = new double*[mpi_JMAX]; Vav[i] = new double*[mpi_JMAX]; Wav[i] = new double*[mpi_JMAX];
		for (j = 0; j < mpi_JMAX; j++){
			Uav[i][j] = new double[mpi_KMAX]; Vav[i][j] = new double[mpi_KMAX]; Wav[i][j] = new double[mpi_KMAX];
		}
	}
	//Initialize the average velocity and Q velocity for PML;
	for (i = 0; i < mpi_IMAX; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				Uav[i][j][k] = 0.0; Vav[i][j][k] = 0.0; Wav[i][j][k] = 0.0;
			}
		}
	}

	// Spcify the variables of PML boundary at the different location;
	cal_PML_var(PML_xbd, PML_ybd, PML_zbd);
	// specify the porous media in the hill;

	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0) cout << "PMLs all set" << endl;

	cal_hill_var(hill,canopy,trunk);

	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0) cout << "Barriers all set" << endl;

}

air::~air()
{
	int i, j;
	for (i = 0; i < mpi_IMAX; i++){
		for (j = 0; j < mpi_JMAX; j++){
			delete[] Uav[i][j]; delete[] Vav[i][j]; delete[] Wav[i][j];
		}
		delete[] Uav[i]; delete[] Vav[i]; delete[] Wav[i];
	}
	delete[] Uav; delete[] Vav; delete[] Wav;
	// Array: hill or barriers
	if (cr_judge != 0 && (istop - istar + 1) != 0 &&
		(jstop - jstar + 1) != 0 && (kstop - kstar + 1) != 0){
		for (i = 0; i < istop - istar + 1; i++){
			for (j = 0; j < jstop - jstar + 1; j++){
				delete[] cr_uc[i][j];
			}
			delete[] cr_uc[i];
		}
		delete[] cr_uc;
	}

	// Array:PML boudnary condition for each direction
	if (PML_judgex != 0){
		switch (PML_judgex){
		case 1:
			if (mpi_IMAX1 != 0){
				for (i = 0; i < mpi_IMAX1; i++){
					for (j = 0; j < mpi_JMAX; j++){
						delete[] Qux1[i][j]; delete[] Qvx1[i][j];
						delete[] Qwx1[i][j]; delete[] Qpx1[i][j];
					}
					delete[] Qux1[i]; delete[] Qvx1[i];
					delete[] Qwx1[i]; delete[] Qpx1[i];
				}
				delete[] Qux1; delete[] Qvx1;
				delete[] Qwx1; delete[] Qpx1;
			}
			break;
		case 2:
			if (mpi_IMAX2 != 0){
				for (i = 0; i < mpi_IMAX2; i++){
					for (j = 0; j < mpi_JMAX; j++){
						delete[] Qux2[i][j]; delete[] Qvx2[i][j];
						delete[] Qwx2[i][j]; delete[] Qpx2[i][j];
					}
					delete[] Qux2[i]; delete[] Qvx2[i];
					delete[] Qwx2[i]; delete[] Qpx2[i];
				}
				delete[] Qux2; delete[] Qvx2;
				delete[] Qwx2; delete[] Qpx2;
			}
			break;
		case 3:
			if (mpi_IMAX1 != 0){
				for (i = 0; i < mpi_IMAX1; i++){
					for (j = 0; j < mpi_JMAX; j++){
						delete[] Qux1[i][j]; delete[] Qvx1[i][j];
						delete[] Qwx1[i][j]; delete[] Qpx1[i][j];
					}
					delete[] Qux1[i]; delete[] Qvx1[i];
					delete[] Qwx1[i]; delete[] Qpx1[i];
				}
				delete[] Qux1; delete[] Qvx1;
				delete[] Qwx1; delete[] Qpx1;
			}
			if (mpi_IMAX2 != 0){
				for (i = 0; i < mpi_IMAX2; i++){
					for (j = 0; j < mpi_JMAX; j++){
						delete[] Qux2[i][j]; delete[] Qvx2[i][j];
						delete[] Qwx2[i][j]; delete[] Qpx2[i][j];
					}
					delete[] Qux2[i]; delete[] Qvx2[i];
					delete[] Qwx2[i]; delete[] Qpx2[i];
				}
				delete[]Qux2; delete[] Qvx2;
				delete[] Qwx2; delete[] Qpx2;
			}
			break;
		default:
			cout << "please entrance your choice again!";
		}
	}

	if (PML_judgey != 0){
		switch (PML_judgey){
		case 1:
			if (mpi_JMAX1 != 0){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX1; j++){
						delete[] Quy1[i][j]; delete[] Qvy1[i][j];
						delete[] Qwy1[i][j]; delete[] Qpy1[i][j];
					}
					delete[] Quy1[i]; delete[] Qvy1[i];
					delete[] Qwy1[i]; delete[] Qpy1[i];
				}
				delete[] Quy1; delete[] Qvy1;
				delete[] Qwy1; delete[] Qpy1;
			}
			break;
		case 2:
			if (mpi_JMAX2 != 0){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX2; j++){
						delete[] Quy2[i][j]; delete[] Qvy2[i][j];
						delete[] Qwy2[i][j]; delete[] Qpy2[i][j];
					}
					delete[] Quy2[i]; delete[] Qvy2[i];
					delete[] Qwy2[i]; delete[] Qpy2[i];
				}
				delete[] Quy2; delete[] Qvy2;
				delete[] Qwy2; delete[] Qpy2;
			}
			break;
		case 3:
			if (mpi_JMAX1 != 0){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX1; j++){
						delete[] Quy1[i][j]; delete[] Qvy1[i][j];
						delete[] Qwy1[i][j]; delete[] Qpy1[i][j];
					}
					delete[] Quy1[i]; delete[] Qvy1[i];
					delete[] Qwy1[i]; delete[] Qpy1[i];
				}
				delete[] Quy1; delete[] Qvy1;
				delete[] Qwy1; delete[] Qpy1;
			}
			if (mpi_JMAX2 != 0){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX2; j++){
						delete[] Quy2[i][j]; delete[] Qvy2[i][j];
						delete[] Qwy2[i][j]; delete[] Qpy2[i][j];
					}
					delete[] Quy2[i]; delete[] Qvy2[i];
					delete[] Qwy2[i]; delete[] Qpy2[i];
				}
				delete[] Quy2; delete[] Qvy2;
				delete[] Qwy2; delete[] Qpy2;
			}
			break;
		default:
			cout << "please entrance your choice again!";
		}
	}

	if (PML_judgez != 0){
		switch (PML_judgez){
		case 1:
			if (mpi_KMAX1 != 0){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX; j++){
						delete[] Quz1[i][j]; delete[] Qvz1[i][j];
						delete[] Qwz1[i][j]; delete[] Qpz1[i][j];
					}
					delete[] Quz1[i]; delete[] Qvz1[i];
					delete[] Qwz1[i]; delete[] Qpz1[i];
				}
				delete[] Quz1; delete[] Qvz1;
				delete[] Qwz1; delete[] Qpz1;
			}
			break;
		case 2:
			if (mpi_KMAX2 != 0){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX; j++){
						delete[] Quz2[i][j]; delete[] Qvz2[i][j];
						delete[] Qwz2[i][j]; delete[] Qpz2[i][j];
					}
					delete[] Quz2[i]; delete[] Qvz2[i];
					delete[] Qwz2[i]; delete[] Qpz2[i];
				}
				delete[] Quz2; delete[] Qvz2;
				delete[] Qwz2; delete[] Qpz2;
			}
			break;
		case 3:
			if (mpi_KMAX1 != 0){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX; j++){
						delete[] Quz1[i][j]; delete[] Qvz1[i][j];
						delete[] Qwz1[i][j]; delete[] Qpz1[i][j];
					}
					delete[] Quz1[i]; delete[] Qvz1[i];
					delete[] Qwz1[i]; delete[] Qpz1[i];
				}
				delete[] Quz1; delete[] Qvz1;
				delete[] Qwz1; delete[] Qpz1;
			}
			if (mpi_KMAX2 != 0){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX; j++){
						delete[] Quz2[i][j]; delete[] Qvz2[i][j];
						delete[] Qwz2[i][j]; delete[] Qpz2[i][j];
					}
					delete[] Quz2[i]; delete[] Qvz2[i];
					delete[] Qwz2[i]; delete[] Qpz2[i];
				}
				delete[] Quz2; delete[] Qvz2;
				delete[] Qwz2; delete[] Qpz2;
			}
			break;
		default:
			cout << "please entrance your choice again!";
		}
	}
	delete SP;

}
void air::cal_PML_var(PML_boundary PML_xbd, PML_boundary PML_ybd, PML_boundary PML_zbd)
{
	int i, j;
	PML_judgex = PML_xbd.judge;
	PML_judgey = PML_ybd.judge;
	PML_judgez = PML_zbd.judge;
	mpi_IMAX1 = 0;
	mpi_IMAX2 = 0;
	mpi_JMAX1 = 0;
	mpi_JMAX2 = 0;
	mpi_KMAX1 = 0;
	mpi_KMAX2 = 0;
	PML_x1 = 0; PML_x2 = 0;
	PML_y1 = 0; PML_y2 = 0;
	PML_z1 = 0; PML_z2 = 0;
	PML_Width1 = 0;  PML_Width2 = 0;
	PML_Length1 = 0;  PML_Length2 = 0;
	PML_Height1 = 0; PML_Height2 = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	if (PML_judgex != 0){
		//if (mpi_rank == 0) cout << "MPI_judgex = " << PML_judgex << endl;
		PML_AbsorbXmax = PML_xbd.resistivity;
		PML_alphax = PML_xbd.alpha;
		switch (PML_judgex){
		case 1:
			PML_Width1 = 1.0 / fabs(axisbound_X[IMAX1 - 1] - axisbound_X[0]);
			PML_x1 = axisbound_X[IMAX1 - 1];
			IMAX2 = 0;

			//MPI_Barrier(MPI_COMM_WORLD);
			//if (mpi_rank == 0) cout << "PML_Width1 = " << PML_Width1 << ", PML_x1 = " << PML_x1 << ", IMAX1 = " << IMAX1 << ", mpi_IMAX = " << mpi_IMAX << endl;

			if (mpi_i1 >= 0 && mpi_i1 < IMAX1 + 1){
				if (mpi_i2 >= IMAX1) mpi_IMAX1 = IMAX1 - mpi_i1 + 1;
				else mpi_IMAX1 = mpi_IMAX;
			}
			else mpi_IMAX1 = 0;

			//MPI_Barrier(MPI_COMM_WORLD);
			//if (mpi_rank == 0) cout << "mpi_i1 = " << mpi_i1 << ", mpi_i2 = " << mpi_i2 << ", mpi_IMAX1 = " << mpi_IMAX1 << ", mpi_JMAX = " << mpi_JMAX << ", mpi_KMAX = " << mpi_KMAX << endl;

			if (mpi_IMAX1 != 0){
				Qux1 = new double**[mpi_IMAX1];
				Qvx1 = new double**[mpi_IMAX1];
				Qwx1 = new double**[mpi_IMAX1];
				Qpx1 = new double**[mpi_IMAX1];
				for (i = 0; i < mpi_IMAX1; i++){
					Qux1[i] = new double*[mpi_JMAX];
					Qvx1[i] = new double*[mpi_JMAX];
					Qwx1[i] = new double*[mpi_JMAX];
					Qpx1[i] = new double*[mpi_JMAX];
					for (j = 0; j < mpi_JMAX; j++){
						Qux1[i][j] = new double[mpi_KMAX];
						Qvx1[i][j] = new double[mpi_KMAX];
						Qwx1[i][j] = new double[mpi_KMAX];
						Qpx1[i][j] = new double[mpi_KMAX];
					}

				}
			}
			//MPI_Barrier(MPI_COMM_WORLD);
			//cout << "mpi rank = " << mpi_rank << ", memory allocated! " << endl;
			break;
		case 2:
			PML_Width2 = 1.0 / fabs(axisbound_X[IMAX - 1] - axisbound_X[IMAX - IMAX2]);
			PML_x2 = axisbound_X[IMAX - IMAX2];
			IMAX1 = 0;
			if (mpi_i2 >= IMAX - IMAX2 - 1 && mpi_i2 < IMAX){
				if (mpi_i1 <= IMAX - IMAX2 - 1) mpi_IMAX2 = mpi_i2 - (IMAX - IMAX2 - 1) + 1;
				else mpi_IMAX2 = mpi_IMAX;
			}
			else mpi_IMAX2 = 0;
			if (mpi_IMAX2 != 0){
				Qux2 = new double **[mpi_IMAX2];
				Qvx2 = new double **[mpi_IMAX2];
				Qwx2 = new double **[mpi_IMAX2];
				Qpx2 = new double **[mpi_IMAX2];
				for (i = 0; i < mpi_IMAX2; i++){
					Qux2[i] = new double*[mpi_JMAX];
					Qvx2[i] = new double*[mpi_JMAX];
					Qwx2[i] = new double*[mpi_JMAX];
					Qpx2[i] = new double*[mpi_JMAX];
					for (j = 0; j < mpi_JMAX; j++){
						Qux2[i][j] = new double[mpi_KMAX];
						Qvx2[i][j] = new double[mpi_KMAX];
						Qwx2[i][j] = new double[mpi_KMAX];
						Qpx2[i][j] = new double[mpi_KMAX];
					}

				}
			}
			break;
		case 3:
			PML_Width1 = 1.0 / fabs(axisbound_X[IMAX1 - 1] - axisbound_X[0]);
			PML_x1 = axisbound_X[IMAX1 - 1];
			PML_Width2 = 1.0 / fabs(axisbound_X[IMAX - 1] - axisbound_X[IMAX - IMAX2]);
			PML_x2 = axisbound_X[IMAX - IMAX2];
			if (mpi_i1 >= 0 && mpi_i1 < IMAX1 + 1){
				if (mpi_i2 >= IMAX1) mpi_IMAX1 = IMAX1 - mpi_i1 + 1;
				else mpi_IMAX1 = mpi_IMAX;
			}
			else mpi_IMAX1 = 0;

			if (mpi_i2 >= IMAX - IMAX2 - 1 && mpi_i2 < IMAX){
				if (mpi_i1 <= IMAX - IMAX2 - 1) mpi_IMAX2 = mpi_i2 - (IMAX - IMAX2 - 1) + 1;
				else mpi_IMAX2 = mpi_IMAX;
			}
			else mpi_IMAX2 = 0;

			if (mpi_IMAX1 != 0){
				Qux1 = new double**[mpi_IMAX1];
				Qvx1 = new double**[mpi_IMAX1];
				Qwx1 = new double**[mpi_IMAX1];
				Qpx1 = new double**[mpi_IMAX1];
				for (i = 0; i < mpi_IMAX1; i++){
					Qux1[i] = new double*[mpi_JMAX];
					Qvx1[i] = new double*[mpi_JMAX];
					Qwx1[i] = new double*[mpi_JMAX];
					Qpx1[i] = new double*[mpi_JMAX];
					for (j = 0; j < mpi_JMAX; j++){
						Qux1[i][j] = new double[mpi_KMAX];
						Qvx1[i][j] = new double[mpi_KMAX];
						Qwx1[i][j] = new double[mpi_KMAX];
						Qpx1[i][j] = new double[mpi_KMAX];
					}

				}
			}
			if (mpi_IMAX2 != 0){
				Qux2 = new double **[mpi_IMAX2];
				Qvx2 = new double **[mpi_IMAX2];
				Qwx2 = new double **[mpi_IMAX2];
				Qpx2 = new double **[mpi_IMAX2];
				for (i = 0; i < mpi_IMAX2; i++){
					Qux2[i] = new double*[mpi_JMAX];
					Qvx2[i] = new double*[mpi_JMAX];
					Qwx2[i] = new double*[mpi_JMAX];
					Qpx2[i] = new double*[mpi_JMAX];
					for (j = 0; j < mpi_JMAX; j++){
						Qux2[i][j] = new double[mpi_KMAX];
						Qvx2[i][j] = new double[mpi_KMAX];
						Qwx2[i][j] = new double[mpi_KMAX];
						Qpx2[i][j] = new double[mpi_KMAX];
					}

				}
			}
			break;
		default:
			cout << "please entrance your choice again!";
		}
	}
	else{
		IMAX1 = 0;
		IMAX2 = 0;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (PML_judgey != 0){
		//if (mpi_rank == 0) cout << "MPI_judgex = " << PML_judgex << endl;
		PML_AbsorbYmax = PML_ybd.resistivity;
		PML_alphay = PML_ybd.alpha;
		switch (PML_judgey){
		case 1:
			PML_Length1 = 1.0 / fabs(axisbound_Y[JMAX1 - 1] - axisbound_Y[0]);
			PML_y1 = axisbound_Y[JMAX1 - 1];
			PML_JMAX1_time = PML_ybd.time_step;
			PML_JMAX1 = 0;
			JMAX2 = 0;
			if (mpi_j1 >= 0 && mpi_j1 < JMAX1 + 1){
				if (mpi_j2 >= JMAX1) mpi_JMAX1 = JMAX1 - mpi_j1 + 1;
				else mpi_JMAX1 = mpi_JMAX;
			}
			else mpi_JMAX1 = 0;
			if (mpi_JMAX1 != 0){
				Quy1 = new double**[mpi_IMAX];
				Qvy1 = new double**[mpi_IMAX];
				Qwy1 = new double**[mpi_IMAX];
				Qpy1 = new double**[mpi_IMAX];
				for (i = 0; i < mpi_IMAX; i++){
					Quy1[i] = new double*[mpi_JMAX1];
					Qvy1[i] = new double*[mpi_JMAX1];
					Qwy1[i] = new double*[mpi_JMAX1];
					Qpy1[i] = new double*[mpi_JMAX1];
					for (j = 0; j < mpi_JMAX1; j++){
						Quy1[i][j] = new double[mpi_KMAX];
						Qvy1[i][j] = new double[mpi_KMAX];
						Qwy1[i][j] = new double[mpi_KMAX];
						Qpy1[i][j] = new double[mpi_KMAX];
					}

				}
			}
			break;
		case 2:
			PML_Length2 = 1.0 / fabs(axisbound_Y[JMAX - 1] - axisbound_Y[JMAX - JMAX2]);
			PML_y2 = axisbound_Y[JMAX - JMAX2];
			JMAX1 = 0;
			if (mpi_j2 >= JMAX - JMAX2 - 1 && mpi_j2 < JMAX){
				if (mpi_j1 <= JMAX - JMAX2 - 1) mpi_JMAX2 = mpi_j2 - (JMAX - JMAX2 - 1) + 1;
				else mpi_JMAX2 = mpi_JMAX;
			}
			else mpi_JMAX2 = 0;

			if (mpi_JMAX2 != 0){
				Quy2 = new double**[mpi_IMAX];
				Qvy2 = new double**[mpi_IMAX];
				Qwy2 = new double**[mpi_IMAX];
				Qpy2 = new double**[mpi_IMAX];
				for (i = 0; i < mpi_IMAX; i++){
					Quy2[i] = new double*[mpi_JMAX2];
					Qvy2[i] = new double*[mpi_JMAX2];
					Qwy2[i] = new double*[mpi_JMAX2];
					Qpy2[i] = new double*[mpi_JMAX2];
					for (j = 0; j < mpi_JMAX2; j++){
						Quy2[i][j] = new double[mpi_KMAX];
						Qvy2[i][j] = new double[mpi_KMAX];
						Qwy2[i][j] = new double[mpi_KMAX];
						Qpy2[i][j] = new double[mpi_KMAX];
					}

				}
			}
			break;
		case 3:
			PML_Length1 = 1.0 / fabs(axisbound_Y[JMAX1 - 1] - axisbound_Y[0]);
			PML_y1 = axisbound_Y[JMAX1 - 1];
			PML_JMAX1_time = PML_ybd.time_step;
			PML_JMAX1 = 0;
			PML_Length2 = 1.0 / fabs(axisbound_Y[JMAX - 1] - axisbound_Y[JMAX - JMAX2]);
			PML_y2 = axisbound_Y[JMAX - JMAX2];

			if (mpi_j1 >= 0 && mpi_j1 < JMAX1 + 1){
				if (mpi_j2 >= JMAX1) mpi_JMAX1 = JMAX1 - mpi_j1 + 1;
				else mpi_JMAX1 = mpi_JMAX;
			}
			else mpi_JMAX1 = 0;

			if (mpi_j2 >= JMAX - JMAX2 - 1 && mpi_j2 < JMAX){
				if (mpi_j1 <= JMAX - JMAX2 - 1)mpi_JMAX2 = mpi_j2 - (JMAX - JMAX2 - 1) + 1;
				else mpi_JMAX2 = mpi_JMAX;
			}
			else mpi_JMAX2 = 0;

			if (mpi_JMAX1 != 0){
				Quy1 = new double**[mpi_IMAX];
				Qvy1 = new double**[mpi_IMAX];
				Qwy1 = new double**[mpi_IMAX];
				Qpy1 = new double**[mpi_IMAX];
				for (i = 0; i < mpi_IMAX; i++){
					Quy1[i] = new double*[mpi_JMAX1];
					Qvy1[i] = new double*[mpi_JMAX1];
					Qwy1[i] = new double*[mpi_JMAX1];
					Qpy1[i] = new double*[mpi_JMAX1];
					for (j = 0; j < mpi_JMAX1; j++){
						Quy1[i][j] = new double[mpi_KMAX];
						Qvy1[i][j] = new double[mpi_KMAX];
						Qwy1[i][j] = new double[mpi_KMAX];
						Qpy1[i][j] = new double[mpi_KMAX];
					}

				}
			}
			if (mpi_JMAX2 != 0){
				Quy2 = new double**[mpi_IMAX];
				Qvy2 = new double**[mpi_IMAX];
				Qwy2 = new double**[mpi_IMAX];
				Qpy2 = new double**[mpi_IMAX];
				for (i = 0; i < mpi_IMAX; i++){
					Quy2[i] = new double*[mpi_JMAX2];
					Qvy2[i] = new double*[mpi_JMAX2];
					Qwy2[i] = new double*[mpi_JMAX2];
					Qpy2[i] = new double*[mpi_JMAX2];
					for (j = 0; j < mpi_JMAX2; j++){
						Quy2[i][j] = new double[mpi_KMAX];
						Qvy2[i][j] = new double[mpi_KMAX];
						Qwy2[i][j] = new double[mpi_KMAX];
						Qpy2[i][j] = new double[mpi_KMAX];
					}

				}
			}
			break;
		default:
			cout << "please entrance your choice again!";
		}
	}
	else{
		JMAX1 = 0;
		JMAX2 = 0;
	}

	if (PML_judgez != 0){
		PML_AbsorbZmax = PML_zbd.resistivity;
		PML_alphaz = PML_zbd.alpha;
		switch (PML_judgez){
		case 1:
			PML_Height1 = 1.0 / fabs(axisbound_Z[KMAX1 - 1] - axisbound_Z[0]);
			PML_z1 = axisbound_Z[KMAX1 - 1];
			KMAX2 = 0;

			if (mpi_k1 >= 0 && mpi_k1 < KMAX1 + 1){
				if (mpi_k2 >= KMAX1) mpi_KMAX1 = KMAX1 - mpi_k1 + 1;
				else mpi_KMAX1 = mpi_KMAX;
			}
			else mpi_KMAX1 = 0;
			if (mpi_KMAX1 != 0){
				Quz1 = new double**[mpi_IMAX];
				Qvz1 = new double**[mpi_IMAX];
				Qwz1 = new double**[mpi_IMAX];
				Qpz1 = new double**[mpi_IMAX];
				for (i = 0; i < mpi_IMAX; i++){
					Quz1[i] = new double*[mpi_JMAX];
					Qvz1[i] = new double*[mpi_JMAX];
					Qwz1[i] = new double*[mpi_JMAX];
					Qpz1[i] = new double*[mpi_JMAX];
					for (j = 0; j < mpi_JMAX; j++){
						Quz1[i][j] = new double[mpi_KMAX1];
						Qvz1[i][j] = new double[mpi_KMAX1];
						Qwz1[i][j] = new double[mpi_KMAX1];
						Qpz1[i][j] = new double[mpi_KMAX1];
					}

				}
			}
			break;
		case 2:
			PML_Height2 = 1.0 / fabs(axisbound_Z[KMAX - 1] - axisbound_Z[KMAX - KMAX2]);
			PML_z2 = axisbound_Z[KMAX - KMAX2];
			KMAX1 = 0;

			if (mpi_k2 >= KMAX - KMAX2 - 1 && mpi_k2 < KMAX){
				if (mpi_k1 <= KMAX - KMAX2 - 1)mpi_KMAX2 = mpi_k2 - (KMAX - KMAX2 - 1) + 1;
				else mpi_KMAX2 = mpi_KMAX;
			}
			else mpi_KMAX2 = 0;

			if (mpi_KMAX2 != 0){
				Quz2 = new double**[mpi_IMAX];
				Qvz2 = new double**[mpi_IMAX];
				Qwz2 = new double**[mpi_IMAX];
				Qpz2 = new double**[mpi_IMAX];
				for (i = 0; i < mpi_IMAX; i++){
					Quz2[i] = new double *[mpi_JMAX];
					Qvz2[i] = new double *[mpi_JMAX];
					Qwz2[i] = new double *[mpi_JMAX];
					Qpz2[i] = new double *[mpi_JMAX];
					for (j = 0; j < mpi_JMAX; j++){
						Quz2[i][j] = new double[mpi_KMAX2];
						Qvz2[i][j] = new double[mpi_KMAX2];
						Qwz2[i][j] = new double[mpi_KMAX2];
						Qpz2[i][j] = new double[mpi_KMAX2];
					}
				}
			}
			break;
		case 3:
			PML_Height1 = 1.0 / fabs(axisbound_Z[KMAX1 - 1] - axisbound_Z[0]);
			PML_z1 = axisbound_Z[KMAX1 - 1];
			PML_Height2 = 1.0 / fabs(axisbound_Z[KMAX - 1] - axisbound_Z[KMAX - KMAX2]);
			PML_z2 = axisbound_Z[KMAX - KMAX2];

			if (mpi_k1 >= 0 && mpi_k1 < KMAX1 + 1){
				if (mpi_k2 >= KMAX1) mpi_KMAX1 = KMAX1 - mpi_k1 + 1;
				else mpi_KMAX1 = mpi_KMAX;
			}
			else mpi_KMAX1 = 0;

			if (mpi_k2 >= KMAX - KMAX2 - 1 && mpi_k2 < KMAX){
				if (mpi_k1 <= KMAX - KMAX2 - 1)mpi_KMAX2 = mpi_k2 - (KMAX - KMAX2 - 1) + 1;
				else mpi_KMAX2 = mpi_KMAX;
			}
			else mpi_KMAX2 = 0;

			if (mpi_KMAX1 != 0){
				Quz1 = new double**[mpi_IMAX];
				Qvz1 = new double**[mpi_IMAX];
				Qwz1 = new double**[mpi_IMAX];
				Qpz1 = new double**[mpi_IMAX];
				for (i = 0; i < mpi_IMAX; i++){
					Quz1[i] = new double*[mpi_JMAX];
					Qvz1[i] = new double*[mpi_JMAX];
					Qwz1[i] = new double*[mpi_JMAX];
					Qpz1[i] = new double*[mpi_JMAX];
					for (j = 0; j < mpi_JMAX; j++){
						Quz1[i][j] = new double[mpi_KMAX1];
						Qvz1[i][j] = new double[mpi_KMAX1];
						Qwz1[i][j] = new double[mpi_KMAX1];
						Qpz1[i][j] = new double[mpi_KMAX1];
					}

				}
			}

			if (mpi_KMAX2 != 0){
				Quz2 = new double**[mpi_IMAX];
				Qvz2 = new double**[mpi_IMAX];
				Qwz2 = new double**[mpi_IMAX];
				Qpz2 = new double**[mpi_IMAX];
				for (i = 0; i < mpi_IMAX; i++){
					Quz2[i] = new double *[mpi_JMAX];
					Qvz2[i] = new double *[mpi_JMAX];
					Qwz2[i] = new double *[mpi_JMAX];
					Qpz2[i] = new double *[mpi_JMAX];
					for (j = 0; j < mpi_JMAX; j++){
						Quz2[i][j] = new double[mpi_KMAX2];
						Qvz2[i][j] = new double[mpi_KMAX2];
						Qwz2[i][j] = new double[mpi_KMAX2];
						Qpz2[i][j] = new double[mpi_KMAX2];
					}
				}
			}
			break;
		default:
			cout << "please entrance your choice again!";
		}
	}
	else{
		KMAX1 = 0;
		KMAX2 = 0;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0) {
		cout << "PML_x1 = " << PML_x1 << ", PML_x2 = " << PML_x2 << ", PML_Width1 = " << PML_Width1 << ", PML_Width2 = " << PML_Width2 << endl;
		cout << "PML_y1 = " << PML_y1 << ", PML_y2 = " << PML_y2 << ", PML_Length1 = " << PML_Length1 << ", PML_Length2 = " << PML_Length2 << endl;
		cout << "PML_z1 = " << PML_z1 << ", PML_z2 = " << PML_z2 << ", PML_Height1 = " << PML_Height1 << ", PML_Height2 = " << PML_Height2 << endl;
	}

}
void air::cal_hill_var(hillpore hill, hillpore canopy, hillpore trunk)
{
	int i, j, k,ii, jj;
	cr_judge = hill.curve_judge;
	if (cr_judge == 0) return;
	cr_a = hill.curve_coefa;
	cr_b = hill.curve_coefb;
	cr_c = hill.curve_coefc;
	cr_d = hill.curve_coefd;
	cr_e = hill.curve_coefe;
	cr_f = hill.curve_coeff;
	cr_g = hill.curve_coefg;
	cr_h = hill.curve_coefh;
	cr_i = hill.curve_coefi;
	cr_j = hill.curve_coefj;
	cr_k = hill.curve_coefk;
	cr_l = hill.curve_coefl;
	cr_points = hill.curve_points;
	double pi, cr_height, ang, cr_length;
	double k_slope1, k_slope2, k_slope3, k_slope4, k_slope5, k_slope6;
	crpo_eff_density = 1.0 / hill.eff_density;
	crpo_porosity = hill.porosity;
	crpo_resistivity = hill.resistivity;
	crpo_Kp = hill.Kp;
	crpo_Beta = crpo_resistivity*diff_t*crpo_eff_density;
	crpo_Gama = diff_t / crpo_Kp / crpo_porosity;
	crpo_invBeta = 1.0 / (1.0 + crpo_Beta);

	crpo_eff_density1 = 1.0 / canopy.eff_density;
	crpo_porosity1 = canopy.porosity;
	crpo_resistivity1 = canopy.resistivity;
	crpo_Kp1 = canopy.Kp;
	crpo_Beta1 = crpo_resistivity1*diff_t*crpo_eff_density1;
	crpo_Gama1 = diff_t / crpo_Kp1 / crpo_porosity1;
	crpo_invBeta1 = 1.0 / (1.0 + crpo_Beta1);

	crpo_eff_density2 = 1.0 / trunk.eff_density;
	crpo_porosity2 = trunk.porosity;
	crpo_resistivity2 = trunk.resistivity;
	crpo_Kp2 = trunk.Kp;
	crpo_Beta2 = crpo_resistivity2*diff_t*crpo_eff_density2;
	crpo_Gama2 = diff_t / crpo_Kp2 / crpo_porosity2;
	crpo_invBeta2 = 1.0 / (1.0 + crpo_Beta2);

	istar = 0;	istop = 0;	jstar = 0;
	jstop = 0;	kstar = 0;	kstop = 0;

	for (i = 0; i < IMAX; i++){
		if (fabs(axisbound_X[i] - hill.curve_xstar) <= diff_x1 / 8.0) istar = i;
		if (fabs(axisbound_X[i] - hill.curve_xstop) <= diff_x1 / 8.0){
			istop = i;
			break;
		}
	}
	for (j = 0; j < JMAX; j++){
		if (fabs(axisbound_Y[j] - hill.curve_ystar) <= diff_y1 / 8.0) jstar = j;
		if (fabs(axisbound_Y[j] - hill.curve_ystop) <= diff_y1 / 8.0){
			jstop = j;
			break;
		}
	}
	for (k = 0; k < KMAX; k++){
		if (fabs(axisbound_Z[k] - hill.curve_zstar) <= diff_z1 / 8.0) kstar = k;
		if (fabs(axisbound_Z[k] - hill.curve_zstop) <= diff_z1 / 8.0){
			kstop = k;
			break;
		}
	}

	if (mpi_rank == 0) cout << "istar = " << istar << ", istop = " << istop << ", jstar = " << jstar << ", jstop = " << jstop << ", kstar = " << kstar << ", kstop = " << kstop << endl;

	if ((istop - istar + 1) != 0 && (jstop - jstar + 1) != 0 && (kstop - kstar + 1) != 0){
		cr_uc = new int**[istop - istar + 1];
		for (i = 0; i < istop - istar + 1; i++){
			cr_uc[i] = new int*[jstop - jstar + 1];
			for (j = 0; j < jstop - jstar + 1; j++){
				cr_uc[i][j] = new int[kstop - kstar + 1];
			}
		}
	}

	switch (cr_judge){
	case 1: 
		if ((istop - istar + 1) != 0 && (jstop - jstar + 1) != 0 && (kstop - kstar + 1) != 0){
			for (i = istar; i < istop + 1; i++){
				for (j = jstar; j < jstop + 1; j++){
					for (k = kstar; k < kstop + 1; k++){
						cr_uc[i - istar][j - jstar][k - kstar] = 0;
						cr_length = pow((axisbound_X[i] - cr_c), 2) + pow((axisbound_Y[j] - cr_d), 2)
							+ pow((axisbound_Z[k] - cr_e), 2);
						if (sqrt(cr_length) <= cr_b) cr_uc[i - istar][j - jstar][k - kstar] = 1;
					}
				}
			}
		}
		break;
	case 2:
		if ((istop - istar + 1) != 0 && (jstop - jstar + 1) != 0 && (kstop - kstar + 1) != 0){
			for (i = istar; i < istop + 1; i++){
				for (j = jstar; j < jstop + 1; j++){
					for (k = kstar; k < kstop + 1; k++){
						cr_uc[i - istar][j - jstar][k - kstar] = 0;
						
						if (axisbound_X[i] >= cr_a && axisbound_X[i] <= cr_b && axisbound_Y[j] >= cr_c && axisbound_Y[j] <= cr_d && axisbound_Z[k] >= cr_e && axisbound_Z[k] <= cr_f)
							cr_uc[i - istar][j - jstar][k - kstar] = 1;
						//cr_length = pow((axisbound_X[i] - cr_c), 2) + pow((axisbound_Y[j] - cr_d), 2);
						//if (sqrt(cr_length) <= cr_b && axisbound_Z[k] <= cr_e) cr_uc[i - istar][j - jstar][k - kstar] = 1;
					}
				}
			}
		}
		break;
	case 3:
		if ((istop - istar + 1) != 0 && (jstop - jstar + 1) != 0 && (kstop - kstar + 1) != 0){
			for (i = istar; i < istop + 1; i++){
				for (j = jstar; j < jstop + 1; j++){
					for (k = kstar; k < kstop + 1; k++){
						cr_uc[i - istar][j - jstar][k - kstar] = 0;
						for (jj = 0; jj < cr_points; jj++){
							cr_length = pow((axisbound_X[i] - cr_c), 2) + pow((axisbound_Y[j] - (cr_d + jj*cr_b)), 2);
							if (sqrt(cr_length) <= cr_a && axisbound_Z[k] <= cr_e)cr_uc[i - istar][j - jstar][k - kstar] = 1;

							cr_length = pow((axisbound_X[i] - cr_f), 2) + pow((axisbound_Y[j] - (cr_g + jj*cr_b)), 2);
							if (sqrt(cr_length) <= cr_a && axisbound_Z[k] <= cr_h)cr_uc[i - istar][j - jstar][k - kstar] = 1;

							cr_length = pow((axisbound_X[i] - cr_i), 2) + pow((axisbound_Y[j] - (cr_j + jj*cr_b)), 2);
							if (sqrt(cr_length) <= cr_a && axisbound_Z[k] <= cr_k)cr_uc[i - istar][j - jstar][k - kstar] = 1;
						}
					}
				}
			}
		}
		break;
	case 4: // cylindrical array
		if ((istop - istar + 1) != 0 && (jstop - jstar + 1) != 0 && (kstop - kstar + 1) != 0){
			for (i = istar; i < istop + 1; i++){
				for (j = jstar; j < jstop + 1; j++){
					for (k = kstar; k < kstop + 1; k++){
						cr_uc[i - istar][j - jstar][k - kstar] = 0;
						for (ii = 0; ii < cr_g; ii++){
							for (jj = 0; jj < cr_f; jj++){
								cr_length = pow((axisbound_X[i] - (cr_b + ii * cr_e)), 2) + pow((axisbound_Y[j] - (cr_c + jj*cr_d)), 2);
								if (sqrt(cr_length) <= cr_a ) 
									cr_uc[i - istar][j - jstar][k - kstar] = 1;
							}
						}
					}
				}
			}
		}
		break;
	case 5:
		if ((istop - istar + 1) != 0 && (jstop - jstar + 1) != 0 && (kstop - kstar + 1) != 0){
			for (i = istar; i < istop + 1; i++){
				for (j = jstar; j < jstop + 1; j++){
					for (k = kstar; k < kstop + 1; k++){
						cr_uc[i - istar][j - jstar][k - kstar] = 0;
						for (ii = 0; ii < 1; ii++){
							for (jj = 0; jj < 1; jj++){
								cr_length = pow((axisbound_X[i] - (0.0 + ii * 2.0)), 2) + pow((axisbound_Y[j] - (3.0 + jj*2.0)), 2);
								cr_height = -25 /(9-1e-50) * sqrt(cr_length) + 3.0;
								if (sqrt(cr_length) <= 0.9 && axisbound_Z[k] <= cr_height && axisbound_Z[k] >= 0.5) cr_uc[i - istar][j - jstar][k - kstar] = 1;

								if (sqrt(cr_length) <= 0.1 && axisbound_Z[k] <= 1.5) cr_uc[i - istar][j - jstar][k - kstar] = 2;
							}
						}
					}
				}
			}
		}
		if ((istop - istar + 1) != 0 && (jstop - jstar + 1) != 0 && (kstop - kstar + 1) != 0){
			for (i = istar; i < istop + 1; i++){
				for (j = jstar; j < jstop + 1; j++){
					for (k = kstar; k < kstop + 1; k++){
						cr_length = pow((axisbound_X[i] - (0.0 + ii * 2.0)), 2) + pow((axisbound_Y[j] - (3.0 + jj*2.0)), 2);

						if (sqrt(cr_length) <= 0.1 && axisbound_Z[k] <= 1.5) cr_uc[i - istar][j - jstar][k - kstar] = 2;
					}
				}
			}
		}	
		break;
	default:
		cout << "please entrance your choice again!";
	}

	/*
		cr_y=new double [cr_points]; cr_z=new double [cr_points]; cr_uc=new int** [IMAX];
		for (i=0;i<IMAX;i++){
		cr_uc[i]=new int* [JMAX];
		for (j=0;j<JMAX;j++){
		cr_uc[i][j]=new int [KMAX];
		}
		}

		for (i=0;i<IMAX;i++){
		for (j=0;j<JMAX;j++){
		for (k=0;k<KMAX;k++){
		cr_uc[i][j][k]=0;
		}
		}
		}
		if (cr_judge!=0){
		istar=561;istop=1441;
		jstar=161;jstop=1041;
		kstar=1;kstop=12;

		ofstream outfile_circle("curve_xy.dat",ios::out | ios::binary);
		outfile_circle.setf(ios::scientific,ios::floatfield);
		outfile_circle.precision(6);
		pi=3.1415926;
		cr_height=-999.0;
		if (cr_judge==1){
		ang=(pi*2)/cr_points;
		for (i=0;i<cr_points;i++){
		cr_y[i]=cr_c+cr_b*cos(i*ang);    // (cr_c,cr_d)=(xc,yc), R_radius=cr_b;
		cr_z[i]=cr_d+cr_b*sin(i*ang);
		outfile_circle<<cr_y[i]<<" "<<cr_z[i]<<endl;
		}
		for (j=jstar;j<=jstop;j++){
		for (i=istar;i<=istop;i++){
		cr_length=pow((whole_Y[i][j]-cr_c),2)+pow((whole_Z[i][j]-cr_d),2);
		if (sqrt(cr_length)<cr_b){
		cr_uc[i][j]=1;
		}
		}
		}
		}else if(cr_judge==2){
		for (i=0;i<cr_points;i++){
		cr_y[i]=(i-1)*6.0/(cr_points-1);
		cr_z[i]=cr_a*pow(cr_y[i],3)+cr_b*pow(cr_y[i],2)+
		cr_c*cr_y[i]+cr_d;
		cr_z[i]=sqrt(cr_z[i])-4.0;
		outfile_circle<<cr_y[i]<<" "<<cr_z[i]<<endl;
		}

		for (j=jstar;j<=jstop;j++){
		for (i=istar;i<=istop;i++){
		cr_height=cr_a*pow(whole_Y[i][j],3)
		+cr_b*pow(whole_Y[i][j],2)+
		cr_c*whole_Y[i][j]+cr_d;
		cr_height=sqrt(cr_height)-4.0;
		if (cr_height>whole_Z[i][j]){
		cr_uc[i][j]=1;
		}
		}

		}
		}else if(cr_judge==3){ // triangular
		k_slope1=cr_c/(cr_b-cr_a);
		k_slope2=cr_c/(cr_b-cr_d);
		for (j=jstar;j<=jstop;j++){
		for(i=istar;i<=istop;i++){
		if (whole_Y[i][j]>=cr_a && whole_Y[i][j]<=cr_b){
		cr_height=k_slope1*(whole_Y[i][j]-cr_a);
		if (cr_height>=whole_Z[i][j]){
		cr_uc[i][j]=1;
		}
		}else if(whole_Y[i][j]>=cr_b && whole_Y[i][j]<=cr_d){
		cr_height=k_slope2*(whole_Y[i][j]-cr_d);
		if (cr_height>=whole_Z[i][j]){
		cr_uc[i][j]=1;
		}
		}
		}
		}
		}else if (cr_judge==4){// rectangular
		cr_height=cr_c;
		for (j=jstar;j<=jstop;j++){
		for(i=istar;i<=istop;i++){
		if (whole_Y[i][j]>=cr_b && whole_Y[i][j]<=cr_d &&
		whole_Z[i][j]<=cr_height){
		cr_uc[i][j]=1;
		}
		}
		}
		}else if (cr_judge==5){// two wedges
		k_slope1=cr_c/(cr_b-cr_a);
		k_slope2=cr_c/(cr_b-cr_d);
		k_slope3=cr_g/(cr_f-cr_e);
		k_slope4=cr_g/(cr_f-cr_h);
		for (j=jstar;j<=jstop;j++){
		for(i=istar;i<=istop;i++){
		if (whole_Y[i][j]>=cr_a && whole_Y[i][j]<=cr_b){
		cr_height=k_slope1*(whole_Y[i][j]-cr_a);
		if (cr_height>=whole_Z[i][j]){
		cr_uc[i][j]=1;
		}
		}else if(whole_Y[i][j]>=cr_b && whole_Y[i][j]<=cr_d){
		cr_height=k_slope2*(whole_Y[i][j]-cr_d);
		if (cr_height>=whole_Z[i][j]){
		cr_uc[i][j]=1;
		}
		}else if(whole_Y[i][j]>=cr_e && whole_Y[i][j]<=cr_f){
		cr_height=k_slope3*(whole_Y[i][j]-cr_e);
		if (cr_height>=whole_Z[i][j]){
		cr_uc[i][j]=1;
		}
		}else if(whole_Y[i][j]>=cr_f && whole_Y[i][j]<=cr_h){
		cr_height=k_slope4*(whole_Y[i][j]-cr_h);
		if (cr_height>=whole_Z[i][j]){
		cr_uc[i][j]=1;
		}
		}
		}
		}
		}else if (cr_judge==6){ // two rectangulars
		//	        cr_height=cr_c;
		for (j=jstar;j<=jstop;j++){
		for(i=istar;i<=istop;i++){
		if (whole_Y[i][j]>=cr_b && whole_Y[i][j]<=cr_d &&
		whole_Z[i][j]<=cr_c){
		cr_uc[i][j]=1;
		}else if (whole_Y[i][j]>=cr_f && whole_Y[i][j]<=cr_h &&
		whole_Z[i][j]<=cr_g){
		cr_uc[i][j]=1;
		}
		}
		}
		}else if (cr_judge==7){
		k_slope1=cr_c/(cr_b-cr_a);
		k_slope2=cr_c/(cr_b-cr_d);
		for (j=jstar;j<=jstop;j++){
		for(i=istar;i<=istop;i++){
		if (whole_Y[i][j]>=cr_a && whole_Y[i][j]<=cr_b){
		cr_height=k_slope1*(whole_Y[i][j]-cr_a);
		if (cr_height>=whole_Z[i][j]){
		cr_uc[i][j]=1;
		}
		}else if(whole_Y[i][j]>=cr_b && whole_Y[i][j]<=cr_d){
		cr_height=k_slope2*(whole_Y[i][j]-cr_d);
		if (cr_height>=whole_Z[i][j]){
		cr_uc[i][j]=1;
		}
		}else if (whole_Y[i][j]>=cr_f && whole_Y[i][j]<=cr_h &&
		whole_Z[i][j]<=cr_g){
		cr_uc[i][j]=1;
		}
		}
		}

		}
		outfile_circle.close();
		}
		*/

}

void air::SetWindProfile(int N)
{
	int i,j,k;
	if (N!=0){
		PML_y1=axisbound_Y[JMAX1-1]+N*(move_frame.trail_DJ-1)*diff_y1;
		PML_y2=axisbound_Y[JMAX-JMAX2]+N*(move_frame.trail_DJ-1)*diff_y1;
	}
	for(i=0;i<mpi_IMAX;i++){
		for(j=0;j<mpi_JMAX;j++){
			for(k=0;k<mpi_KMAX;k++){
					double Z0;
					Z0=0.0;
					Uav[i][j][k]=SP->Uav(X[i],Y[j],Z[k]-Z0);
					Vav[i][j][k]=SP->Vav(X[i],Y[j],Z[k]-Z0);
					Wav[i][j][k]=SP->Wav(X[i],Y[j],Z[k]-Z0);
			}
		}
	}
}
void air::save_restart_air(char *restartfile)
{
	int i, j, k;
	ofstream outfile11(restartfile, ios::app | ios::binary);
	outfile11.setf(ios::scientific, ios::floatfield);
	outfile11.setf(ios::showpos);
	outfile11.precision(10);
	outfile11.width(18);

	for (i = 0; i < mpi_IMAX; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				outfile11 << Uav[i][j][k] << " ";
				outfile11 << Vav[i][j][k] << " ";
				outfile11 << Wav[i][j][k] << " ";
			}
		}
	}

	for (i = 0; i < mpi_IMAX1; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				outfile11 << Qux1[i][j][k] << " ";
				outfile11 << Qvx1[i][j][k] << " ";
				outfile11 << Qwx1[i][j][k] << " ";
				outfile11 << Qpx1[i][j][k] << " ";
			}
		}
	}
	for (i = 0; i < mpi_IMAX2; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				outfile11 << Qux2[i][j][k] << " ";
				outfile11 << Qvx2[i][j][k] << " ";
				outfile11 << Qwx2[i][j][k] << " ";
				outfile11 << Qpx2[i][j][k] << " ";
			}
		}
	}

	for (j = 0; j < mpi_JMAX1; j++){
		for (k = 0; k < mpi_KMAX; k++){
			for (i = 0; i < mpi_IMAX; i++){
				outfile11 << Quy1[i][j][k] << " ";
				outfile11 << Qvy1[i][j][k] << " ";
				outfile11 << Qwy1[i][j][k] << " ";
				outfile11 << Qpy1[i][j][k] << " ";
			}
		}
	}

	for (j = 0; j < mpi_JMAX2; j++){
		for (k = 0; k < mpi_KMAX; k++){
			for (i = 0; i < mpi_IMAX; i++){
				outfile11 << Quy2[i][j][k] << " ";
				outfile11 << Qvy2[i][j][k] << " ";
				outfile11 << Qwy2[i][j][k] << " ";
				outfile11 << Qpy2[i][j][k] << " ";
			}
		}
	}

	for (k = 0; k < mpi_KMAX1; k++){
		for (i = 0; i < mpi_IMAX; i++){
			for (j = 0; j < mpi_JMAX; j++){
				outfile11 << Quz1[i][j][k] << " ";
				outfile11 << Qvz1[i][j][k] << " ";
				outfile11 << Qwz1[i][j][k] << " ";
				outfile11 << Qpz1[i][j][k] << " ";
			}
		}
	}

	for (k = 0; k < mpi_KMAX2; k++){
		for (i = 0; i < mpi_IMAX; i++){
			for (j = 0; j < mpi_JMAX; j++){
				outfile11 << Quz2[i][j][k] << " ";
				outfile11 << Qvz2[i][j][k] << " ";
				outfile11 << Qwz2[i][j][k] << " ";
				outfile11 << Qpz2[i][j][k] << " ";
			}
		}
	}
	outfile11 << endl;
	outfile11.close();
}

void air::input_restart_air(char *restart_infile, int *point_pois)
{
	int i, j, k, int_pois;
	int_pois = *point_pois;
	ifstream infile(restart_infile, ios::in | ios::binary);
	infile.seekg(int_pois);
	for (i = 0; i < mpi_IMAX; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				infile >> Uav[i][j][k];
				infile >> Vav[i][j][k];
				infile >> Wav[i][j][k];
			}
		}
	}
	for (i = 0; i < mpi_IMAX1; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				infile >> Qux1[i][j][k];
				infile >> Qvx1[i][j][k];
				infile >> Qwx1[i][j][k];
				infile >> Qpx1[i][j][k];
			}
		}
	}
	for (i = 0; i < mpi_IMAX2; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				infile >> Qux2[i][j][k];
				infile >> Qvx2[i][j][k];
				infile >> Qwx2[i][j][k];
				infile >> Qpx2[i][j][k];
			}
		}
	}

	for (j = 0; j < mpi_JMAX1; j++){
		for (k = 0; k < mpi_KMAX; k++){
			for (i = 0; i < mpi_IMAX; i++){
				infile >> Quy1[i][j][k];
				infile >> Qvy1[i][j][k];
				infile >> Qwy1[i][j][k];
				infile >> Qpy1[i][j][k];
			}
		}
	}
	for (j = 0; j < mpi_JMAX2; j++){
		for (k = 0; k < mpi_KMAX; k++){
			for (i = 0; i < mpi_IMAX; i++){
				infile >> Quy2[i][j][k];
				infile >> Qvy2[i][j][k];
				infile >> Qwy2[i][j][k];
				infile >> Qpy2[i][j][k];
			}
		}
	}

	for (k = 0; k < mpi_KMAX1; k++){
		for (i = 0; i < mpi_IMAX; i++){
			for (j = 0; j < mpi_JMAX; j++){
				infile >> Quz1[i][j][k];
				infile >> Qvz1[i][j][k];
				infile >> Qwz1[i][j][k];
				infile >> Qpz1[i][j][k];
			}
		}
	}

	for (k = 0; k < mpi_KMAX2; k++){
		for (i = 0; i < mpi_IMAX; i++){
			for (j = 0; j < mpi_JMAX; j++){
				infile >> Quz2[i][j][k];
				infile >> Qvz2[i][j][k];
				infile >> Qwz2[i][j][k];
				infile >> Qpz2[i][j][k];
			}
		}
	}
	infile.ignore(100, '\n');
	int_pois = infile.tellg();
	*point_pois = int_pois;
	infile.close();
}
//--------------------------calculate velocity of V and W for new time (nn)----------------//
void air::cal_fuvw(double ***temp_u, double ***temp_v, double ***temp_w,
	double ***temp_p)
{
	int i, j, k, ii, jj, kk;
	int in, jn, kn;
	//  PML boundary for z direction
	double absorb_sigmaz1, absorb_sigmaz2, temp_z, absorb_betaz1, absorb_betaz2, absorb_coefzx, absorb_coefzy;
	//  PML boundary for x direction
	double absorb_sigmax1, absorb_sigmax2, temp_x, absorb_betax1, absorb_betax2, absorb_coefxz, absorb_coefxy;
	//  PML boundary for y direction
	double absorb_sigmay1, absorb_sigmay2, temp_y, absorb_betay1, absorb_betay2, absorb_coefyx, absorb_coefyz;

	if (scheme == FB_v_p && (PML_judgey == 1 || PML_judgey == 3)){
		PML_JMAX1 = PML_JMAX1 + 1;
		if (PML_JMAX1 >= PML_JMAX1_time)JMAX1 = JMAX3;
		else JMAX1 = 0;
	}
	for (i = 1; i < mpi_IMAX - 1; i++){
		for (j = 1; j < mpi_JMAX - 1; j++){
			for (k = 1; k < mpi_KMAX - 1; k++){
				in = mpi_i1 + i; jn = mpi_j1 + j; kn = mpi_k1 + k;
				x_over_X = Geo_over_var(Y_over_y[j], Z_over_z[k], Y_over_z,
					Z_over_y, Jacobian[i][j][k]);
				x_over_Y = Geo_over_var(X_over_z, Z_over_y, X_over_y,
					Z_over_z[k], Jacobian[i][j][k]);
				x_over_Z = Geo_over_var(X_over_y, Y_over_z, X_over_z,
					Y_over_y[j], Jacobian[i][j][k]);
				y_over_X = Geo_over_var(Y_over_z, Z_over_x, Y_over_x,
					Z_over_z[k], Jacobian[i][j][k]);
				y_over_Y = Geo_over_var(X_over_x[i], Z_over_z[k], X_over_z,
					Z_over_x, Jacobian[i][j][k]);
				y_over_Z = Geo_over_var(X_over_z, Y_over_x, X_over_x[i],
					Y_over_z, Jacobian[i][j][k]);
				z_over_X = Geo_over_var(Y_over_x, Z_over_y, Y_over_y[j],
					Z_over_x, Jacobian[i][j][k]);
				z_over_Y = Geo_over_var(X_over_y, Z_over_x, X_over_x[i],
					Z_over_y, Jacobian[i][j][k]);
				z_over_Z = Geo_over_var(X_over_x[i], Y_over_y[j], X_over_y,
					Y_over_x, Jacobian[i][j][k]);
				coef1 = (Uav[i][j][k] * x_over_X + Vav[i][j][k] * x_over_Y
					+ Wav[i][j][k] * x_over_Z)*diff_x*diff_t;
				coef2 = (Uav[i][j][k] * y_over_X + Vav[i][j][k] * y_over_Y
					+ Wav[i][j][k] * y_over_Z)*diff_y*diff_t;
				coef3 = (Uav[i][j][k] * z_over_X + Vav[i][j][k] * z_over_Y
					+ Wav[i][j][k] * z_over_Z)*diff_z*diff_t;
				Uav_over_X = x_over_X*(Uav[i + 1][j][k] - Uav[i - 1][j][k])*diff_x*0.5 +
					y_over_X*(Uav[i][j + 1][k] - Uav[i][j - 1][k])*diff_y*0.5 +
					z_over_X*(Uav[i][j][k + 1] - Uav[i][j][k - 1])*diff_z*0.5;
				Uav_over_Y = x_over_Y*(Uav[i + 1][j][k] - Uav[i - 1][j][k])*diff_x*0.5 +
					y_over_Y*(Uav[i][j + 1][k] - Uav[i][j - 1][k])*diff_y*0.5 +
					z_over_Y*(Uav[i][j][k + 1] - Uav[i][j][k - 1])*diff_z*0.5;
				Uav_over_Z = x_over_Z*(Uav[i + 1][j][k] - Uav[i - 1][j][k])*diff_x*0.5 +
					y_over_Z*(Uav[i][j + 1][k] - Uav[i][j - 1][k])*diff_y*0.5 +
					z_over_Z*(Uav[i][j][k + 1] - Uav[i][j][k - 1])*diff_z*0.5;
				Vav_over_X = x_over_X*(Vav[i + 1][j][k] - Vav[i - 1][j][k])*diff_x*0.5 +
					y_over_X*(Vav[i][j + 1][k] - Vav[i][j - 1][k])*diff_y*0.5 +
					z_over_X*(Vav[i][j][k + 1] - Vav[i][j][k - 1])*diff_z*0.5;
				Vav_over_Y = x_over_Y*(Vav[i + 1][j][k] - Vav[i - 1][j][k])*diff_x*0.5 +
					y_over_Y*(Vav[i][j + 1][k] - Vav[i][j - 1][k])*diff_y*0.5 +
					z_over_Y*(Vav[i][j][k + 1] - Vav[i][j][k - 1])*diff_z*0.5;
				Vav_over_Z = x_over_Z*(Vav[i + 1][j][k] - Vav[i - 1][j][k])*diff_x*0.5 +
					y_over_Z*(Vav[i][j + 1][k] - Vav[i][j - 1][k])*diff_y*0.5 +
					z_over_Z*(Vav[i][j][k + 1] - Vav[i][j][k - 1])*diff_z*0.5;
				Wav_over_X = x_over_X*(Wav[i + 1][j][k] - Wav[i - 1][j][k])*diff_x*0.5 +
					y_over_X*(Wav[i][j + 1][k] - Wav[i][j - 1][k])*diff_y*0.5 +
					z_over_X*(Wav[i][j][k + 1] - Wav[i][j][k - 1])*diff_z*0.5;
				Wav_over_Y = x_over_Y*(Wav[i + 1][j][k] - Wav[i - 1][j][k])*diff_x*0.5 +
					y_over_Y*(Wav[i][j + 1][k] - Wav[i][j - 1][k])*diff_y*0.5 +
					z_over_Y*(Wav[i][j][k + 1] - Wav[i][j][k - 1])*diff_z*0.5;
				Wav_over_Z = x_over_Z*(Wav[i + 1][j][k] - Wav[i - 1][j][k])*diff_x*0.5 +
					y_over_Z*(Wav[i][j + 1][k] - Wav[i][j - 1][k])*diff_y*0.5 +
					z_over_Z*(Wav[i][j][k + 1] - Wav[i][j][k - 1])*diff_z*0.5;

				u_nn[i][j][k] = -coef1*(temp_u[i + 1][j][k] - temp_u[i][j][k])
					- coef2*(temp_u[i][j + 1][k] - temp_u[i][j][k])
					- coef3*(temp_u[i][j][k + 1] - temp_u[i][j][k]);
				u_nn[i][j][k] += -temp_u[i][j][k] * Uav_over_X*diff_t
					- temp_v[i][j][k] * Uav_over_Y*diff_t
					- temp_w[i][j][k] * Uav_over_Z*diff_t;
				u_nn[i][j][k] += -Aav*diff_t*(x_over_X*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
					+ y_over_X*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
					+ z_over_X*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z);

				v_nn[i][j][k] = -coef1*(temp_v[i + 1][j][k] - temp_v[i][j][k])
					- coef2*(temp_v[i][j + 1][k] - temp_v[i][j][k])
					- coef3*(temp_v[i][j][k + 1] - temp_v[i][j][k]);
				v_nn[i][j][k] += -temp_u[i][j][k] * Vav_over_X*diff_t
					- temp_v[i][j][k] * Vav_over_Y*diff_t
					- temp_w[i][j][k] * Vav_over_Z*diff_t;
				v_nn[i][j][k] += -Aav*diff_t*(x_over_Y*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
					+ y_over_Y*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
					+ z_over_Y*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z);

				w_nn[i][j][k] = -coef1*(temp_w[i + 1][j][k] - temp_w[i][j][k])
					- coef2*(temp_w[i][j + 1][k] - temp_w[i][j][k])
					- coef3*(temp_w[i][j][k + 1] - temp_w[i][j][k]);
				w_nn[i][j][k] += -temp_u[i][j][k] * Wav_over_X*diff_t
					- temp_v[i][j][k] * Wav_over_Y*diff_t
					- temp_w[i][j][k] * Wav_over_Z*diff_t;
				w_nn[i][j][k] += -Aav*diff_t*(x_over_Z*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
					+ y_over_Z*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
					+ z_over_Z*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z);
				// PML X1:
				if (in < IMAX1){//&& kn>=KMAX1 && kn<KMAX-KMAX2
					absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);//sigma_x1
					absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

					u_nn[i][j][k] += -absorb_sigmax1*(temp_u[i][j][k])*diff_t + (1.0 - 1.0 / absorb_betax1)
						*(Aav*(x_over_X*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x +
						y_over_X*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y +
						z_over_X*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z))*diff_t;

					v_nn[i][j][k] += -absorb_sigmax1*(temp_v[i][j][k]
						+ Aav*(x_over_Y*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x +
						y_over_Y*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y +
						z_over_Y*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;

					w_nn[i][j][k] += -absorb_sigmax1*(temp_w[i][j][k]
						+ Aav*(x_over_Z*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x +
						y_over_Z*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y +
						z_over_Z*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;

				}
				// PML X2:
				if (in >= IMAX - IMAX2){//&& kn>=KMAX1 && kn<KMAX-KMAX2
					ii = i - (mpi_IMAX - mpi_IMAX2);
					absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
					absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

					u_nn[i][j][k] += -absorb_sigmax2*(temp_u[i][j][k])*diff_t + (1.0 - 1.0 / absorb_betax2)
						*(Aav*(x_over_X*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x +
						y_over_X*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y +
						z_over_X*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z))*diff_t;

					v_nn[i][j][k] += -absorb_sigmax2*(temp_v[i][j][k]
						+ Aav*(x_over_Y*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x +
						y_over_Y*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y +
						z_over_Y*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;

					w_nn[i][j][k] += -absorb_sigmax2*(temp_w[i][j][k]
						+ Aav*(x_over_Z*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x +
						y_over_Z*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y +
						z_over_Z*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;
				}
				// PML Y1 , Y1+X
				if (jn < JMAX1){
					absorb_sigmay1 = PML_AbsorbYmax*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);
					absorb_betay1 = 1.0 + 25.0*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);

					u_nn[i][j][k] += -absorb_sigmay1*(temp_u[i][j][k] +
						Aav*(x_over_X*(Qpy1[i][j][k] - Qpy1[i - 1][j][k])*diff_x +
						y_over_X*(Qpy1[i][j][k] - Qpy1[i][j - 1][k])*diff_y +
						z_over_X*(Qpy1[i][j][k] - Qpy1[i][j][k - 1])*diff_z))*diff_t;

					v_nn[i][j][k] += (1.0 - 1.0 / absorb_betay1)*Aav*(x_over_Y*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_Y*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_Y*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z)*diff_t
						- absorb_sigmay1*temp_v[i][j][k] * diff_t;

					w_nn[i][j][k] += -absorb_sigmay1*(temp_w[i][j][k] +
						Aav*(x_over_Z*(Qpy1[i][j][k] - Qpy1[i - 1][j][k])*diff_x +
						y_over_Z*(Qpy1[i][j][k] - Qpy1[i][j - 1][k])*diff_y +
						z_over_Z*(Qpy1[i][j][k] - Qpy1[i][j][k - 1])*diff_z))*diff_t;

					if (in < IMAX1){
						absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);//sigma_x1
						absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

						u_nn[i][j][k] += -(Qux1[i][j][k] * absorb_sigmax1*absorb_sigmay1
							+ (1 / absorb_betax1 - 1)*absorb_sigmay1*Aav*(x_over_X*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x
							+ y_over_X*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y
							+ z_over_X*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;

						v_nn[i][j][k] += -(Qvx1[i][j][k] * absorb_sigmax1*absorb_sigmay1
							+ (1 / absorb_betay1 - 1)*absorb_sigmax1*Aav*(x_over_Y*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x
							+ y_over_Y*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y
							+ z_over_Y*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;

						w_nn[i][j][k] += -Qwx1[i][j][k] * absorb_sigmax1*absorb_sigmay1*diff_t;
					}
					else if (in >= IMAX - IMAX2){
						ii = i - (mpi_IMAX - mpi_IMAX2);
						absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
						absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

						u_nn[i][j][k] += -(Qux2[ii][j][k] * absorb_sigmax2*absorb_sigmay1
							+ (1 / absorb_betax2 - 1)*absorb_sigmay1*Aav*(x_over_X*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x
							+ y_over_X*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y
							+ z_over_X*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;

						v_nn[i][j][k] += -(Qvx2[ii][j][k] * absorb_sigmax2*absorb_sigmay1
							+ (1 / absorb_betay1 - 1)*absorb_sigmax2*Aav*(x_over_Y*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x +
							y_over_Y*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y +
							z_over_Y*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;

						w_nn[i][j][k] += -Qwx2[ii][j][k] * absorb_sigmax2*absorb_sigmay1*diff_t;
					}
				}
				// PML Y2 , Y2+X
				if (jn >= JMAX - JMAX2){
					jj = j - (mpi_JMAX - mpi_JMAX2);
					absorb_sigmay2 = PML_AbsorbYmax*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);
					absorb_betay2 = 1.0 + 25.0*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);

					u_nn[i][j][k] += -absorb_sigmay2*(temp_u[i][j][k] +
						Aav*(x_over_X*(Qpy2[i][jj][k] - Qpy2[i - 1][jj][k])*diff_x +
						y_over_X*(Qpy2[i][jj][k] - Qpy2[i][jj - 1][k])*diff_y +
						z_over_X*(Qpy2[i][jj][k] - Qpy2[i][jj][k - 1])*diff_z))*diff_t;

					v_nn[i][j][k] += (1 - 1 / absorb_betay2)*Aav*(x_over_Y*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_Y*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_Y*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z)*diff_t
						- absorb_sigmay2*temp_v[i][j][k] * diff_t;

					w_nn[i][j][k] += -absorb_sigmay2*(temp_w[i][j][k] +
						Aav*(x_over_Z*(Qpy2[i][jj][k] - Qpy2[i - 1][jj][k])*diff_x +
						y_over_Z*(Qpy2[i][jj][k] - Qpy2[i][jj - 1][k])*diff_y +
						z_over_Z*(Qpy2[i][jj][k] - Qpy2[i][jj][k - 1])*diff_z))*diff_t;

					if (in < IMAX1){
						absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);//sigma_x1
						absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

						u_nn[i][j][k] += -(Qux1[i][j][k] * absorb_sigmax1*absorb_sigmay2
							+ (1.0 / absorb_betax1 - 1.0)*absorb_sigmay2*Aav*(x_over_X*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x
							+ y_over_X*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y
							+ z_over_X*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;

						v_nn[i][j][k] += -(Qvx1[i][j][k] * absorb_sigmax1*absorb_sigmay2
							+ (1.0 / absorb_betay2 - 1.0)*absorb_sigmax1*Aav*(x_over_Y*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x
							+ y_over_Y*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y
							+ z_over_Y*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;

						w_nn[i][j][k] += -Qwx1[i][j][k] * absorb_sigmax1*absorb_sigmay2*diff_t;
					}
					else if (in >= IMAX - IMAX2){
						ii = i - (mpi_IMAX - mpi_IMAX2);
						absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
						absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

						u_nn[i][j][k] += -(Qux2[ii][j][k] * absorb_sigmax2*absorb_sigmay2
							+ (1.0 / absorb_betax2 - 1.0)*absorb_sigmay2*Aav*(x_over_X*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x
							+ y_over_X*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y
							+ z_over_X*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;

						v_nn[i][j][k] += -(Qvx2[ii][j][k] * absorb_sigmax2*absorb_sigmay2
							+ (1.0 / absorb_betay2 - 1.0)*absorb_sigmax2*Aav*(x_over_Y*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x
							+ y_over_Y*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y
							+ z_over_Y*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;

						w_nn[i][j][k] += -Qwx2[ii][j][k] * absorb_sigmax2*absorb_sigmay2*diff_t;
					}/**/
				}
				//PML Z1
				if (kn < KMAX1){
					absorb_sigmaz1 = PML_AbsorbZmax*pow((PML_z1 - Z[k])*PML_Height1, PML_alphaz);
					absorb_betaz1 = 1.0 + 25.0*pow((PML_z1 - Z[k])*PML_Height1, PML_alphaz);

					u_nn[i][j][k] += -absorb_sigmaz1*(temp_u[i][j][k]
						+ Aav*(x_over_X*(Qpz1[i][j][k] - Qpz1[i - 1][j][k])*diff_x
						+ y_over_X*(Qpz1[i][j][k] - Qpz1[i][j - 1][k])*diff_y
						+ z_over_X*(Qpz1[i][j][k] - Qpz1[i][j][k - 1])*diff_z))*diff_t;

					v_nn[i][j][k] += -absorb_sigmaz1*(temp_v[i][j][k]
						+ Aav*(x_over_Y*(Qpz1[i][j][k] - Qpz1[i - 1][j][k])*diff_x
						+ y_over_Y*(Qpz1[i][j][k] - Qpz1[i][j - 1][k])*diff_y
						+ z_over_Y*(Qpz1[i][j][k] - Qpz1[i][j][k - 1])*diff_z))*diff_t;

					w_nn[i][j][k] += (1.0 - 1.0 / absorb_betaz1)*Aav*(x_over_Z*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_Z*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_Z*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z)*diff_t
						- absorb_sigmaz1*temp_w[i][j][k] * diff_t;

					if (in < IMAX1) {
						absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);//sigma_x1
						absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

						u_nn[i][j][k] += -(Qux1[i][j][k] * absorb_sigmax1*absorb_sigmaz1
							+ (1.0 / absorb_betax1 - 1.0)*absorb_sigmaz1*Aav*(x_over_X*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x
							+ y_over_X*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y
							+ z_over_X*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;

						v_nn[i][j][k] += -Qvx1[i][j][k] * absorb_sigmax1*absorb_sigmaz1*diff_t;

						w_nn[i][j][k] += -(Qwx1[i][j][k] * absorb_sigmax1*absorb_sigmaz1
							+ (1.0 / absorb_betaz1 - 1.0)*absorb_sigmax1*Aav*(x_over_Z*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x
							+ y_over_Z*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y
							+ z_over_Z*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;
					}
					else if (in >= IMAX - IMAX2){
						ii = i - (mpi_IMAX - mpi_IMAX2);
						absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
						absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

						u_nn[i][j][k] += -(Qux2[ii][j][k] * absorb_sigmax2*absorb_sigmaz1
							+ (1.0 / absorb_betax2 - 1.0)*absorb_sigmaz1*Aav*(x_over_X*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x
							+ y_over_X*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y
							+ z_over_X*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;

						v_nn[i][j][k] += -Qvx2[ii][j][k] * absorb_sigmax2*absorb_sigmaz1*diff_t;

						w_nn[i][j][k] += -(Qwx2[ii][j][k] * absorb_sigmax2*absorb_sigmaz1
							+ (1.0 / absorb_betaz1 - 1.0)*absorb_sigmax2*Aav*(x_over_Z*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x
							+ y_over_Z*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y
							+ z_over_Z*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;
					}

					if (jn < JMAX1){
						absorb_sigmay1 = PML_AbsorbYmax*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);
						absorb_betay1 = 1.0 + 25.0*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);

						u_nn[i][j][k] += -Quy1[i][j][k] * absorb_sigmay1*absorb_sigmaz1*diff_t;

						v_nn[i][j][k] += -(Qvy1[i][j][k] * absorb_sigmay1*absorb_sigmaz1
							+ (1.0 / absorb_betay1 - 1.0)*absorb_sigmaz1*Aav*(x_over_Y*(Qpy1[i][j][k] - Qpy1[i - 1][j][k])*diff_x
							+ y_over_Y*(Qpy1[i][j][k] - Qpy1[i][j - 1][k])*diff_y
							+ z_over_Y*(Qpy1[i][j][k] - Qpy1[i][j][k - 1])*diff_z))*diff_t;

						w_nn[i][j][k] += -(Qwy1[i][j][k] * absorb_sigmay1*absorb_sigmaz1
							+ (1.0 / absorb_betaz1 - 1.0)*absorb_sigmay1*Aav*(x_over_Z*(Qpy1[i][j][k] - Qpy1[i - 1][j][k])*diff_x
							+ y_over_Z*(Qpy1[i][j][k] - Qpy1[i][j - 1][k])*diff_y
							+ z_over_Z*(Qpy1[i][j][k] - Qpy1[i][j][k - 1])*diff_z))*diff_t;
					}
					else if (jn >= JMAX - JMAX2){
						jj = j - (mpi_JMAX - mpi_JMAX2);
						absorb_sigmay2 = PML_AbsorbYmax*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);
						absorb_betay2 = 1.0 + 25.0*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);

						u_nn[i][j][k] += -Quy2[i][jj][k] * absorb_sigmay2*absorb_sigmaz1*diff_t;

						v_nn[i][j][k] += -(Qvy2[i][jj][k] * absorb_sigmay2*absorb_sigmaz1
							+ (1.0 / absorb_betay2 - 1.0)*absorb_sigmaz1*Aav*(x_over_Y*(Qpy2[i][jj][k] - Qpy2[i - 1][jj][k])*diff_x
							+ y_over_Y*(Qpy2[i][jj][k] - Qpy2[i][jj - 1][k])*diff_y
							+ z_over_Y*(Qpy2[i][jj][k] - Qpy2[i][jj][k - 1])*diff_z))*diff_t;

						w_nn[i][j][k] += -(Qwy2[i][jj][k] * absorb_sigmay2*absorb_sigmaz1
							+ (1.0 / absorb_betaz1 - 1.0)*absorb_sigmay2*Aav*(x_over_Z*(Qpy2[i][jj][k] - Qpy2[i - 1][jj][k])*diff_x
							+ y_over_Z*(Qpy2[i][jj][k] - Qpy2[i][jj - 1][k])*diff_y
							+ z_over_Z*(Qpy2[i][jj][k] - Qpy2[i][jj][k - 1])*diff_z))*diff_t;/**/
					}
				}
				//PML Z2
				if (kn >= KMAX - KMAX2){
					kk = k - (mpi_KMAX - mpi_KMAX2);
					absorb_sigmaz2 = PML_AbsorbZmax*pow((Z[k] - PML_z2)*PML_Height2, PML_alphaz);
					absorb_betaz2 = 1.0 + 25.0*pow((Z[k] - PML_z2)*PML_Height2, PML_alphaz);

					u_nn[i][j][k] += -absorb_sigmaz2*(temp_u[i][j][k]
						+ Aav*(x_over_X*(Qpz2[i][j][kk] - Qpz2[i - 1][j][kk])*diff_x
						+ y_over_X*(Qpz2[i][j][kk] - Qpz2[i][j - 1][kk])*diff_y
						+ z_over_X*(Qpz2[i][j][kk] - Qpz2[i][j][kk - 1])*diff_z))*diff_t;

					v_nn[i][j][k] += -absorb_sigmaz2*(temp_v[i][j][k]
						+ Aav*(x_over_Y*(Qpz2[i][j][kk] - Qpz2[i - 1][j][kk])*diff_x
						+ y_over_Y*(Qpz2[i][j][kk] - Qpz2[i][j - 1][kk])*diff_y
						+ z_over_Y*(Qpz2[i][j][kk] - Qpz2[i][j][kk - 1])*diff_z))*diff_t;

					w_nn[i][j][k] += (1.0 - 1.0 / absorb_betaz2)*Aav*(x_over_Z*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_Z*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_Z*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z)*diff_t
						- absorb_sigmaz2*temp_w[i][j][k] * diff_t;

					if (in < IMAX1) {
						absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);//sigma_x1
						absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

						u_nn[i][j][k] += -(Qux1[i][j][k] * absorb_sigmax1*absorb_sigmaz2
							+ (1.0 / absorb_betax1 - 1.0)*absorb_sigmaz2*Aav*(x_over_X*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x
							+ y_over_X*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y
							+ z_over_X*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;

						v_nn[i][j][k] += -Qvx1[i][j][k] * absorb_sigmax1*absorb_sigmaz2*diff_t;

						w_nn[i][j][k] += -(Qwx1[i][j][k] * absorb_sigmax1*absorb_sigmaz2
							+ (1.0 / absorb_betaz2 - 1.0)*absorb_sigmax1*Aav*(x_over_Z*(Qpx1[i][j][k] - Qpx1[i - 1][j][k])*diff_x
							+ y_over_Z*(Qpx1[i][j][k] - Qpx1[i][j - 1][k])*diff_y
							+ z_over_Z*(Qpx1[i][j][k] - Qpx1[i][j][k - 1])*diff_z))*diff_t;
					}
					else if (in >= IMAX - IMAX2){
						ii = i - (mpi_IMAX - mpi_IMAX2);
						absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
						absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

						u_nn[i][j][k] += -(Qux2[ii][j][k] * absorb_sigmax2*absorb_sigmaz2
							+ (1.0 / absorb_betax2 - 1.0)*absorb_sigmaz2*Aav*(x_over_X*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x
							+ y_over_X*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y
							+ z_over_X*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;

						v_nn[i][j][k] += -Qvx2[ii][j][k] * absorb_sigmax2*absorb_sigmaz2*diff_t;

						w_nn[i][j][k] += -(Qwx2[ii][j][k] * absorb_sigmax2*absorb_sigmaz2
							+ (1.0 / absorb_betaz2 - 1.0)*absorb_sigmax2*Aav*(x_over_Z*(Qpx2[ii][j][k] - Qpx2[ii - 1][j][k])*diff_x
							+ y_over_Z*(Qpx2[ii][j][k] - Qpx2[ii][j - 1][k])*diff_y
							+ z_over_Z*(Qpx2[ii][j][k] - Qpx2[ii][j][k - 1])*diff_z))*diff_t;
					}

					if (jn < JMAX1){
						absorb_sigmay1 = PML_AbsorbYmax*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);
						absorb_betay1 = 1.0 + 25.0*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);

						u_nn[i][j][k] += -Quy1[i][j][k] * absorb_sigmay1*absorb_sigmaz2*diff_t;

						v_nn[i][j][k] += -(Qvy1[i][j][k] * absorb_sigmay1*absorb_sigmaz2
							+ (1.0 / absorb_betay1 - 1.0)*absorb_sigmaz2*Aav*(x_over_Y*(Qpy1[i][j][k] - Qpy1[i - 1][j][k])*diff_x
							+ y_over_Y*(Qpy1[i][j][k] - Qpy1[i][j - 1][k])*diff_y
							+ z_over_Y*(Qpy1[i][j][k] - Qpy1[i][j][k - 1])*diff_z))*diff_t;

						w_nn[i][j][k] += -(Qwy1[i][j][k] * absorb_sigmay1*absorb_sigmaz2
							+ (1.0 / absorb_betaz2 - 1.0)*absorb_sigmay1*Aav*(x_over_Z*(Qpy1[i][j][k] - Qpy1[i - 1][j][k])*diff_x
							+ y_over_Z*(Qpy1[i][j][k] - Qpy1[i][j - 1][k])*diff_y
							+ z_over_Z*(Qpy1[i][j][k] - Qpy1[i][j][k - 1])*diff_z))*diff_t;
					}
					else if (jn >= JMAX - JMAX2){
						jj = j - (mpi_JMAX - mpi_JMAX2);
						absorb_sigmay2 = PML_AbsorbYmax*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);
						absorb_betay2 = 1.0 + 25.0*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);

						u_nn[i][j][k] += -Quy2[i][jj][k] * absorb_sigmay2*absorb_sigmaz2*diff_t;

						v_nn[i][j][k] += -(Qvy2[i][jj][k] * absorb_sigmay2*absorb_sigmaz2
							+ (1.0 / absorb_betay2 - 1.0)*absorb_sigmaz2*Aav*(x_over_Y*(Qpy2[i][jj][k] - Qpy2[i - 1][jj][k])*diff_x
							+ y_over_Y*(Qpy2[i][jj][k] - Qpy2[i][jj - 1][k])*diff_y
							+ z_over_Y*(Qpy2[i][jj][k] - Qpy2[i][jj][k - 1])*diff_z))*diff_t;

						w_nn[i][j][k] += -(Qwy2[i][jj][k] * absorb_sigmay2*absorb_sigmaz2
							+ (1.0 / absorb_betaz2 - 1.0)*absorb_sigmay2*Aav*(x_over_Z*(Qpy2[i][jj][k] - Qpy2[i - 1][jj][k])*diff_x
							+ y_over_Z*(Qpy2[i][jj][k] - Qpy2[i][jj - 1][k])*diff_y
							+ z_over_Z*(Qpy2[i][jj][k] - Qpy2[i][jj][k - 1])*diff_z))*diff_t;
					}/**/
				}

				//calculation inside hill part 
				if (cr_judge != 0 && (in >= istar && in <= istop) && (jn >= jstar && jn <= jstop)
					&& (kn >= kstar && kn <= kstop) && cr_uc[in - istar][jn - jstar][kn - kstar] == 1){
					u_nn[i][j][k] = -1.0*crpo_invBeta1*diff_t*crpo_eff_density1*
						(x_over_X*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_X*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_X*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z) +
						temp_u[i][j][k] * crpo_invBeta1 - temp_u[i][j][k];
					v_nn[i][j][k] = -1.0*crpo_invBeta1*diff_t*crpo_eff_density1*
						(x_over_Y*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_Y*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_Y*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z) +
						temp_v[i][j][k] * crpo_invBeta1 - temp_v[i][j][k];

					w_nn[i][j][k] = -1.0*crpo_invBeta1*diff_t*crpo_eff_density1*
						(x_over_Z*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_Z*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_Z*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z) +
						temp_w[i][j][k] * crpo_invBeta1 - temp_w[i][j][k];
				}
				if (cr_judge != 0 && (in >= istar && in <= istop) && (jn >= jstar && jn <= jstop)
					&& (kn >= kstar && kn <= kstop) && cr_uc[in - istar][jn - jstar][kn - kstar] == 2){
					u_nn[i][j][k] = -1.0*crpo_invBeta2*diff_t*crpo_eff_density2*
						(x_over_X*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_X*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_X*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z) +
						temp_u[i][j][k] * crpo_invBeta2 - temp_u[i][j][k];
					v_nn[i][j][k] = -1.0*crpo_invBeta2*diff_t*crpo_eff_density2*
						(x_over_Y*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_Y*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_Y*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z) +
						temp_v[i][j][k] * crpo_invBeta2 - temp_v[i][j][k];

					w_nn[i][j][k] = -1.0*crpo_invBeta2*diff_t*crpo_eff_density2*
						(x_over_Z*(temp_p[i][j][k] - temp_p[i - 1][j][k])*diff_x
						+ y_over_Z*(temp_p[i][j][k] - temp_p[i][j - 1][k])*diff_y
						+ z_over_Z*(temp_p[i][j][k] - temp_p[i][j][k - 1])*diff_z) +
						temp_w[i][j][k] * crpo_invBeta2 - temp_w[i][j][k];
				}
			}
		}
	}
}


//------------------------------calculate pressure (p) for new time(nn)------------//
void air::cal_fp(double ***temp_u, double ***temp_v, double ***temp_w,
	double ***temp_p)
{
	int i, j, k, ii, jj, kk;
	int in, jn, kn;
	//  PML boundary for z direction
	double absorb_sigmaz1, absorb_sigmaz2, temp_z, absorb_betaz1, absorb_betaz2, absorb_coefzx, absorb_coefzy;
	//  PML boundary for x direction
	double absorb_sigmax1, absorb_sigmax2, temp_x, absorb_betax1, absorb_betax2, absorb_coefxz, absorb_coefxy;
	//  PML boundary for y direction
	double absorb_sigmay1, absorb_sigmay2, temp_y, absorb_betay1, absorb_betay2, absorb_coefyx, absorb_coefyz;

	if (scheme == FB_p_v && (PML_judgey == 1 || PML_judgey == 3)){
		PML_JMAX1 = PML_JMAX1 + 1;
		if (PML_JMAX1 >= PML_JMAX1_time)JMAX1 = JMAX3;
		else JMAX1 = 0;
	}

	for (i = 1; i < mpi_IMAX - 1; i++){
		for (j = 1; j < mpi_JMAX - 1; j++){
			for (k = 1; k < mpi_KMAX - 1; k++){
				in = mpi_i1 + i; jn = mpi_j1 + j; kn = mpi_k1 + k;
				x_over_X = Geo_over_var(Y_over_y[j], Z_over_z[k], Y_over_z,
					Z_over_y, Jacobian[i][j][k]);
				x_over_Y = Geo_over_var(X_over_z, Z_over_y, X_over_y,
					Z_over_z[k], Jacobian[i][j][k]);
				x_over_Z = Geo_over_var(X_over_y, Y_over_z, X_over_z,
					Y_over_y[j], Jacobian[i][j][k]);
				y_over_X = Geo_over_var(Y_over_z, Z_over_x, Y_over_x,
					Z_over_z[k], Jacobian[i][j][k]);
				y_over_Y = Geo_over_var(X_over_x[i], Z_over_z[k], X_over_z,
					Z_over_x, Jacobian[i][j][k]);
				y_over_Z = Geo_over_var(X_over_z, Y_over_x, X_over_x[i],
					Y_over_z, Jacobian[i][j][k]);
				z_over_X = Geo_over_var(Y_over_x, Z_over_y, Y_over_y[j],
					Z_over_x, Jacobian[i][j][k]);
				z_over_Y = Geo_over_var(X_over_y, Z_over_x, X_over_x[i],
					Z_over_y, Jacobian[i][j][k]);
				z_over_Z = Geo_over_var(X_over_x[i], Y_over_y[j], X_over_y,
					Y_over_x, Jacobian[i][j][k]);
				coef1 = (Uav[i][j][k] * x_over_X + Vav[i][j][k] * x_over_Y
					+ Wav[i][j][k] * x_over_Z)*diff_x*diff_t;
				coef2 = (Uav[i][j][k] * y_over_X + Vav[i][j][k] * y_over_Y
					+ Wav[i][j][k] * y_over_Z)*diff_y*diff_t;
				coef3 = (Uav[i][j][k] * z_over_X + Vav[i][j][k] * z_over_Y
					+ Wav[i][j][k] * z_over_Z)*diff_z*diff_t;
				Uav_over_X = x_over_X*(Uav[i + 1][j][k] - Uav[i - 1][j][k])*diff_x*0.5 +
					y_over_X*(Uav[i][j + 1][k] - Uav[i][j - 1][k])*diff_y*0.5 +
					z_over_X*(Uav[i][j][k + 1] - Uav[i][j][k - 1])*diff_z*0.5;
				Uav_over_Y = x_over_Y*(Uav[i + 1][j][k] - Uav[i - 1][j][k])*diff_x*0.5 +
					y_over_Y*(Uav[i][j + 1][k] - Uav[i][j - 1][k])*diff_y*0.5 +
					z_over_Y*(Uav[i][j][k + 1] - Uav[i][j][k - 1])*diff_z*0.5;
				Uav_over_Z = x_over_Z*(Uav[i + 1][j][k] - Uav[i - 1][j][k])*diff_x*0.5 +
					y_over_Z*(Uav[i][j + 1][k] - Uav[i][j - 1][k])*diff_y*0.5 +
					z_over_Z*(Uav[i][j][k + 1] - Uav[i][j][k - 1])*diff_z*0.5;
				Vav_over_X = x_over_X*(Vav[i + 1][j][k] - Vav[i - 1][j][k])*diff_x*0.5 +
					y_over_X*(Vav[i][j + 1][k] - Vav[i][j - 1][k])*diff_y*0.5 +
					z_over_X*(Vav[i][j][k + 1] - Vav[i][j][k - 1])*diff_z*0.5;
				Vav_over_Y = x_over_Y*(Vav[i + 1][j][k] - Vav[i - 1][j][k])*diff_x*0.5 +
					y_over_Y*(Vav[i][j + 1][k] - Vav[i][j - 1][k])*diff_y*0.5 +
					z_over_Y*(Vav[i][j][k + 1] - Vav[i][j][k - 1])*diff_z*0.5;
				Vav_over_Z = x_over_Z*(Vav[i + 1][j][k] - Vav[i - 1][j][k])*diff_x*0.5 +
					y_over_Z*(Vav[i][j + 1][k] - Vav[i][j - 1][k])*diff_y*0.5 +
					z_over_Z*(Vav[i][j][k + 1] - Vav[i][j][k - 1])*diff_z*0.5;
				Wav_over_X = x_over_X*(Wav[i + 1][j][k] - Wav[i - 1][j][k])*diff_x*0.5 +
					y_over_X*(Wav[i][j + 1][k] - Wav[i][j - 1][k])*diff_y*0.5 +
					z_over_X*(Wav[i][j][k + 1] - Wav[i][j][k - 1])*diff_z*0.5;
				Wav_over_Y = x_over_Y*(Wav[i + 1][j][k] - Wav[i - 1][j][k])*diff_x*0.5 +
					y_over_Y*(Wav[i][j + 1][k] - Wav[i][j - 1][k])*diff_y*0.5 +
					z_over_Y*(Wav[i][j][k + 1] - Wav[i][j][k - 1])*diff_z*0.5;
				Wav_over_Z = x_over_Z*(Wav[i + 1][j][k] - Wav[i - 1][j][k])*diff_x*0.5 +
					y_over_Z*(Wav[i][j + 1][k] - Wav[i][j - 1][k])*diff_y*0.5 +
					z_over_Z*(Wav[i][j][k + 1] - Wav[i][j][k - 1])*diff_z*0.5;

				p_nn[i][j][k] = -coef1*(temp_p[i][j][k] - temp_p[i - 1][j][k])
					- coef2*(temp_p[i][j][k] - temp_p[i][j - 1][k])
					- coef3*(temp_p[i][j][k] - temp_p[i][j][k - 1]);
				p_nn[i][j][k] += -adiabatic_coef*temp_p[i][j][k] * diff_t*(Uav_over_X + Vav_over_Y + Wav_over_Z);

				p_nn[i][j][k] += -adiabatic_coef*Pav*diff_t*(x_over_X*(temp_u[i + 1][j][k] - temp_u[i][j][k])*diff_x
					+ y_over_X*(temp_u[i][j + 1][k] - temp_u[i][j][k])*diff_y
					+ z_over_X*(temp_u[i][j][k + 1] - temp_u[i][j][k])*diff_z
					+ x_over_Y*(temp_v[i + 1][j][k] - temp_v[i][j][k])*diff_x
					+ y_over_Y*(temp_v[i][j + 1][k] - temp_v[i][j][k])*diff_y
					+ z_over_Y*(temp_v[i][j][k + 1] - temp_v[i][j][k])*diff_z
					+ x_over_Z*(temp_w[i + 1][j][k] - temp_w[i][j][k])*diff_x
					+ y_over_Z*(temp_w[i][j + 1][k] - temp_w[i][j][k])*diff_y
					+ z_over_Z*(temp_w[i][j][k + 1] - temp_w[i][j][k])*diff_z);

				// calculate PML boundary for x,y,z
				// PML X1
				if (in < IMAX1){//&& kn>=KMAX1 && kn<KMAX-KMAX2
					absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);
					absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

					p_nn[i][j][k] += (1.0 - 1.0 / absorb_betax1)*(adiabatic_coef*Pav*(x_over_X*(temp_u[i + 1][j][k] - temp_u[i][j][k])*diff_x
						+ y_over_X*(temp_u[i][j + 1][k] - temp_u[i][j][k])*diff_y
						+ z_over_X*(temp_u[i][j][k + 1] - temp_u[i][j][k])*diff_z))*diff_t
						- (absorb_sigmax1 * (temp_p[i][j][k] + adiabatic_coef*Pav*(x_over_Y*(Qvx1[i + 1][j][k] - Qvx1[i][j][k])*diff_x +
						y_over_Y*(Qvx1[i][j + 1][k] - Qvx1[i][j][k])*diff_y +
						z_over_Y*(Qvx1[i][j][k + 1] - Qvx1[i][j][k])*diff_z +
						x_over_Z*(Qwx1[i + 1][j][k] - Qwx1[i][j][k])*diff_x +
						y_over_Z*(Qwx1[i][j + 1][k] - Qwx1[i][j][k])*diff_y +
						z_over_Z*(Qwx1[i][j][k + 1] - Qwx1[i][j][k])*diff_z)))*diff_t;
				}
				// PML X2
				if (in >= IMAX - IMAX2){//&& kn>=KMAX1 && kn<KMAX-KMAX2
					ii = i - (mpi_IMAX - mpi_IMAX2);
					absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
					absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

					p_nn[i][j][k] += (1.0 - 1.0 / absorb_betax2)*(adiabatic_coef*Pav*(x_over_X*(temp_u[i + 1][j][k] - temp_u[i][j][k])*diff_x
						+ y_over_X*(temp_u[i][j + 1][k] - temp_u[i][j][k])*diff_y
						+ z_over_X*(temp_u[i][j][k + 1] - temp_u[i][j][k])*diff_z))*diff_t
						- (absorb_sigmax2 * (temp_p[i][j][k]
						+ adiabatic_coef*Pav*(x_over_Y*(Qvx2[ii + 1][j][k] - Qvx2[ii][j][k])*diff_x +
						y_over_Y*(Qvx2[ii][j + 1][k] - Qvx2[ii][j][k])*diff_y +
						z_over_Y*(Qvx2[ii][j][k + 1] - Qvx2[ii][j][k])*diff_z +
						x_over_Z*(Qwx2[ii + 1][j][k] - Qwx2[ii][j][k])*diff_x +
						y_over_Z*(Qwx2[ii][j + 1][k] - Qwx2[ii][j][k])*diff_y +
						z_over_Z*(Qwx2[ii][j][k + 1] - Qwx2[ii][j][k])*diff_z)))*diff_t;
				}
				// PML Y1 Y1+X
				if (jn < JMAX1){
					absorb_sigmay1 = PML_AbsorbYmax*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);
					absorb_betay1 = 1.0 + 25.0*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);

					p_nn[i][j][k] += (1.0 - 1.0 / absorb_betay1)*(adiabatic_coef*Pav*(x_over_Y*(temp_v[i + 1][j][k] - temp_v[i][j][k])*diff_x
						+ y_over_Y*(temp_v[i][j + 1][k] - temp_v[i][j][k])*diff_y
						+ z_over_Y*(temp_v[i][j][k + 1] - temp_v[i][j][k])*diff_z))*diff_t
						- absorb_sigmay1*(temp_p[i][j][k]
						+ adiabatic_coef*Pav*(x_over_X*(Quy1[i + 1][j][k] - Quy1[i][j][k])*diff_x
						+ y_over_X*(Quy1[i][j + 1][k] - Quy1[i][j][k])*diff_y
						+ z_over_X*(Quy1[i][j][k + 1] - Quy1[i][j][k])*diff_z
						+ x_over_Z*(Qwy1[i + 1][j][k] - Qwy1[i][j][k])*diff_x
						+ y_over_Z*(Qwy1[i][j + 1][k] - Qwy1[i][j][k])*diff_y
						+ z_over_Z*(Qwy1[i][j][k + 1] - Qwy1[i][j][k])*diff_z))*diff_t;

					if (in < IMAX1){
						absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);
						absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

						p_nn[i][j][k] += -(Qpx1[i][j][k] * absorb_sigmax1*absorb_sigmay1 + adiabatic_coef*Pav*(
							(1 / absorb_betax1 - 1)*absorb_sigmay1*(x_over_X*(Qux1[i + 1][j][k] - Qux1[i][j][k])*diff_x
							+ y_over_X*(Qux1[i][j + 1][k] - Qux1[i][j][k])*diff_y
							+ z_over_X*(Qux1[i][j][k + 1] - Qux1[i][j][k])*diff_z)
							+ (1 / absorb_betay1 - 1)*absorb_sigmax1*(x_over_Y*(Qvx1[i + 1][j][k] - Qvx1[i][j][k])*diff_x
							+ y_over_Y*(Qvx1[i][j + 1][k] - Qvx1[i][j][k])*diff_y
							+ z_over_Y*(Qvx1[i][j][k + 1] - Qvx1[i][j][k])*diff_z)))*diff_t;
					}
					else if (in >= IMAX - IMAX2){
						ii = i - (mpi_IMAX - mpi_IMAX2);
						absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
						absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

						p_nn[i][j][k] += -(Qpx2[ii][j][k] * absorb_sigmax2*absorb_sigmay1 + adiabatic_coef*Pav*(
							(1 / absorb_betax2 - 1)*absorb_sigmay1*(x_over_X*(Qux2[ii + 1][j][k] - Qux2[ii][j][k])*diff_x
							+ y_over_X*(Qux2[ii][j + 1][k] - Qux2[ii][j][k])*diff_y
							+ z_over_X*(Qux2[ii][j][k + 1] - Qux2[ii][j][k])*diff_z)
							+ (1 / absorb_betay1 - 1)*absorb_sigmax2*(x_over_Y*(Qvx2[ii + 1][j][k] - Qvx2[ii][j][k])*diff_x
							+ y_over_Y*(Qvx2[ii][j + 1][k] - Qvx2[ii][j][k])*diff_y
							+ z_over_Y*(Qvx2[ii][j][k + 1] - Qvx2[ii][j][k])*diff_z)))*diff_t;
					}/**/
				}
				// PML Y1 Y1+X
				if (jn >= JMAX - JMAX2){
					jj = j - (mpi_JMAX - mpi_JMAX2);
					absorb_sigmay2 = PML_AbsorbYmax*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);
					absorb_betay2 = 1.0 + 25.0*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);

					p_nn[i][j][k] += (1 - 1 / absorb_betay2)*(adiabatic_coef*Pav*(x_over_Y*(temp_v[i + 1][j][k] - temp_v[i][j][k])*diff_x
						+ y_over_Y*(temp_v[i][j + 1][k] - temp_v[i][j][k])*diff_y
						+ z_over_Y*(temp_v[i][j][k + 1] - temp_v[i][j][k])*diff_z))*diff_t
						- absorb_sigmay2*(temp_p[i][j][k] + adiabatic_coef*Pav*(x_over_X*(Quy2[i + 1][jj][k] - Quy2[i][jj][k])*diff_x
						+ y_over_X*(Quy2[i][jj + 1][k] - Quy2[i][jj][k])*diff_y
						+ z_over_X*(Quy2[i][jj][k + 1] - Quy2[i][jj][k])*diff_z
						+ x_over_Z*(Qwy2[i + 1][jj][k] - Qwy2[i][jj][k])*diff_x
						+ y_over_Z*(Qwy2[i][jj + 1][k] - Qwy2[i][jj][k])*diff_y
						+ z_over_Z*(Qwy2[i][jj][k + 1] - Qwy2[i][jj][k])*diff_z))*diff_t;

					if (in < IMAX1){
						absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);
						absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

						p_nn[i][j][k] += -(Qpx1[i][j][k] * absorb_sigmax1*absorb_sigmay2 + adiabatic_coef*Pav*(
							(1 / absorb_betax1 - 1)*absorb_sigmay2*(x_over_X*(Qux1[i + 1][j][k] - Qux1[i][j][k])*diff_x
							+ y_over_X*(Qux1[i][j + 1][k] - Qux1[i][j][k])*diff_y
							+ z_over_X*(Qux1[i][j][k + 1] - Qux1[i][j][k])*diff_z)
							+ (1 / absorb_betay2 - 1)*absorb_sigmax1*(x_over_Y*(Qvx1[i + 1][j][k] - Qvx1[i][j][k])*diff_x
							+ y_over_Y*(Qvx1[i][j + 1][k] - Qvx1[i][j][k])*diff_y
							+ z_over_Y*(Qvx1[i][j][k + 1] - Qvx1[i][j][k])*diff_z)))*diff_t;
					}
					else if (in >= IMAX - IMAX2){
						ii = i - (mpi_IMAX - mpi_IMAX2);
						absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
						absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

						p_nn[i][j][k] += -(Qpx2[ii][j][k] * absorb_sigmax2*absorb_sigmay2 + adiabatic_coef*Pav*(
							(1 / absorb_betax2 - 1)*absorb_sigmay2*(x_over_X*(Qux2[ii + 1][j][k] - Qux2[ii][j][k])*diff_x
							+ y_over_X*(Qux2[ii][j + 1][k] - Qux2[ii][j][k])*diff_y
							+ z_over_X*(Qux2[ii][j][k + 1] - Qux2[ii][j][k])*diff_z)
							+ (1 / absorb_betay2 - 1)*absorb_sigmax2*(x_over_Y*(Qvx2[ii + 1][j][k] - Qvx2[ii][j][k])*diff_x
							+ y_over_Y*(Qvx2[ii][j + 1][k] - Qvx2[ii][j][k])*diff_y
							+ z_over_Y*(Qvx2[ii][j][k + 1] - Qvx2[ii][j][k])*diff_z)))*diff_t;
					}
				}
				// PML Z1 
				if (kn < KMAX1){
					absorb_sigmaz1 = PML_AbsorbZmax*pow((PML_z1 - Z[k])*PML_Height1, PML_alphaz);
					absorb_betaz1 = 1.0 + 25.0*pow((PML_z1 - Z[k])*PML_Height1, PML_alphaz);

					p_nn[i][j][k] += (1.0 - 1.0 / absorb_betaz1)*adiabatic_coef*Pav*(x_over_Z*(temp_w[i + 1][j][k] - temp_w[i][j][k])*diff_x
						+ y_over_Z*(temp_w[i][j + 1][k] - temp_w[i][j][k])*diff_y
						+ z_over_Z*(temp_w[i][j][k + 1] - temp_w[i][j][k])*diff_z)*diff_t
						- absorb_sigmaz1*(temp_p[i][j][k]
						+ adiabatic_coef*Pav*(x_over_X*(Quz1[i + 1][j][k] - Quz1[i][j][k])*diff_x
						+ y_over_X*(Quz1[i][j + 1][k] - Quz1[i][j][k])*diff_y
						+ z_over_X*(Quz1[i][j][k + 1] - Quz1[i][j][k])*diff_z
						+ x_over_Y*(Qvz1[i + 1][j][k] - Qvz1[i][j][k])*diff_x
						+ y_over_Y*(Qvz1[i][j + 1][k] - Qvz1[i][j][k])*diff_y
						+ z_over_Y*(Qvz1[i][j][k + 1] - Qvz1[i][j][k])*diff_z))*diff_t;

					if (in < IMAX1){
						absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);
						absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

						p_nn[i][j][k] += -(Qpx1[i][j][k] * absorb_sigmax1*absorb_sigmaz1 + adiabatic_coef*Pav*(
							(1 / absorb_betax1 - 1)*absorb_sigmaz1*(x_over_X*(Qux1[i + 1][j][k] - Qux1[i][j][k])*diff_x
							+ y_over_X*(Qux1[i][j + 1][k] - Qux1[i][j][k])*diff_y
							+ z_over_X*(Qux1[i][j][k + 1] - Qux1[i][j][k])*diff_z)
							+ (1 / absorb_betaz1 - 1)*absorb_sigmax1*(x_over_Z*(Qwx1[i + 1][j][k] - Qwx1[i][j][k])*diff_x
							+ y_over_Z*(Qwx1[i][j + 1][k] - Qwx1[i][j][k])*diff_y
							+ z_over_Z*(Qwx1[i][j][k + 1] - Qwx1[i][j][k])*diff_z)))*diff_t;
					}
					else if (in >= IMAX - IMAX2){
						ii = i - (mpi_IMAX - mpi_IMAX2);
						absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
						absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

						p_nn[i][j][k] += -(Qpx2[ii][j][k] * absorb_sigmax2*absorb_sigmaz1 + adiabatic_coef*Pav*(
							(1 / absorb_betax2 - 1)*absorb_sigmaz1*(x_over_X*(Qux2[ii + 1][j][k] - Qux2[ii][j][k])*diff_x
							+ y_over_X*(Qux2[ii][j + 1][k] - Qux2[ii][j][k])*diff_y
							+ z_over_X*(Qux2[ii][j][k + 1] - Qux2[ii][j][k])*diff_z)
							+ (1 / absorb_betaz1 - 1)*absorb_sigmax2*(x_over_Z*(Qwx2[ii + 1][j][k] - Qwx2[ii][j][k])*diff_x
							+ y_over_Z*(Qwx2[ii][j + 1][k] - Qwx2[ii][j][k])*diff_y
							+ z_over_Z*(Qwx2[ii][j][k + 1] - Qwx2[ii][j][k])*diff_z)))*diff_t;
					}
					if (jn < JMAX1){
						absorb_sigmay1 = PML_AbsorbYmax*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);
						absorb_betay1 = 1.0 + 25.0*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);

						p_nn[i][j][k] += -(Qpy1[i][j][k] * absorb_sigmay1*absorb_sigmaz1 + adiabatic_coef*Pav*(
							(1 / absorb_betay1 - 1)*absorb_sigmaz1*(x_over_Y*(Qvy1[i + 1][j][k] - Qvy1[i][j][k])*diff_x
							+ y_over_Y*(Qvy1[i][j + 1][k] - Qvy1[i][j][k])*diff_y
							+ z_over_Y*(Qvy1[i][j][k + 1] - Qvy1[i][j][k])*diff_z)
							+ (1 / absorb_betaz1 - 1)*absorb_sigmay1*(x_over_Z*(Qwy1[i + 1][j][k] - Qwy1[i][j][k])*diff_x
							+ y_over_Z*(Qwy1[i][j + 1][k] - Qwy1[i][j][k])*diff_y
							+ z_over_Z*(Qwy1[i][j][k + 1] - Qwy1[i][j][k])*diff_z)))*diff_t;
					}
					else if (jn >= JMAX - JMAX2){
						jj = j - (mpi_JMAX - mpi_JMAX2);
						absorb_sigmay2 = PML_AbsorbYmax*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);
						absorb_betay2 = 1.0 + 25.0*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);

						p_nn[i][j][k] += -(Qpy2[i][jj][k] * absorb_sigmay2*absorb_sigmaz1 + adiabatic_coef*Pav*(
							(1 / absorb_betay2 - 1)*absorb_sigmaz1*(x_over_Y*(Qvy2[i + 1][jj][k] - Qvy2[i][jj][k])*diff_x
							+ y_over_Y*(Qvy2[i][jj + 1][k] - Qvy2[i][jj][k])*diff_y
							+ z_over_Y*(Qvy2[i][jj][k + 1] - Qvy2[i][jj][k])*diff_z)
							+ (1 / absorb_betaz1 - 1)*absorb_sigmay2*(x_over_Z*(Qwy2[i + 1][jj][k] - Qwy2[i][jj][k])*diff_x
							+ y_over_Z*(Qwy2[i][jj + 1][k] - Qwy2[i][jj][k])*diff_y
							+ z_over_Z*(Qwy2[i][jj][k + 1] - Qwy2[i][jj][k])*diff_z)))*diff_t;
					}/**/
				}
				// PML Z2 
				if (kn >= KMAX - KMAX2){
					kk = k - (mpi_KMAX - mpi_KMAX2);
					absorb_sigmaz2 = PML_AbsorbZmax*pow((Z[k] - PML_z2)*PML_Height2, PML_alphaz);
					absorb_betaz2 = 1.0 + 25.0*pow((Z[k] - PML_z2)*PML_Height2, PML_alphaz);

					p_nn[i][j][k] += (1.0 - 1.0 / absorb_betaz2)*adiabatic_coef*Pav*(x_over_Z*(temp_w[i + 1][j][k] - temp_w[i][j][k])*diff_x
						+ y_over_Z*(temp_w[i][j + 1][k] - temp_w[i][j][k])*diff_y
						+ z_over_Z*(temp_w[i][j][k + 1] - temp_w[i][j][k])*diff_z)*diff_t
						- absorb_sigmaz2*(temp_p[i][j][k]
						+ adiabatic_coef*Pav*(x_over_X*(Quz2[i + 1][j][kk] - Quz2[i][j][kk])*diff_x
						+ y_over_X*(Quz2[i][j + 1][kk] - Quz2[i][j][kk])*diff_y
						+ z_over_X*(Quz2[i][j][kk + 1] - Quz2[i][j][kk])*diff_z
						+ x_over_Y*(Qvz2[i + 1][j][kk] - Qvz2[i][j][kk])*diff_x
						+ y_over_Y*(Qvz2[i][j + 1][kk] - Qvz2[i][j][kk])*diff_y
						+ z_over_Y*(Qvz2[i][j][kk + 1] - Qvz2[i][j][kk])*diff_z))*diff_t;
					if (in < IMAX1){
						absorb_sigmax1 = PML_AbsorbXmax*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);
						absorb_betax1 = 1.0 + 25.0*pow((PML_x1 - X[i])*PML_Width1, PML_alphax);

						p_nn[i][j][k] += -(Qpx1[i][j][k] * absorb_sigmax1*absorb_sigmaz2 + adiabatic_coef*Pav*(
							(1 / absorb_betax1 - 1)*absorb_sigmaz2*(x_over_X*(Qux1[i + 1][j][k] - Qux1[i][j][k])*diff_x
							+ y_over_X*(Qux1[i][j + 1][k] - Qux1[i][j][k])*diff_y
							+ z_over_X*(Qux1[i][j][k + 1] - Qux1[i][j][k])*diff_z)
							+ (1 / absorb_betaz2 - 1)*absorb_sigmax1*(x_over_Z*(Qwx1[i + 1][j][k] - Qwx1[i][j][k])*diff_x
							+ y_over_Z*(Qwx1[i][j + 1][k] - Qwx1[i][j][k])*diff_y
							+ z_over_Z*(Qwx1[i][j][k + 1] - Qwx1[i][j][k])*diff_z)))*diff_t;
					}
					else if (in >= IMAX - IMAX2){
						ii = i - (mpi_IMAX - mpi_IMAX2);
						absorb_sigmax2 = PML_AbsorbXmax*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);
						absorb_betax2 = 1.0 + 25.0*pow((X[i] - PML_x2)*PML_Width2, PML_alphax);

						p_nn[i][j][k] += -(Qpx2[ii][j][k] * absorb_sigmax2*absorb_sigmaz2 + adiabatic_coef*Pav*(
							(1 / absorb_betax2 - 1)*absorb_sigmaz2*(x_over_X*(Qux2[ii + 1][j][k] - Qux2[ii][j][k])*diff_x
							+ y_over_X*(Qux2[ii][j + 1][k] - Qux2[ii][j][k])*diff_y
							+ z_over_X*(Qux2[ii][j][k + 1] - Qux2[ii][j][k])*diff_z)
							+ (1 / absorb_betaz2 - 1)*absorb_sigmax2*(x_over_Z*(Qwx2[ii + 1][j][k] - Qwx2[ii][j][k])*diff_x
							+ y_over_Z*(Qwx2[ii][j + 1][k] - Qwx2[ii][j][k])*diff_y
							+ z_over_Z*(Qwx2[ii][j][k + 1] - Qwx2[ii][j][k])*diff_z)))*diff_t;
					}
					if (jn < JMAX1){
						absorb_sigmay1 = PML_AbsorbYmax*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);
						absorb_betay1 = 1.0 + 25.0*pow((PML_y1 - Y[j])*PML_Length1, PML_alphay);

						p_nn[i][j][k] += -(Qpy1[i][j][k] * absorb_sigmay1*absorb_sigmaz2 + adiabatic_coef*Pav*(
							(1 / absorb_betay1 - 1)*absorb_sigmaz2*(x_over_Y*(Qvy1[i + 1][j][k] - Qvy1[i][j][k])*diff_x
							+ y_over_Y*(Qvy1[i][j + 1][k] - Qvy1[i][j][k])*diff_y
							+ z_over_Y*(Qvy1[i][j][k + 1] - Qvy1[i][j][k])*diff_z)
							+ (1 / absorb_betaz2 - 1)*absorb_sigmay1*(x_over_Z*(Qwy1[i + 1][j][k] - Qwy1[i][j][k])*diff_x
							+ y_over_Z*(Qwy1[i][j + 1][k] - Qwy1[i][j][k])*diff_y
							+ z_over_Z*(Qwy1[i][j][k + 1] - Qwy1[i][j][k])*diff_z)))*diff_t;
					}
					else if (jn >= JMAX - JMAX2){
						jj = j - (mpi_JMAX - mpi_JMAX2);
						absorb_sigmay2 = PML_AbsorbYmax*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);
						absorb_betay2 = 1.0 + 25.0*pow((Y[j] - PML_y2)*PML_Length2, PML_alphay);

						p_nn[i][j][k] += -(Qpy2[i][jj][k] * absorb_sigmay2*absorb_sigmaz2 + adiabatic_coef*Pav*(
							(1 / absorb_betay2 - 1)*absorb_sigmaz2*(x_over_Y*(Qvy2[i + 1][jj][k] - Qvy2[i][jj][k])*diff_x
							+ y_over_Y*(Qvy2[i][jj + 1][k] - Qvy2[i][jj][k])*diff_y
							+ z_over_Y*(Qvy2[i][jj][k + 1] - Qvy2[i][jj][k])*diff_z)
							+ (1 / absorb_betaz2 - 1)*absorb_sigmay2*(x_over_Z*(Qwy2[i + 1][jj][k] - Qwy2[i][jj][k])*diff_x
							+ y_over_Z*(Qwy2[i][jj + 1][k] - Qwy2[i][jj][k])*diff_y
							+ z_over_Z*(Qwy2[i][jj][k + 1] - Qwy2[i][jj][k])*diff_z)))*diff_t;
					}/**/
				}

				// calculate inside hill
				if (cr_judge != 0 && (in >= istar && in <= istop) && (jn >= jstar && jn <= jstop)
					&& (kn >= kstar && kn <= kstop) && cr_uc[in - istar][jn - jstar][kn - kstar] ==1 ){
					p_nn[i][j][k] = -crpo_Gama1*(x_over_X*(temp_u[i + 1][j][k] - temp_u[i][j][k])*diff_x
						+ y_over_X*(temp_u[i][j + 1][k] - temp_u[i][j][k])*diff_y
						+ z_over_X*(temp_u[i][j][k + 1] - temp_u[i][j][k])*diff_z
						+ x_over_Y*(temp_v[i + 1][j][k] - temp_v[i][j][k])*diff_x
						+ y_over_Y*(temp_v[i][j + 1][k] - temp_v[i][j][k])*diff_y
						+ z_over_Y*(temp_v[i][j][k + 1] - temp_v[i][j][k])*diff_z
						+ x_over_Z*(temp_w[i + 1][j][k] - temp_w[i][j][k])*diff_x
						+ y_over_Z*(temp_w[i][j + 1][k] - temp_w[i][j][k])*diff_y
						+ z_over_Z*(temp_w[i][j][k + 1] - temp_w[i][j][k])*diff_z);
				}
				if (cr_judge != 0 && (in >= istar && in <= istop) && (jn >= jstar && jn <= jstop)
					&& (kn >= kstar && kn <= kstop) && cr_uc[in - istar][jn - jstar][kn - kstar] ==2){
					p_nn[i][j][k] = -crpo_Gama2*(x_over_X*(temp_u[i + 1][j][k] - temp_u[i][j][k])*diff_x
						+ y_over_X*(temp_u[i][j + 1][k] - temp_u[i][j][k])*diff_y
						+ z_over_X*(temp_u[i][j][k + 1] - temp_u[i][j][k])*diff_z
						+ x_over_Y*(temp_v[i + 1][j][k] - temp_v[i][j][k])*diff_x
						+ y_over_Y*(temp_v[i][j + 1][k] - temp_v[i][j][k])*diff_y
						+ z_over_Y*(temp_v[i][j][k + 1] - temp_v[i][j][k])*diff_z
						+ x_over_Z*(temp_w[i + 1][j][k] - temp_w[i][j][k])*diff_x
						+ y_over_Z*(temp_w[i][j + 1][k] - temp_w[i][j][k])*diff_y
						+ z_over_Z*(temp_w[i][j][k + 1] - temp_w[i][j][k])*diff_z);
				}
			}
		}
	}
}
void air::Update_PML_Qvw()
{
	int i, j, k, ii, jj, kk;
	for (i = 0; i < mpi_IMAX1; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				Qux1[i][j][k] += diff_t*u_nn[i][j][k];
				Qvx1[i][j][k] += diff_t*v_nn[i][j][k];
				Qwx1[i][j][k] += diff_t*w_nn[i][j][k];
			}
		}
	}
	for (i = mpi_IMAX - mpi_IMAX2; i < mpi_IMAX; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				ii = i - (mpi_IMAX - mpi_IMAX2);
				Qux2[ii][j][k] += diff_t*u_nn[i][j][k];
				Qvx2[ii][j][k] += diff_t*v_nn[i][j][k];
				Qwx2[ii][j][k] += diff_t*w_nn[i][j][k];
			}
		}
	}
	if (PML_JMAX1 >= PML_JMAX1_time){
		for (j = 0; j < mpi_JMAX1; j++){
			for (k = 0; k < mpi_KMAX; k++){
				for (i = 0; i < mpi_IMAX; i++){
					Quy1[i][j][k] += diff_t*u_nn[i][j][k];
					Qvy1[i][j][k] += diff_t*v_nn[i][j][k];
					Qwy1[i][j][k] += diff_t*w_nn[i][j][k];
				}
			}
		}
	}
	for (j = mpi_JMAX - mpi_JMAX2; j < mpi_JMAX; j++){
		for (k = 0; k < mpi_KMAX; k++){
			for (i = 0; i < mpi_IMAX; i++){
				jj = j - (mpi_JMAX - mpi_JMAX2);
				Quy2[i][jj][k] += diff_t*u_nn[i][j][k];
				Qvy2[i][jj][k] += diff_t*v_nn[i][j][k];
				Qwy2[i][jj][k] += diff_t*w_nn[i][j][k];
			}
		}
	}

	for (k = 0; k < mpi_KMAX1; k++){
		for (i = 0; i < mpi_IMAX; i++){
			for (j = 0; j < mpi_JMAX; j++){
				Quz1[i][j][k] += diff_t*u_nn[i][j][k];
				Qvz1[i][j][k] += diff_t*v_nn[i][j][k];
				Qwz1[i][j][k] += diff_t*w_nn[i][j][k];
			}
		}
	}
	for (k = mpi_KMAX - mpi_KMAX2; k < mpi_KMAX; k++){
		for (i = 0; i < mpi_IMAX; i++){
			for (j = 0; j < mpi_JMAX; j++){
				kk = k - (mpi_KMAX - mpi_KMAX2);
				Quz2[i][j][kk] += diff_t*u_nn[i][j][k];
				Qvz2[i][j][kk] += diff_t*v_nn[i][j][k];
				Qwz2[i][j][kk] += diff_t*w_nn[i][j][k];
			}
		}
	}
}

void air::Update_PML_Qp()
{
	int i, j, k, ii, jj, kk;
	for (i = 0; i < mpi_IMAX1; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				Qpx1[i][j][k] += diff_t*p_nn[i][j][k];
			}
		}
	}
	for (i = mpi_IMAX - mpi_IMAX2; i < mpi_IMAX; i++){
		for (j = 0; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				ii = i - (mpi_IMAX - mpi_IMAX2);
				Qpx2[ii][j][k] += diff_t*p_nn[i][j][k];
			}
		}
	}
	if (PML_JMAX1 >= PML_JMAX1_time){
		for (j = 0; j < mpi_JMAX1; j++){
			for (k = 0; k < mpi_KMAX; k++){
				for (i = 0; i < mpi_IMAX; i++){
					Qpy1[i][j][k] += diff_t*p_nn[i][j][k];
				}
			}
		}
	}
	for (j = mpi_JMAX - mpi_JMAX2; j < mpi_JMAX; j++){
		for (k = 0; k < mpi_KMAX; k++){
			for (i = 0; i < mpi_IMAX; i++){
				jj = j - (mpi_JMAX - mpi_JMAX2);
				Qpy2[i][jj][k] += diff_t*p_nn[i][j][k];
			}
		}
	}

	for (k = 0; k < mpi_KMAX1; k++){
		for (i = 0; i < mpi_IMAX; i++){
			for (j = 0; j < mpi_JMAX; j++){
				Qpz1[i][j][k] += diff_t*p_nn[i][j][k];
			}
		}
	}
	for (k = mpi_KMAX - mpi_KMAX2; k < mpi_KMAX; k++){
		for (i = 0; i < mpi_IMAX; i++){
			for (j = 0; j < mpi_JMAX; j++){
				kk = k - (mpi_KMAX - mpi_KMAX2);
				Qpz2[i][j][kk] += diff_t*p_nn[i][j][k];
			}
		}
	}
}
void air::Set_Timer(int N)
{
	PML_JMAX1 = N;
	//if (mpi_rank == 0) 
	cout << "rank = "<< mpi_rank << " data loaded, PML_JMAX1 = " << N << endl;
}

void air::SetQMove(int N)
{
	int i, j, k, in, jn, kn, ii, jj, kk;
	int temp_J;
	int NN;

	if (N == 0){
		temp_J = 0;
		for (i = 0; i < mpi_IMAX1; i++){
			for (j = 0; j < mpi_JMAX; j++){
				for (k = 0; k < mpi_KMAX; k++){
					Qux1[i][j][k] = 0.0;
					Qvx1[i][j][k] = 0.0;
					Qwx1[i][j][k] = 0.0;
					Qpx1[i][j][k] = 0.0;
				}
			}
		}
		for (i = mpi_IMAX - mpi_IMAX2; i < mpi_IMAX; i++){
			for (j = 0; j < mpi_JMAX; j++){
				for (k = 0; k < mpi_KMAX; k++){
					ii = i - (mpi_IMAX - mpi_IMAX2);
					Qux2[ii][j][k] = 0.0;
					Qvx2[ii][j][k] = 0.0;
					Qwx2[ii][j][k] = 0.0;
					Qpx2[ii][j][k] = 0.0;
				}
			}
		}
		for (j = 0; j < mpi_JMAX1; j++){
			for (k = 0; k < mpi_KMAX; k++){
				for (i = 0; i < mpi_IMAX; i++){
					Quy1[i][j][k] = 0.0;
					Qvy1[i][j][k] = 0.0;
					Qwy1[i][j][k] = 0.0;
					Qpy1[i][j][k] = 0.0;
				}
			}
		}
		for (j = mpi_JMAX - mpi_JMAX2; j < mpi_JMAX; j++){
			for (k = 0; k < mpi_KMAX; k++){
				for (i = 0; i < mpi_IMAX; i++){
					jj = j - (mpi_JMAX - mpi_JMAX2);
					Quy2[i][jj][k] = 0.0;
					Qvy2[i][jj][k] = 0.0;
					Qwy2[i][jj][k] = 0.0;
					Qpy2[i][jj][k] = 0.0;
				}
			}
		}

		for (k = 0; k < mpi_KMAX1; k++){
			for (i = 0; i < mpi_IMAX; i++){
				for (j = 0; j < mpi_JMAX; j++){
					Quz1[i][j][k] = 0.0;
					Qvz1[i][j][k] = 0.0;
					Qwz1[i][j][k] = 0.0;
					Qpz1[i][j][k] = 0.0;
				}
			}
		}
		for (k = mpi_KMAX - mpi_KMAX2; k < mpi_KMAX; k++){
			for (i = 0; i < mpi_IMAX; i++){
				for (j = 0; j < mpi_JMAX; j++){
					kk = k - (mpi_KMAX - mpi_KMAX2);
					Quz2[i][j][kk] = 0.0;
					Qvz2[i][j][kk] = 0.0;
					Qwz2[i][j][kk] = 0.0;
					Qpz2[i][j][kk] = 0.0;
				}
			}
		}
	}
	else temp_J = move_frame.JMAX - move_frame.lead_DJ;

	if (temp_J != 0){
		if (mpi_IMAX1 != 0){
			NN = 0;
			for (i = 0; i < mpi_IMAX1; i++){
				for (j = 1; j < move_frame.lead_DJ + 1; j++){
					for (k = 0; k < mpi_KMAX; k++){
						uws[NN] = Qux1[i][j][k];
						vws[NN] = Qvx1[i][j][k];
						wws[NN] = Qwx1[i][j][k];
						pws[NN] = Qpx1[i][j][k];
						NN = NN + 1;
					}
				}
			}
			MPI_Sendrecv(uws, NN, MPI_DOUBLE, mpi_nleft, 0,
				uer, NN, MPI_DOUBLE, mpi_nright, 0,mpi_comm3d, &status);
			MPI_Sendrecv(uws, NN, MPI_DOUBLE, mpi_nleft, 0,
				uer, NN, MPI_DOUBLE, mpi_nright, 0, mpi_comm3d, &status);
			MPI_Sendrecv(vws, NN, MPI_DOUBLE, mpi_nleft, 1,
				ver, NN, MPI_DOUBLE, mpi_nright, 1, mpi_comm3d, &status);
			MPI_Sendrecv(wws, NN, MPI_DOUBLE, mpi_nleft, 2,
				wer, NN, MPI_DOUBLE, mpi_nright, 2, mpi_comm3d, &status);
			MPI_Sendrecv(pws, NN, MPI_DOUBLE, mpi_nleft, 3,
				per, NN, MPI_DOUBLE, mpi_nright, 3, mpi_comm3d, &status);

			NN = 0;
			for (i = 0; i < mpi_IMAX1; i++){
				for (j = 0; j < mpi_JMAX; j++){
					for (k = 0; k < mpi_KMAX; k++){
						jn = mpi_j1 + j;
						if (jn < temp_J + 1 && j + move_frame.lead_DJ - 1 < mpi_JMAX - 1){
							Qux1[i][j][k] = Qux1[i][j + move_frame.trail_DJ - 1][k];
							Qvx1[i][j][k] = Qvx1[i][j + move_frame.trail_DJ - 1][k];
							Qwx1[i][j][k] = Qwx1[i][j + move_frame.trail_DJ - 1][k];
							Qpx1[i][j][k] = Qpx1[i][j + move_frame.trail_DJ - 1][k];
						}
						else if (jn >= temp_J + 1){
							Qux1[i][j][k] = 0.0;
							Qvx1[i][j][k] = 0.0;
							Qwx1[i][j][k] = 0.0;
							Qpx1[i][j][k] = 0.0;
						}
						if (j >= mpi_JMAX - move_frame.lead_DJ && mpi_coords[1] != mpi_yarea - 1){
							Qux1[i][j][k] = uer[NN];
							Qvx1[i][j][k] = ver[NN];
							Qwx1[i][j][k] = wer[NN];
							Qpx1[i][j][k] = per[NN];
							NN = NN + 1;
						}
					}
				}
			}
		}
		if (mpi_IMAX2 != 0){
			NN = 0;
			for (i = mpi_IMAX - mpi_IMAX2; i < mpi_IMAX; i++){
				for (j = 1; j < move_frame.lead_DJ + 1; j++){
					for (k = 0; k < mpi_KMAX; k++){
						ii = i - (mpi_IMAX - mpi_IMAX2);
						uws[NN] = Qux2[ii][j][k];
						vws[NN] = Qvx2[ii][j][k];
						wws[NN] = Qwx2[ii][j][k];
						pws[NN] = Qpx2[ii][j][k];
						NN = NN + 1;
					}
				}
			}
			MPI_Sendrecv(uws, NN, MPI_DOUBLE, mpi_nleft, 4,
				uer, NN, MPI_DOUBLE, mpi_nright, 4, mpi_comm3d, &status);
			MPI_Sendrecv(vws, NN, MPI_DOUBLE, mpi_nleft, 5,
				ver, NN, MPI_DOUBLE, mpi_nright, 5, mpi_comm3d, &status);
			MPI_Sendrecv(wws, NN, MPI_DOUBLE, mpi_nleft, 6,
				wer, NN, MPI_DOUBLE, mpi_nright, 6, mpi_comm3d, &status);
			MPI_Sendrecv(pws, NN, MPI_DOUBLE, mpi_nleft, 7,
				per, NN, MPI_DOUBLE, mpi_nright, 7, mpi_comm3d, &status);

			NN = 0;
			for (i = mpi_IMAX - mpi_IMAX2; i < mpi_IMAX; i++){
				for (j = 0; j < mpi_JMAX; j++){
					for (k = 0; k < mpi_KMAX; k++){
						jn = mpi_j1 + j;
						ii = i - (mpi_IMAX - mpi_IMAX2);
						if (jn < temp_J + 1 && j + move_frame.lead_DJ - 1 < mpi_JMAX - 1){
							Qux2[ii][j][k] = Qux2[ii][j + move_frame.trail_DJ - 1][k];
							Qvx2[ii][j][k] = Qvx2[ii][j + move_frame.trail_DJ - 1][k];
							Qwx2[ii][j][k] = Qwx2[ii][j + move_frame.trail_DJ - 1][k];
							Qpx2[ii][j][k] = Qpx2[ii][j + move_frame.trail_DJ - 1][k];
						}
						else if (jn >= temp_J + 1){
							Qux2[ii][j][k] = 0.0;
							Qvx2[ii][j][k] = 0.0;
							Qwx2[ii][j][k] = 0.0;
							Qpx2[ii][j][k] = 0.0;
						}
						if (j >= mpi_JMAX - move_frame.lead_DJ && mpi_coords[1] != mpi_yarea - 1){
							Qux2[ii][j][k] = uer[NN];
							Qvx2[ii][j][k] = ver[NN];
							Qwx2[ii][j][k] = wer[NN];
							Qpx2[ii][j][k] = per[NN];
							NN = NN + 1;
						}
					}
				}
			}
		}
		if (mpi_KMAX1 != 0){
			NN = 0;
			for (k = 0; k < mpi_KMAX1; k++){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 1; j < move_frame.lead_DJ + 1; j++){
						uws[NN] = Quz1[i][j][k];
						vws[NN] = Qvz1[i][j][k];
						wws[NN] = Qwz1[i][j][k];
						pws[NN] = Qpz1[i][j][k];
						NN = NN + 1;
					}
				}
			}

			MPI_Sendrecv(uws, NN, MPI_DOUBLE, mpi_nleft, 0,
				uer, NN, MPI_DOUBLE, mpi_nright, 0, mpi_comm3d, &status);
			MPI_Sendrecv(vws, NN, MPI_DOUBLE, mpi_nleft, 1,
				ver, NN, MPI_DOUBLE, mpi_nright, 1, mpi_comm3d, &status);
			MPI_Sendrecv(wws, NN, MPI_DOUBLE, mpi_nleft, 2,
				wer, NN, MPI_DOUBLE, mpi_nright, 2, mpi_comm3d, &status);
			MPI_Sendrecv(pws, NN, MPI_DOUBLE, mpi_nleft, 3,
				per, NN, MPI_DOUBLE, mpi_nright, 3, mpi_comm3d, &status);

			NN = 0;
			for (k = 0; k < mpi_KMAX1; k++){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX; j++){
						jn = mpi_j1 + j;
						if (jn < temp_J + 1 && j + move_frame.lead_DJ - 1 < mpi_JMAX - 1){
							Quz1[i][j][k] = Quz1[i][j + move_frame.trail_DJ - 1][k];
							Qvz1[i][j][k] = Qvz1[i][j + move_frame.trail_DJ - 1][k];
							Qwz1[i][j][k] = Qwz1[i][j + move_frame.trail_DJ - 1][k];
							Qpz1[i][j][k] = Qpz1[i][j + move_frame.trail_DJ - 1][k];
						}
						else if (jn >= temp_J + 1){
							Quz1[i][j][k] = 0.0;
							Qvz1[i][j][k] = 0.0;
							Qwz1[i][j][k] = 0.0;
							Qpz1[i][j][k] = 0.0;
						}
						if (j >= mpi_JMAX - move_frame.lead_DJ && mpi_coords[1] != mpi_yarea - 1){
							Quz1[i][j][k] = uer[NN];
							Qvz1[i][j][k] = ver[NN];
							Qwz1[i][j][k] = wer[NN];
							Qpz1[i][j][k] = per[NN];
							NN = NN + 1;
						}
					}
				}
			}
		}
		if (mpi_KMAX2 != 0){
			NN = 0;
			for (k = mpi_KMAX - mpi_KMAX2; k < mpi_KMAX; k++){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 1; j < move_frame.lead_DJ + 1; j++){
						kk = k - (mpi_KMAX - mpi_KMAX2);
						uws[NN] = Quz2[i][j][kk];
						vws[NN] = Qvz2[i][j][kk];
						wws[NN] = Qwz2[i][j][kk];
						pws[NN] = Qpz2[i][j][kk];
						NN = NN + 1;
					}
				}
			}
			MPI_Sendrecv(uws, NN, MPI_DOUBLE, mpi_nleft, 4,
				uer, NN, MPI_DOUBLE, mpi_nright, 4, mpi_comm3d, &status);
			MPI_Sendrecv(vws, NN, MPI_DOUBLE, mpi_nleft, 5,
				ver, NN, MPI_DOUBLE, mpi_nright, 5, mpi_comm3d, &status);
			MPI_Sendrecv(wws, NN, MPI_DOUBLE, mpi_nleft, 6,
				wer, NN, MPI_DOUBLE, mpi_nright, 6, mpi_comm3d, &status);
			MPI_Sendrecv(pws, NN, MPI_DOUBLE, mpi_nleft, 7,
				per, NN, MPI_DOUBLE, mpi_nright, 7, mpi_comm3d, &status);

			NN = 0;
			for (k = mpi_KMAX - mpi_KMAX2; k < mpi_KMAX; k++){
				for (i = 0; i < mpi_IMAX; i++){
					for (j = 0; j < mpi_JMAX; j++){
						jn = mpi_j1 + j;
						kk = k - (mpi_KMAX - mpi_KMAX2);
						if (jn < temp_J + 1 && j + move_frame.lead_DJ - 1 < mpi_JMAX - 1){
							Quz2[i][j][kk] = Quz2[i][j + move_frame.trail_DJ - 1][kk];
							Qvz2[i][j][kk] = Qvz2[i][j + move_frame.trail_DJ - 1][kk];
							Qwz2[i][j][kk] = Qwz2[i][j + move_frame.trail_DJ - 1][kk];
							Qpz2[i][j][kk] = Qpz2[i][j + move_frame.trail_DJ - 1][kk];
						}
						else if (jn >= temp_J + 1){
							Quz2[i][j][kk] = 0.0;
							Qvz2[i][j][kk] = 0.0;
							Qwz2[i][j][kk] = 0.0;
							Qpz2[i][j][kk] = 0.0;
						}
						if (j >= mpi_JMAX - move_frame.lead_DJ && mpi_coords[1] != mpi_yarea - 1){
							Quz2[i][j][kk] = uer[NN];
							Qvz2[i][j][kk] = ver[NN];
							Qwz2[i][j][kk] = wer[NN];
							Qpz2[i][j][kk] = per[NN];
							NN = NN + 1;
						}
					}
				}
			}
		}
	}
}