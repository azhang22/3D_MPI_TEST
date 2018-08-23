//+++++++++++++++++++++++filename: calculation_2D.cpp +++++++++++++++++++++++++++++//

//-----------------------fathe class for air and porous class------------------//
/*******************************************************************************/
/* air and porous class will inherit from the class calculation_2D.
/*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "calculation_2D.h"

calculation_2D::calculation_2D(scheme_list scheme1,DifferenceStep Step,char *coordi,
			MovingFrame MF,const double gauss_width,MPI_info mpi_var,int judge_porous)
{
	int i,j,k,N;
	char coordi1[20]="",coordi2[20]="",coordi3[20]="";
	scheme=scheme1;
	move_frame=MF;
	diff_t=Step.diff_t;
	diff_x=1.0/Step.diff_x;
	diff_y=1.0/Step.diff_y;
	diff_z=1.0/Step.diff_z;
	diff_x1=Step.diff_x;
	diff_y1=Step.diff_y;
	diff_z1=Step.diff_z;
	gaussian_coef=4.0*log(2.0)/pow(gauss_width,2);
	judge_porous1=judge_porous;
	// define each subdomain information
	mpi_coords=new int[3];
	mpi_rank=mpi_var.mpi_rank;
	mpi_size=mpi_var.mpi_size;
	mpi_xarea=mpi_var.mpi_xarea;
	mpi_yarea=mpi_var.mpi_yarea;
	mpi_zarea=mpi_var.mpi_zarea;
	mpi_nleft=mpi_var.mpi_nleft;
	mpi_nright=mpi_var.mpi_nright;
	mpi_nbottom=mpi_var.mpi_nbottom;
	mpi_ntop=mpi_var.mpi_ntop;
	mpi_nfront=mpi_var.mpi_nfront;
	mpi_nback=mpi_var.mpi_nback;
	mpi_comm3d=mpi_var.mpi_comm3d;
	for (i=0;i<3;i++)mpi_coords[i]=mpi_var.mpi_coords[i];

	strcat_s(coordi1,coordi);
	strcat_s(coordi1,"1");
	strcat_s(coordi1,".dat");
	//read data from files
	ifstream myfile1(coordi1,ios::in|ios::binary);
	if(!myfile1){
		cout<<"cannot open file:"<<coordi1<<" for read!!!"<<endl;
		return;
	}
	myfile1>>IMAX>>IMAX1>>IMAX2>>JMAX>>JMAX1>>JMAX2>>KMAX>>KMAX1>>KMAX2>>IJKMAX;
	JMAX3=JMAX1;
	//the following setting is for one moving frame++++
	// if judge of moving frame is on, let IJKMAX is the times of move_frame.JMAX;
	// and JMAX of geometry should be equal to move_frame.JMAX
	int N_limit;
	if(move_frame.Judge==1){
		if ((IJKMAX-move_frame.JMAX)%(move_frame.lead_DJ-1)!=0){
			N_limit=(IJKMAX-move_frame.JMAX)/(move_frame.lead_DJ-1);
			IJKMAX=move_frame.JMAX+(N_limit+1)*(move_frame.lead_DJ-1);
		}
		JMAX=move_frame.JMAX;
	}else{
		move_frame.JMAX=JMAX;
		IJKMAX=JMAX;
	}
	// domain decomposition
	mpe_decomp1d(IMAX,mpi_xarea,mpi_coords[0],mpi_i1,mpi_i2);
	mpe_decomp1d(JMAX,mpi_yarea,mpi_coords[1],mpi_j1,mpi_j2);
	mpe_decomp1d(KMAX,mpi_zarea,mpi_coords[2],mpi_k1,mpi_k2);
	mpi_i1=mpi_i1-1;
	mpi_i2=mpi_i2-1;
	if (mpi_coords[0]!=0) mpi_i1=mpi_i1-1;
	if (mpi_coords[0]!=mpi_xarea-1)mpi_i2=mpi_i2+1;
	mpi_IMAX=mpi_i2-mpi_i1+1;

	mpi_j1=mpi_j1-1;
	mpi_j2=mpi_j2-1;
	if (mpi_coords[1]!=0) mpi_j1=mpi_j1-1;
	if (mpi_coords[1]!=mpi_yarea-1)mpi_j2=mpi_j2+1;
	mpi_JMAX=mpi_j2-mpi_j1+1;

	if (judge_porous1==0){
		mpi_k1=mpi_k1-1;
		mpi_k2=mpi_k2-1;
		if (mpi_coords[2]!=0) mpi_k1=mpi_k1-1;
		if (mpi_coords[2]!=mpi_zarea-1)mpi_k2=mpi_k2+1;
		mpi_KMAX=mpi_k2-mpi_k1+1;
	}else{
		mpi_k1=0;
		mpi_k2=KMAX-1;
		mpi_KMAX=mpi_k2-mpi_k1+1;
	}

	// allocate memory for geometry
	X=new double [mpi_IMAX];Y=new double [mpi_JMAX];
	Z=new double [mpi_KMAX];
	X_over_x=new double[mpi_IMAX];Y_over_y=new double [mpi_JMAX];
	Z_over_z=new double[mpi_KMAX];
	X_over_y=0.0;X_over_z=0.0;
	Y_over_x=0.0;Y_over_z=0.0;
	Z_over_x=0.0;Z_over_y=0.0;

	//allocate memory for jacobi,velocity, pressure
	Jacobian=new double** [mpi_IMAX];
	u_n=new double** [mpi_IMAX];v_n=new double** [mpi_IMAX];w_n=new double** [mpi_IMAX];p_n=new double** [mpi_IMAX];
	u_nn=new double** [mpi_IMAX];v_nn=new double** [mpi_IMAX];w_nn=new double** [mpi_IMAX];p_nn=new double** [mpi_IMAX];
	for(i=0;i<mpi_IMAX;i++){
		Jacobian[i]=new double* [mpi_JMAX];
		u_n[i]=new double* [mpi_JMAX];v_n[i]=new double* [mpi_JMAX];
		w_n[i]=new double* [mpi_JMAX];p_n[i]=new double* [mpi_JMAX];
		u_nn[i]=new double* [mpi_JMAX];v_nn[i]=new double* [mpi_JMAX];
		w_nn[i]=new double* [mpi_JMAX];p_nn[i]=new double* [mpi_JMAX];
		for (j=0;j<mpi_JMAX;j++){
			Jacobian[i][j]=new double [mpi_KMAX];
			u_n[i][j]=new double [mpi_KMAX];v_n[i][j]=new double [mpi_KMAX];
			w_n[i][j]=new double [mpi_KMAX];p_n[i][j]=new double [mpi_KMAX];
			u_nn[i][j]=new double [mpi_KMAX];v_nn[i][j]=new double [mpi_KMAX];
			w_nn[i][j]=new double [mpi_KMAX];p_nn[i][j]=new double [mpi_KMAX];
		}
	}
	if(move_frame.Judge==1){
		N=move_frame.lead_DJ*mpi_IMAX*mpi_KMAX;
		uer=new double [N];
		uws=new double [N];
		ver=new double [N];
		vws=new double [N];
		wer=new double [N];
		wws=new double [N];
		per=new double [N];
		pws=new double [N];
	}

	int ii,jj,kk;
	double temp_X,temp_Y,temp_Z;
	axisbound_X=new double [IMAX];
	axisbound_Y=new double [JMAX];
	axisbound_Z=new double [KMAX];
	// Input 3D geometric data for X,Y,Z;
	for(i=0;i<IMAX;i++){
		myfile1>>temp_X;
		if (i>=mpi_i1 && i<=mpi_i2){
			ii=i-mpi_i1;
			X[ii]=temp_X;
		}
		axisbound_X[i]=temp_X;
	}
	myfile1.close();

	strcat_s(coordi2,coordi);
	strcat_s(coordi2,"2");
	strcat_s(coordi2,".dat");
	//read data from files
	ifstream myfile2(coordi2,ios::in|ios::binary);
	if(!myfile2){
		cout<<"cannot open file:"<<coordi2<<" for read!!!"<<endl;
		return;
	}
	for(j=0;j<JMAX;j++){
		myfile2>>temp_Y;
		if (j>=mpi_j1 && j<=mpi_j2 ){
			jj=j-mpi_j1;
			Y[jj]=temp_Y;
		}
		axisbound_Y[j]=temp_Y;
	}
	myfile2.close();

	strcat_s(coordi3,coordi);
	strcat_s(coordi3,"3");
	strcat_s(coordi3,".dat");
	ifstream myfile3(coordi3,ios::in|ios::binary);
	if(!myfile3){
		cout<<"cannot open file:"<<coordi3<<" for read!!!"<<endl;
		return;
	}
	// The position of orignal point in Y direction;
	for(k=0;k<KMAX;k++){
		myfile3>>temp_Z;
		if (k>=mpi_k1 && k<=mpi_k2){
			kk=k-mpi_k1;
			Z[kk]=temp_Z;
		}
		axisbound_Z[k]=temp_Z;
	}
	myfile3.close();

	// out put x,y,z,pressure for pressure contour
	out_n=4;
	if ((IMAX-1)%out_n==0)IMAX_out=(IMAX-1)/out_n+1;
	else IMAX_out=(IMAX-1)/out_n+2;
	if ((JMAX-1)%out_n==0)JMAX_out=(JMAX-1)/out_n+1;
	else JMAX_out=(JMAX-1)/out_n+2;
	if ((KMAX-1)%out_n==0)KMAX_out=(KMAX-1)/out_n+1;
	else KMAX_out=(KMAX-1)/out_n+2;
	X_out=new double [IMAX_out];
	Y_out=new double [JMAX_out];
	Z_out=new double [KMAX_out];
	p_in=new double [IMAX_out*JMAX_out*KMAX_out];
	p_out=new double [IMAX_out*JMAX_out*KMAX_out];
	
	int source_judge;
	ps = new double[100000];
	source_judge = 0;
	if (source_judge == 1) {
		ifstream mysource("source_input.dat", ios::in | ios::binary);
		if (!mysource){
			cout << "cannot open file:" << coordi << " for read!!!" << endl;
			return;
		}
		for (i = 0; i <= 100000; i++){
			mysource >> ps[i];
			//if (mpi_rank == 0) cout << "ps[" << i << "] = " << ps[i] << endl;
		}
		mysource.close();
		//cout << "Source file loaded!" << endl;
	}
}

calculation_2D::~calculation_2D()
{
	int i,j;
	for(i=0;i<mpi_IMAX;i++){
		for(j=0;j<mpi_JMAX;j++){
			delete[] p_n[i][j];delete[] w_n[i][j];delete[] v_n[i][j];delete[] u_n[i][j];
			delete[] p_nn[i][j];delete[] w_nn[i][j];delete[] v_nn[i][j];delete[] u_nn[i][j];
			delete[] Jacobian[i][j];
		}
		delete[] p_n[i];delete[] w_n[i];delete[] v_n[i];delete[] u_n[i];
		delete[] p_nn[i];delete[] w_nn[i];delete[] v_nn[i];delete[] u_nn[i];
		delete[] Jacobian[i];
	}
	delete[] p_n;delete[] w_n;delete[] v_n;delete[] u_n;
	delete[] p_nn;delete[] w_nn;delete[] v_nn;delete[] u_nn;
	delete[] X;delete[] Y;delete[] Z;delete[] Jacobian;
	delete[] X_over_x;delete[] Y_over_y;delete[] Z_over_z;
	delete[] axisbound_X;delete[] axisbound_Y;delete[] axisbound_Z;
	/*
	for (i=0;i<IMAX_out;i++){
		for (j=0;j<JMAX_out;j++){
			delete[] p_out[i][j];
			delete[] p_in[i][j];
		}
		delete[] p_out[i];
		delete[] p_in[i];
	}
	*/
	delete[] X_out;delete[] Y_out;delete[] Z_out;
	delete[] p_out;
	delete[] p_in;

	if(move_frame.Judge==1){
		delete[] uer;delete[]uws;delete[] ver;delete[]vws;
		delete[] wer;delete[]wws;delete[] per;delete[]pws;
	}
	delete[] mpi_coords;
}


double calculation_2D::cal_InitialPressure(double x,double y,double z)
{
	double r;
	r=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	//this is from Salomons paper(2002) p=exp(-40*r*r);
	return (1.0*exp(-gaussian_coef*r*r));
}

void calculation_2D::Set_InitialCond(Position source,int judge)
{
	int i,j,k;
	for( i=0;i<mpi_IMAX;i++){
		for(j=0;j<mpi_JMAX;j++){
			for(k=0;k<mpi_KMAX;k++){
				u_n[i][j][k]=0;
				v_n[i][j][k]=0;
				w_n[i][j][k]=0;
				if (judge == 1)		p_n[i][j][k]=cal_InitialPressure(X[i]-source.x,	Y[j]-source.y,Z[k]-source.z);
				else p_n[i][j][k] = 0;
			}
		}
	}
}


void calculation_2D::get_position(int &ii,int &jj,int &kk,Position receiver)
{
	int i,j,k,N,NN;
	double Y_frame;
	Y_frame=axisbound_Y[JMAX-1]-axisbound_Y[0];
	// If N>0 then it is movingframe;
	// If N=0 then it is in the first compuational domain;
	N=int(receiver.y/Y_frame+1.0e-6);
	// return the position of receiver in the physical domain;
	
	for (i=0;i<IMAX;i++){
		if (receiver.x-axisbound_X[i]<=diff_x1/8.0){
			ii=i;
			break;
		}
	}
	NN=N*(JMAX-1);
	if (NN>=IJKMAX){
		ii=0;jj=0;kk=0;
		return;
	}else{
		for (j=NN;j<IJKMAX;j++){
			if (receiver.y-N*Y_frame-axisbound_Y[j-NN]<=diff_y1/8.0){
				jj=j;
				break;
			}
		}
	}
	for (k=0;k<KMAX;k++){
		if (receiver.z-axisbound_Z[k]<=diff_z1/8.0){
			kk=k;
			return;
		}
	}
	ii=0;jj=0;kk=0;
}

void calculation_2D::UpdateInitialCond(int N)
{
	//transfer data under time n+1 to data under time n
	//update u,v,w,p in the computational domain
	int i,j,k;
	for(i=0;i<mpi_IMAX;i++){
		for(j=0;j<mpi_JMAX;j++){
			for (k=0;k<mpi_KMAX;k++){
				u_n[i][j][k]=u_nn[i][j][k];
				v_n[i][j][k]=v_nn[i][j][k];
				w_n[i][j][k]=w_nn[i][j][k];
				p_n[i][j][k]=p_nn[i][j][k];
			}
		}
	}
}

void calculation_2D::UpdateSource(double f0, int N, Position source)
{
	int i, j, k;
	double mvspd = 0; // object moving speed
	//double omega = 2 * acos(-1)*f0;
	for (i = 0; i<mpi_IMAX; i++){
		for (j = 0; j<mpi_JMAX; j++){
			for (k = 0; k<mpi_KMAX; k++){
				//p_nn[i][j][k] = p_nn[i][j][k] + sin(omega*N*diff_t)*(1 / pow(sqrt(acos(-1)*0.36), 3))*cal_InitialPressure(X[i] - source.x, Y[j] - (source.y + mvspd*N*diff_t), Z[k] - source.z);
				if (N%23 == 0) 	p_nn[i][j][k] = p_nn[i][j][k] + ps[N/23] * cal_InitialPressure(X[i] - source.x, Y[j] - (source.y + mvspd*N*diff_t), Z[k] - source.z);
			}
		}
	}
	if (mpi_rank == 0) cout << "Source strength " << ps[N] << endl;
}

double calculation_2D::get_pressure(int ii,int jj,int kk,int N)
{
	int jl,in,jn,kn;
	double temp_pvalue,temp_pvalue1;
	// Get the receiver's pressure from the computational domain;
	if ((jj>=N*(move_frame.trail_DJ-1))&&(jj<JMAX+N*(move_frame.trail_DJ-1))){
		jl=jj-N*(move_frame.trail_DJ-1);
		if (jl>mpi_j1 && jl<mpi_j2 && ii>mpi_i1 && ii<mpi_i2
									&& kk>mpi_k1 && kk<mpi_k2){
			in=ii-mpi_i1;jn=jl-mpi_j1;kn=kk-mpi_k1;
			temp_pvalue=p_n[in][jn][kn];
		}else temp_pvalue=0.0;
	}else temp_pvalue=0.0;
	//mpi_comm3d.Allreduce(&temp_pvalue,&temp_pvalue1,1,MPI_DOUBLE,MPI::SUM);
	MPI_Allreduce(&temp_pvalue, &temp_pvalue1, 1, MPI_DOUBLE, MPI_SUM, mpi_comm3d);
	return temp_pvalue1;
}

void calculation_2D::SetMovingFrame(int N)
//N is total frame number starting from 0
{
	int i,j,k,ii,jj,kk,temp_J;
	//transfer coordinate
	if (N!=0){
		for (j=0;j<mpi_JMAX;j++){
			jj=mpi_j1+j;
			Y[j]=axisbound_Y[jj]+N*(move_frame.trail_DJ-1)*diff_y1;
		}
	}
	//calculate partial differential of coordinate transformation using central difference
	for(i=1;i<mpi_IMAX-1;i++)X_over_x[i]=(X[i+1]-X[i-1])*diff_x*0.5;
	for(j=1;j<mpi_JMAX-1;j++)Y_over_y[j]=(Y[j+1]-Y[j-1])*diff_y*0.5;
	for(k=1;k<mpi_KMAX-1;k++)Z_over_z[k]=(Z[k+1]-Z[k-1])*diff_z*0.5;

	for(i=1;i<mpi_IMAX-1;i++){
		for(j=1;j<mpi_JMAX-1;j++){
			for(k=1;k<mpi_KMAX-1;k++){
				Jacobian[i][j][k]=X_over_x[i]*Y_over_y[j]*Z_over_z[k]+
									Z_over_x*X_over_y*Y_over_z+
									Y_over_x*Z_over_y*X_over_z;
				Jacobian[i][j][k]+=-X_over_x[i]*Z_over_y*Y_over_z-
										Y_over_x*X_over_y*Z_over_z[k]-
										Z_over_x*Y_over_y[j]*X_over_z;
				Jacobian[i][j][k]=1.0/Jacobian[i][j][k];
			}
		}
	}
	int NN;
	// get the data of westboundary condition for movingframe correspoing
	if (N!=0){
		NN=0;
		for (i=0;i<mpi_IMAX;i++){
			for (j=1;j<move_frame.lead_DJ+1;j++){
				for (k=0;k<mpi_KMAX;k++){
					uws[NN]=u_n[i][j][k];
					vws[NN]=v_n[i][j][k];
					wws[NN]=w_n[i][j][k];
					pws[NN]=p_n[i][j][k];
					NN=NN+1;
				}
			}
		}
		MPI_Sendrecv(uws,NN,MPI_DOUBLE,mpi_nleft,0,
							uer,NN,MPI_DOUBLE,mpi_nright,0,mpi_comm3d, &status);
		MPI_Sendrecv(vws,NN,MPI_DOUBLE,mpi_nleft,1,
							ver,NN,MPI_DOUBLE,mpi_nright,1,mpi_comm3d, &status);
		MPI_Sendrecv(wws,NN,MPI_DOUBLE,mpi_nleft,2,
							wer,NN,MPI_DOUBLE,mpi_nright,2,mpi_comm3d, &status);
		MPI_Sendrecv(pws,NN,MPI_DOUBLE,mpi_nleft,3,
							per,NN,MPI_DOUBLE,mpi_nright,3,mpi_comm3d, &status);
	}
	if (N==0) temp_J=0;
   	else temp_J=move_frame.JMAX-move_frame.lead_DJ;
	//update u v,w,p for the moving frame
	if (temp_J!=0){
		for (k=0;k<mpi_KMAX;k++){
			for (i=0;i<mpi_IMAX;i++){
				for (j=0;j<mpi_JMAX;j++){
					jj=mpi_j1+j;
					if (jj<temp_J+1 && j+move_frame.lead_DJ-1<mpi_JMAX-1){
						u_n[i][j][k]=u_n[i][j+move_frame.trail_DJ-1][k];
						v_n[i][j][k]=v_n[i][j+move_frame.trail_DJ-1][k];
						w_n[i][j][k]=w_n[i][j+move_frame.trail_DJ-1][k];
						p_n[i][j][k]=p_n[i][j+move_frame.trail_DJ-1][k];
					}else if (jj>=temp_J+1){
						u_n[i][j][k]=0.0;
						v_n[i][j][k]=0.0;
						w_n[i][j][k]=0.0;
						p_n[i][j][k]=0.0;
					}
				}
				if (mpi_coords[1]==0)p_n[i][0][k]=p_n[i][1][k];
			}
		}
	}
	if (N!=0 && mpi_coords[1]!=mpi_yarea-1){
		NN=0;
		for (i=0;i<mpi_IMAX;i++){
			for (j=mpi_JMAX-move_frame.lead_DJ;j<mpi_JMAX;j++){
				for (k=0;k<mpi_KMAX;k++){
					u_n[i][j][k]=uer[NN];
					v_n[i][j][k]=ver[NN];
					w_n[i][j][k]=wer[NN];
					p_n[i][j][k]=per[NN];
					NN=NN+1;
				}
			}
		}
	}
}
void calculation_2D::save_restart_cal(char *restartfile)
{
	int i,j,k;
	ofstream outfile11(restartfile,ios::app|ios::binary);
	outfile11.setf(ios::scientific,ios::floatfield);
	outfile11.setf(ios::showpos);
	outfile11.precision(10);
	outfile11.width(18);
	//for(i=0;i<mpi_IMAX;i++) outfile11<<X[i]<<" ";
	//for(j=0;j<mpi_JMAX;j++) outfile11<<Y[j]<<" ";
	//for(k=0;k<mpi_KMAX;k++) outfile11<<Z[k]<<" ";
	for(i=0;i<mpi_IMAX;i++){
		for(j=0;j<mpi_JMAX;j++){
			for(k=0;k<mpi_KMAX;k++){
				outfile11<<u_n[i][j][k]<<" ";
				outfile11<<v_n[i][j][k]<<" ";
				outfile11<<w_n[i][j][k]<<" ";
				outfile11<<p_n[i][j][k]<<" ";
			}
		}
	}
	outfile11<<endl;
	outfile11.close();
}

void calculation_2D::input_restart_cal(char *restart_infile,int *point_pois)
{
	int i,j,k,int_pois;
	int_pois=*point_pois;
	ifstream infile(restart_infile,ios::in|ios::binary);
	infile.seekg(int_pois);
	//for(i=0;i<mpi_IMAX;i++)infile>>X[i];
	//for(j=0;j<mpi_JMAX;j++)infile>>Y[j];
	//for(k=0;k<mpi_KMAX;k++)infile>>Z[k];
	for (i=0;i<mpi_IMAX;i++){
		for (j=0;j<mpi_JMAX;j++){
			for (k=0;k<mpi_KMAX;k++){
				infile>>u_n[i][j][k];
				infile>>v_n[i][j][k];
				infile>>w_n[i][j][k];
				infile>>p_n[i][j][k];
			}
		}
	}
	infile.ignore(100,'\n');
	int_pois=infile.tellg();
	*point_pois=int_pois;
	infile.close();
}

//cal_fuvw() and cal_fp() to get v_nn=fuvw*diff_t;p_nn=fp*diff_t//
void calculation_2D::Cal_velocity(void)
{
	int i,j,k;
	switch(scheme){
	case FB_v_p://forward - backward scheme (1st)
		cal_fuvw(u_n,v_n,w_n,p_n);
		for(i=1;i<mpi_IMAX-1;i++){
			for(j=1;j<mpi_JMAX-1;j++){
				for(k=1;k<mpi_KMAX-1;k++){
					u_nn[i][j][k]=u_n[i][j][k]+u_nn[i][j][k];
					v_nn[i][j][k]=v_n[i][j][k]+v_nn[i][j][k];
					w_nn[i][j][k]=w_n[i][j][k]+w_nn[i][j][k];
				}
			}
		}
		break;
	case FB_p_v://forward - backward scheme (1st)
		cal_fuvw(u_n,v_n,w_n,p_nn);
		for(i=1;i<mpi_IMAX-1;i++){
			for(j=1;j<mpi_JMAX-1;j++){
				for(k=1;k<mpi_KMAX-1;k++){
					u_nn[i][j][k]=u_n[i][j][k]+u_nn[i][j][k];
					v_nn[i][j][k]=v_n[i][j][k]+v_nn[i][j][k];
					w_nn[i][j][k]=w_n[i][j][k]+w_nn[i][j][k];
				}
			}
		}
		break;
	default:
		printf("wrong with scheme name");
	}//end switch
}

void calculation_2D::Cal_pressure(void)
{
	int i,j,k;
	//double temp_p;
	switch(scheme){
	case FB_v_p://forward - backward scheme (1st)
		{
		cal_fp(u_nn,v_nn,w_nn,p_n);
		for(i=1;i<mpi_IMAX-1;i++){
			for(j=1;j<mpi_JMAX-1;j++){
				for(k=1;k<mpi_KMAX-1;k++){
					p_nn[i][j][k]=p_n[i][j][k]+p_nn[i][j][k];
				}
			}
		}
		}
		break;
	case FB_p_v://forward - backward scheme (1st)
		{
		cal_fp(u_n,v_n,w_n,p_n);
		for(i=1;i<mpi_IMAX-1;i++){
			for(j=1;j<mpi_JMAX-1;j++){
				for(k=1;k<mpi_KMAX-1;k++){
					p_nn[i][j][k]=p_n[i][j][k]+p_nn[i][j][k];
				}
			}
		}
		}
		break;
	default:
		printf("wrong with scheme name");
	}//end switch
}

void calculation_2D:: MPIexchange_pressure(){
	int i,j,k;
	for(i=0;i<mpi_IMAX;i++){
		for(k=0;k<mpi_KMAX;k++){
			MPI_Sendrecv(&p_nn[i][1][k],1,MPI_DOUBLE,mpi_nleft,0,
										&p_nn[i][mpi_JMAX-1][k],1,MPI_DOUBLE,mpi_nright,0,mpi_comm3d, &status);
			MPI_Sendrecv(&p_nn[i][mpi_JMAX-2][k],1,MPI_DOUBLE,mpi_nright,1,
										&p_nn[i][0][k],1,MPI_DOUBLE,mpi_nleft,1,mpi_comm3d, &status);
										
		}
	}
	for(j=0;j<mpi_JMAX;j++){
		for(k=0;k<mpi_KMAX;k++){
			MPI_Sendrecv(&p_nn[1][j][k],1,MPI_DOUBLE,mpi_nback,2,
										&p_nn[mpi_IMAX-1][j][k],1,MPI_DOUBLE,mpi_nfront,2,mpi_comm3d, &status);
			MPI_Sendrecv(&p_nn[mpi_IMAX-2][j][k],1,MPI_DOUBLE,mpi_nfront,3,
										&p_nn[0][j][k],1,MPI_DOUBLE,mpi_nback,3,mpi_comm3d, &status);
		}
	}
	if (judge_porous1==0){
		for(j=0;j<mpi_JMAX;j++){
			for(i=0;i<mpi_IMAX;i++){
				MPI_Sendrecv(&p_nn[i][j][1],1,MPI_DOUBLE,mpi_nbottom,4,
										&p_nn[i][j][mpi_KMAX-1],1,MPI_DOUBLE,mpi_ntop,4,mpi_comm3d, &status);
				MPI_Sendrecv(&p_nn[i][j][mpi_KMAX-2],1,MPI_DOUBLE,mpi_ntop,5,
										&p_nn[i][j][0],1,MPI_DOUBLE,mpi_nbottom,5,mpi_comm3d, &status);
			}
		}
	}
}
//-----------------------------update boundary condtions of pressure--------------------//
void calculation_2D::UpdateBC_pressure(boundary_location BC)
{
	int i,j,k;
	switch(BC){
	case WestBC:
		{
			if (mpi_coords[1]==0){
				for (i=0;i<mpi_IMAX;i++){
					for(k=0;k<mpi_KMAX;k++){
						p_nn[i][0][k]=p_nn[i][1][k];//left
					}
				}
			}
		}
		break;
	case EastBC:
		{
			if (mpi_coords[1]==mpi_yarea-1){
				for (i=0;i<mpi_IMAX;i++){
					for(k=0;k<mpi_KMAX;k++){
						p_nn[i][mpi_JMAX-1][k]=p_nn[i][mpi_JMAX-2][k];//right
					}
				}
			}
		}
		break;
	case SouthBC:
		{
			if (mpi_coords[2]==0){
				for (i=0;i<mpi_IMAX;i++){
					for(j=0;j<mpi_JMAX;j++){
						p_nn[i][j][0]=p_nn[i][j][1]; // Bottom boundary
					}
				}
			}
		}
		break;
	case NorthBC:
		{
			if (mpi_coords[2]==mpi_zarea-1){
				for (i=0;i<mpi_IMAX;i++){
					for(j=0;j<mpi_JMAX;j++){
						p_nn[i][j][mpi_KMAX-1]=p_nn[i][j][mpi_KMAX-2]; // Top boundary
					}
				}
			}
		}
		break;
	case FrontBC:
		{
			if (mpi_coords[0]==0){
				for (j=0;j<mpi_JMAX;j++){
					for(k=0;k<mpi_KMAX;k++){
						p_nn[0][j][k]=p_nn[1][j][k];
					}
				}
			}
		}
		break;
	case BackBC:
		{
			if (mpi_coords[0]==mpi_xarea-1){
				for (j=0;j<mpi_JMAX;j++){
					for(k=0;k<mpi_KMAX;k++){
						p_nn[mpi_IMAX-1][j][k]=p_nn[mpi_IMAX-2][j][k];
					}
				}
			}
		}
		break;
	default:
		cout << "it is wrong to entrance boundary, Pressure BC = " << BC << endl;
	}
}
void calculation_2D:: MPIexchange_velocity()
{
	int i,j,k;
	for(i=0;i<mpi_IMAX;i++){
		for(k=0;k<mpi_KMAX;k++){
			MPI_Sendrecv(&u_nn[i][1][k],1,MPI_DOUBLE,mpi_nleft,0,
										&u_nn[i][mpi_JMAX-1][k],1,MPI_DOUBLE,mpi_nright,0,mpi_comm3d, &status);
			MPI_Sendrecv(&u_nn[i][mpi_JMAX-2][k],1,MPI_DOUBLE,mpi_nright,1,
										&u_nn[i][0][k],1,MPI_DOUBLE,mpi_nleft,1,mpi_comm3d, &status);
			MPI_Sendrecv(&v_nn[i][1][k],1,MPI_DOUBLE,mpi_nleft,2,
										&v_nn[i][mpi_JMAX-1][k],1,MPI_DOUBLE,mpi_nright,2,mpi_comm3d, &status);
			MPI_Sendrecv(&v_nn[i][mpi_JMAX-2][k],1,MPI_DOUBLE,mpi_nright,3,
										&v_nn[i][0][k],1,MPI_DOUBLE,mpi_nleft,3,mpi_comm3d, &status);
			MPI_Sendrecv(&w_nn[i][1][k],1,MPI_DOUBLE,mpi_nleft,4,
										&w_nn[i][mpi_JMAX-1][k],1,MPI_DOUBLE,mpi_nright,4,mpi_comm3d, &status);
			MPI_Sendrecv(&w_nn[i][mpi_JMAX-2][k],1,MPI_DOUBLE,mpi_nright,5,
										&w_nn[i][0][k],1,MPI_DOUBLE,mpi_nleft,5,mpi_comm3d, &status);
		}
	}
	for(j=0;j<mpi_JMAX;j++){
		for(k=0;k<mpi_KMAX;k++){
			MPI_Sendrecv(&u_nn[1][j][k],1,MPI_DOUBLE,mpi_nback,6,
										&u_nn[mpi_IMAX-1][j][k],1,MPI_DOUBLE,mpi_nfront,6,mpi_comm3d, &status);
			MPI_Sendrecv(&u_nn[mpi_IMAX-2][j][k],1,MPI_DOUBLE,mpi_nfront,7,
										&u_nn[0][j][k],1,MPI_DOUBLE,mpi_nback,7,mpi_comm3d, &status);
			MPI_Sendrecv(&v_nn[1][j][k],1,MPI_DOUBLE,mpi_nback,8,
										&v_nn[mpi_IMAX-1][j][k],1,MPI_DOUBLE,mpi_nfront,8,mpi_comm3d, &status);
			MPI_Sendrecv(&v_nn[mpi_IMAX-2][j][k],1,MPI_DOUBLE,mpi_nfront,9,
										&v_nn[0][j][k],1,MPI_DOUBLE,mpi_nback,9,mpi_comm3d, &status);
			MPI_Sendrecv(&w_nn[1][j][k],1,MPI_DOUBLE,mpi_nback,10,
										&w_nn[mpi_IMAX-1][j][k],1,MPI_DOUBLE,mpi_nfront,10,mpi_comm3d, &status);
			MPI_Sendrecv(&w_nn[mpi_IMAX-2][j][k],1,MPI_DOUBLE,mpi_nfront,11,
										&w_nn[0][j][k],1,MPI_DOUBLE,mpi_nback,11,mpi_comm3d, &status);
		}
	}
	if (judge_porous1==0){
		for(j=0;j<mpi_JMAX;j++){
			for(i=0;i<mpi_IMAX;i++){
				MPI_Sendrecv(&u_nn[i][j][1],1,MPI_DOUBLE,mpi_nbottom,12,
										&u_nn[i][j][mpi_KMAX-1],1,MPI_DOUBLE,mpi_ntop,12,mpi_comm3d, &status);
				MPI_Sendrecv(&u_nn[i][j][mpi_KMAX-2],1,MPI_DOUBLE,mpi_ntop,13,
										&u_nn[i][j][0],1,MPI_DOUBLE,mpi_nbottom,13,mpi_comm3d, &status);
				MPI_Sendrecv(&v_nn[i][j][1],1,MPI_DOUBLE,mpi_nbottom,14,
										&v_nn[i][j][mpi_KMAX-1],1,MPI_DOUBLE,mpi_ntop,14,mpi_comm3d, &status);
				MPI_Sendrecv(&v_nn[i][j][mpi_KMAX-2],1,MPI_DOUBLE,mpi_ntop,15,
										&v_nn[i][j][0],1,MPI_DOUBLE,mpi_nbottom,15,mpi_comm3d, &status);
				MPI_Sendrecv(&w_nn[i][j][1],1,MPI_DOUBLE,mpi_nbottom,16,
										&w_nn[i][j][mpi_KMAX-1],1,MPI_DOUBLE,mpi_ntop,16,mpi_comm3d, &status);
				MPI_Sendrecv(&w_nn[i][j][mpi_KMAX-2],1,MPI_DOUBLE,mpi_ntop,17,
										&w_nn[i][j][0],1,MPI_DOUBLE,mpi_nbottom,17,mpi_comm3d, &status);
			}
		}
	}
}
//--------------------------update boundary conditions of velocity(V and W)------------------//
void calculation_2D::UpdateBC_velocity(boundary_location BC)
{
	int i,j,k;
	
	switch(BC){
	case WestBC:
		{// left
			if (mpi_coords[1]==0){
				for (i=0;i<mpi_IMAX;i++){
					for(k=0;k<mpi_KMAX;k++){
						u_nn[i][0][k]=u_nn[i][1][k];
						v_nn[i][0][k]=0.0;
						w_nn[i][0][k]=w_nn[i][1][k];
					}
				}
			}
		}
		break;
	case EastBC:
		{// right
			if (mpi_coords[1]==mpi_yarea-1){
				for (i=0;i<mpi_IMAX;i++){
					for(k=0;k<mpi_KMAX;k++){
						u_nn[i][mpi_JMAX-1][k]=u_nn[i][mpi_JMAX-2][k];
						v_nn[i][mpi_JMAX-1][k]=0.0;
						w_nn[i][mpi_JMAX-1][k]=w_nn[i][mpi_JMAX-2][k];
					}
				}
			}
		}

		break;
	case SouthBC:
		{ // bottom
			if (mpi_coords[2]==0){
				for (i=0;i<mpi_IMAX;i++){
					for(j=0;j<mpi_JMAX;j++){
						u_nn[i][j][0]=u_nn[i][j][1];
						v_nn[i][j][0]=v_nn[i][j][1];
						w_nn[i][j][0]=0.0;
					}
				}
			}
		}
		break;
	case NorthBC:
		{// top
			if (mpi_coords[2]==mpi_zarea-1){
				for (i=0;i<mpi_IMAX;i++){
					for(j=0;j<mpi_JMAX;j++){
						u_nn[i][j][mpi_KMAX-1]=u_nn[i][j][mpi_KMAX-2];
						v_nn[i][j][mpi_KMAX-1]=v_nn[i][j][mpi_KMAX-2];
						w_nn[i][j][mpi_KMAX-1]=0.0;
					}
				}
			}
		}
		break;
	case FrontBC:
		{// Front
			if (mpi_coords[0]==0){
				for (j=0;j<mpi_JMAX;j++){
					for(k=0;k<mpi_KMAX;k++){
						u_nn[0][j][k]=0.0;
						v_nn[0][j][k]=v_nn[1][j][k];
						w_nn[0][j][k]=w_nn[1][j][k];
					}
				}
			}
		}
		break;
	case BackBC:
		{// Back
			if (mpi_coords[0]==mpi_xarea-1){
				for (j=0;j<mpi_JMAX;j++){
					for(k=0;k<mpi_KMAX;k++){
						u_nn[mpi_IMAX-1][j][k]=0.0;
						v_nn[mpi_IMAX-1][j][k]=v_nn[mpi_IMAX-2][j][k];
						w_nn[mpi_IMAX-1][j][k]=w_nn[mpi_IMAX-2][j][k];
					}
				}
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary, velcotiy BC = " << BC << endl;
	}
}
void calculation_2D::UpdateBC_velocity(boundary_location BC,calculation_2D& porous)
{
	//update absorbing boundary
	int i,j;
	double temp_u,temp_v,temp_w;
	int AirDx_CouplingDx,AirDy_CouplingDy,AirDz_CouplingDz;
	if (mpi_coords[2]!=0) return;
	AirDx_CouplingDx=(int)(diff_x1/porous.diff_x1);
	AirDy_CouplingDy=(int)(diff_y1/porous.diff_y1);
	AirDz_CouplingDz=(int)(diff_z1/porous.diff_z1);
	switch(BC){
	case SouthBC:
		if(AirDz_CouplingDz<5.0){
			for(i=0;i<mpi_IMAX;i++){
				for (j=0;j<mpi_JMAX;j++){
					u_nn[i][j][0]=porous.u_nn[i][j][porous.mpi_KMAX-2];
					porous.u_nn[i][j][porous.mpi_KMAX-1]=u_nn[i][j][1];
					v_nn[i][j][0]=porous.v_nn[i][j][porous.mpi_KMAX-2];
					porous.v_nn[i][j][porous.mpi_KMAX-1]=v_nn[i][j][1];
					w_nn[i][j][0]=porous.w_nn[i][j][porous.mpi_KMAX-2];
					porous.w_nn[i][j][porous.mpi_KMAX-1]=w_nn[i][j][1];
				}
			}
		}else{
		/*
			//transient interface
			for(i=0;i<(IMAX-1);i++){
				for (j=0;j<(JMAX-1);j++){
					for(int k=0;k<(AirDx_CouplingDx);k++){
						temp_u=(u_nn[i+1][j][1]-u_nn[i][j][1])/(AirDx_CouplingDx)*k+u_nn[i][j][1];
						porous.u_nn[i*AirDx_CouplingDx+k][j][porous.KMAX-1]=
							(temp_u-porous.u_nn[i*AirDx_CouplingDx+k][j][porous.KMAX-2])/
							(AirDz_CouplingDz+1)+porous.u_nn[i*AirDx_CouplingDx+k][j][porous.KMAX-2];

						temp_v=(v_nn[i+1][j][1]-v_nn[i][j][1])/(AirDx_CouplingDx)*k+v_nn[i][j][1];
						porous.v_nn[i*AirDx_CouplingDx+k][j][porous.KMAX-1]=
							(temp_v-porous.v_nn[i*AirDx_CouplingDx+k][j][porous.KMAX-2])/
							(AirDz_CouplingDz+1)+porous.v_nn[i*AirDx_CouplingDx+k][j][porous.KMAX-2];

						temp_w=(w_nn[i+1][j][1]-w_nn[i][j][1])/(AirDx_CouplingDx)*k+w_nn[i][j][1];
						porous.w_nn[i*AirDx_CouplingDx+k][j][porous.KMAX-1]=
							(temp_w-porous.w_nn[i*AirDx_CouplingDx+k][j][porous.KMAX-2])/
							(AirDz_CouplingDz+1)+porous.w_nn[i*AirDx_CouplingDx+k][j][porous.KMAX-2];
					}

					for(int k=0;k<(AirDy_CouplingDy);k++){
						temp_u=(u_nn[i][j+1][1]-u_nn[i][j][1])/(AirDy_CouplingDy)*k+u_nn[i][j][1];
						porous.u_nn[i][j*AirDy_CouplingDy+k][porous.KMAX-1]=
							(temp_u-porous.u_nn[i][j*AirDy_CouplingDy+k][porous.KMAX-2])/
							(AirDz_CouplingDz+1)+porous.u_nn[i][j*AirDy_CouplingDy+k][porous.KMAX-2];

						temp_v=(v_nn[i+1][j][1]-v_nn[i][j][1])/(AirDy_CouplingDy)*k+v_nn[i][j][1];
						porous.v_nn[i][j*AirDy_CouplingDy+k][porous.KMAX-1]=
							(temp_v-porous.v_nn[i][j*AirDy_CouplingDy+k][porous.KMAX-2])/
							(AirDz_CouplingDz+1)+porous.v_nn[i][j*AirDy_CouplingDy+k][porous.KMAX-2];

						temp_w=(w_nn[i+1][j][1]-w_nn[i][j][1])/(AirDy_CouplingDy)*k+w_nn[i][j][1];
						porous.w_nn[i][j*AirDy_CouplingDy+k][porous.KMAX-1]=
							(temp_w-porous.w_nn[i][j*AirDy_CouplingDy+k][porous.KMAX-2])/
							(AirDz_CouplingDz+1)+porous.w_nn[i][j*AirDy_CouplingDy+k][porous.KMAX-2];
					}
				}
			}
			for(i=1;i<(IMAX-1);i++){
				for (j=0;j<(JMAX-1);j++){
					u_nn[i][j][0]=porous.u_nn[i*AirDx_CouplingDx][j][porous.KMAX-1]
					v_nn[i][j][0]=porous.v_nn[i*AirDx_CouplingDx][j][porous.KMAX-1];
					w_nn[i][j][0]=porous.w_nn[i*AirDx_CouplingDx][j][porous.KMAX-1];
				}
			}
		*/
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary, when UpdateBC_velocity, BC = " << BC << endl;
	}
}

void calculation_2D::UpdateBC_pressure(boundary_location BC,calculation_2D& porous)
{
	//update absorbing boundary
	int i,j;
	double temp_p;
	int AirDx_CouplingDx,AirDy_CouplingDy,AirDz_CouplingDz;
	if (mpi_coords[2]!=0) return;
	AirDx_CouplingDx=(int)(diff_x1/porous.diff_x1);
	AirDy_CouplingDy=(int)(diff_y1/porous.diff_y1);
	AirDz_CouplingDz=(int)(diff_z1/porous.diff_z1);

	switch(BC){
	case SouthBC:
		if(AirDz_CouplingDz<5.0){
			for(i=0;i<mpi_IMAX;i++){
				for (j=0;j<mpi_JMAX;j++){
					p_nn[i][j][0]=porous.p_nn[i][j][porous.mpi_KMAX-2];
					porous.p_nn[i][j][porous.mpi_KMAX-1]=p_nn[i][j][1];
				}
			}
		}else{
			//transient interface
/*
			for(i=0;i<(IMAX-1);i++){
				for(int k=0;k<(AirDy_CouplingDy);k++){
					temp_p=(p_nn[i+1][1]-p_nn[i][1])/(AirDy_CouplingDy)*k+
						p_nn[i][1];
					porous.p_nn[i*AirDy_CouplingDy+k][porous.JMAX-1]=
						(temp_p-porous.p_nn[i*AirDy_CouplingDy+k][porous.JMAX-2])/
						(AirDz_CouplingDz+1)+porous.p_nn[i*AirDy_CouplingDy+k][porous.JMAX-2];
				}
			}
			for(i=1;i<(IMAX-1);i++){
				p_nn[i][0]=porous.p_nn[i*AirDy_CouplingDy][porous.JMAX-1];
			}
*/
		}

		break;
	default:
		cout<<"it is wrong to entrance boundary when UpdateBC_pressure, BC = " << BC <<endl;
	}
}
void calculation_2D::mpi_send_data(int N)
{
	int i,j,k,i1,i2,j1,j2,k1,k2;
	int in,jn,kn,ii,jj,kk,record_n;

	if (mpi_rank==0){
		record_n=0;
		for(i=0;i<IMAX;i++){
			if (i%out_n==0||i==IMAX-1){
				X_out[record_n]=axisbound_X[i];
				record_n=record_n+1;
			}
		}
		record_n=0;
		for(j=0;j<JMAX;j++){
			if (j%out_n==0||j==JMAX-1){
				Y_out[record_n]=axisbound_Y[j]+N*(move_frame.trail_DJ-1)*diff_y1;
				record_n=record_n+1;
			}
		}
		record_n=0;
		for(k=0;k<KMAX;k++){
			if (k%out_n==0||k==KMAX-1){
				Z_out[record_n]=axisbound_Z[k];
				record_n=record_n+1;
			}
		}
	}
	record_n=0;
	for(i=0;i<IMAX_out;i++){
		for(j=0;j<JMAX_out;j++){
			for(k=0;k<KMAX_out;k++){
				p_in[record_n]=0.0;
				record_n=record_n+1;
			}
		}
	}

	if (judge_porous1==0||(judge_porous1==1 && mpi_coords[2]==0)){
		if (mpi_coords[0]!=0)i1=1;
		else i1=0;
		if (mpi_coords[0]!=mpi_xarea-1)i2=1;
		else i2=0;
	
		if (mpi_coords[1]!=0)j1=1;
		else j1=0;
		if (mpi_coords[1]!=mpi_yarea-1)j2=1;
		else j2=0;

		if (mpi_coords[2]!=0)k1=1;
		else k1=0;
		if (mpi_coords[2]!=mpi_zarea-1)k2=1;
		else k2=0;

		record_n=0;
		for(k=0;k<KMAX;k++ ){
			if (k%out_n==0||k==KMAX-1){
				for(j=0;j<JMAX;j++ ){
					if (j%out_n==0||j==JMAX-1){
						for(i=0;i<IMAX;i++ ){
							if (i%out_n==0||i==IMAX-1){
								if (i>=mpi_i1+i1 && i<=mpi_i2-i2 &&
									j>=mpi_j1+j1 && j<=mpi_j2-j2 &&
									k>=mpi_k1+k1 && k<=mpi_k2-k2){
									ii=i-mpi_i1;jj=j-mpi_j1;
									kk=k-mpi_k1;
									p_in[record_n]=p_nn[ii][jj][kk];
								}else p_in[record_n]=0.0;
								record_n=record_n+1;
							}
						}
					}
				}
			}
		}
	}
	MPI_Reduce(p_in,p_out,IMAX_out*JMAX_out*KMAX_out,MPI_DOUBLE,MPI_SUM,0,mpi_comm3d);

	/*
	for(i=0;i<IMAX;i++ ){
		if (i%out_n==0||i==IMAX-1){
			if (i!=IMAX-1)in=i/out_n;
			else in=IMAX_out-1;
			for(j=0;j<JMAX;j++ ){
				if (j%out_n==0||j==JMAX-1){
					if (j!=JMAX-1)jn=j/out_n;
					else jn=JMAX_out-1;
					for(k=0;k<KMAX;k++ ){
						if (k%out_n==0||k==KMAX-1){
							if (k!=KMAX-1)kn=k/out_n;
							else kn=KMAX_out-1;
							if (i>=mpi_i1+i1 && i<=mpi_i2-i2 &&
								j>=mpi_j1+j1 && j<=mpi_j2-j2 &&
								k>=mpi_k1+k1 && k<=mpi_k2-k2){
								ii=i-mpi_i1;jj=j-mpi_j1;
								kk=k-mpi_k1;
								p_in[in][jn][kn]=p_nn[ii][jj][kk];
							}else p_in[in][jn][kn]=0.0;
						}
					}
				}
			}
		}
	}
	mpi_comm3d.Reduce(p_in,p_out,IMAX_out*JMAX_out*KMAX_out,MPI_DOUBLE,MPI::SUM,0);
	*/
	/*
	for (i=0;i<IMAX_out;i++){
		for (j=0;j<JMAX_out;j++){
			for (k=0;k<KMAX_out;k++){
				p_out[i][j][k]=0.0;
			}
		}
	}
	for(i=i1;i<mpi_IMAX-i2;i++){
		ii=mpi_i1+i;
		if (ii%out_n==0||ii==IMAX-1){
			if (ii!=IMAX-1)in=ii/out_n;
			else in=IMAX_out-1;
			for(j=j1;j<mpi_JMAX-j2;j++){
				jj=mpi_j1+j;
				if (jj%out_n==0||jj==JMAX-1){
					if (jj!=JMAX-1)jn=jj/out_n;
					else jn=JMAX_out-1;
					for(k=k1;k<mpi_KMAX-k2;k++){
						kk=mpi_k1+k;
						if (kk%out_n==0||kk==KMAX-1){
							if (kk!=KMAX-1)kn=kk/out_n;
							else kn=KMAX_out-1;
							mpi_comm3d.Reduce(&p_nn[i][j][k],&p_out[in][jn][kn],1,MPI_DOUBLE,MPI::SUM,0);
						}
					}
				}
			}
		}
	}
	*/
}
void getOstream(ostream& output_stream, calculation_2D& my_object)
{
	int IMAX,JMAX,KMAX;

	IMAX=my_object.IMAX_out;JMAX=my_object.JMAX_out;KMAX=my_object.KMAX_out;
	output_stream << "I=" << IMAX << ",J=" << JMAX << ",K=" << KMAX << ",F=POINT" << endl;
	int record_n=0;
	for(int k=0;k<KMAX;k++){
		for(int j=0;j<JMAX;j++){
			for(int i=0;i<IMAX;i++){
				output_stream << my_object.X_out[i] << " " << my_object.Y_out[j] << " "
					<< my_object.Z_out[k]<< " " << my_object.p_out[record_n] << " " << endl;
				record_n=record_n+1;			
			}
		}
	}
}

ostream& operator<<(ostream& output_stream, calculation_2D& my_object)
{
	getOstream(output_stream,my_object);
	return output_stream;
}
