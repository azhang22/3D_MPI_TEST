//++++++++++++++++++++++++++++filename: airpore.cpp ++++++++++++++++++++++++++//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "airpore.h"
#include "air.h"
#include "porous.h"
#include "time.h"
#include <direct.h> // for windows


//****************************airpore code*******************************//
airpore::airpore(char *input_file)
{
	ifstream infile(input_file);
	infile>>gauss_width;infile.ignore(100,'\n');
	infile>>move_frame.JMAX>>move_frame.lead_DJ>>move_frame.trail_DJ>>move_frame.Judge>>move_frame.ratio;
	infile.ignore(100,'\n');
	int scheme_type;
	infile>>scheme_type;infile.ignore(100,'\n');
	switch(scheme_type)
	{
	case 1:
		scheme=FB_v_p;break;
	case 2:
		scheme=FB_p_v;break;
	case 3:
		scheme=FB_vp;break;
	case 4:
		scheme=LeapTrap;break;
	case 5:
		scheme=RK4;break;
	case 6:
		scheme=RK2;break;
	case 7:
		scheme=ABM;break;
	default:
		scheme=FB_p_v;
	}
	int *boundary_choice,*num_PML;
	boundary_choice=new int [6];
	num_PML=new int [4];
	int i;
	for(i=0;i<6;i++) infile>>boundary_choice[i];
	infile.ignore(100,'\n');
	infile.ignore(100,'\n');
	for(i=0;i<3;i++)num_PML[i]=0;
	for(i=0;i<6;i++){
		switch(boundary_choice[i]){
		case 1://rigid
		{
			switch(i){
			case 0:
				air_boundary.air_north=rigid;break;
			case 1:
				air_boundary.air_south=rigid;break;
			case 2:
				air_boundary.air_east=rigid;break;
			case 3:
				air_boundary.air_west=rigid;break;
			case 4:
				air_boundary.air_front=rigid;break;
			case 5:
				air_boundary.air_back=rigid;break;
			default:
				cout<<"please entrance your choice again!";
			}
		}
		break;
		case 2://absorbing layer
		{
			switch(i){
			case 0:
				air_boundary.air_north=absorbing;
				num_PML[0]=num_PML[0]+1;
				break;
			case 1:
				air_boundary.air_south=absorbing;
				num_PML[0]=num_PML[0]+10;
				break;
			case 2:
				air_boundary.air_east=absorbing;
				num_PML[1]=num_PML[1]+1;
				break;
			case 3:
				air_boundary.air_west=absorbing;
				num_PML[1]=num_PML[1]+10;
				break;
			case 4:
				air_boundary.air_front=absorbing;
				num_PML[2]=num_PML[2]+1;
				break;
			case 5:
				air_boundary.air_back=absorbing;
				num_PML[2]=num_PML[2]+10;
				break;
			default:
				cout<<"please entrance your choice again!";
			}
		}
		break;
		case 3://porous media
		{
			switch(i){
			case 0:
				air_boundary.air_north=porous_media;break;
			case 1:
				air_boundary.air_south=porous_media;break;
			case 2:
				air_boundary.air_east=porous_media;break;
			case 3:
				air_boundary.air_west=porous_media;break;
			case 4:
				air_boundary.air_front=porous_media;break;
			case 5:
				air_boundary.air_back=porous_media;break;
			default:
				cout<<"please entrance your choice again!";
			}
		}
		break;
		case 4: //radiation
		{
			switch(i){
			case 0:
				air_boundary.air_north=radiation;break;
			case 1:
				air_boundary.air_south=radiation;break;
			case 2:
				air_boundary.air_west=radiation;break;
			case 3:
				air_boundary.air_east=radiation;break;
			case 4:
				air_boundary.air_back=radiation;break;
			case 5:
				air_boundary.air_front=radiation;break;
			default:
				cout<<"please entrance your choice again!";
			}
		}
		break;
		default:
			cout<<"please entrance your choice again!";
		}
	}
	switch(air_boundary.air_north){
	case absorbing:case rigid:
		CaseNo=1;break;
	case porous_media:
		CaseNo=2;break;
	case radiation:
		CaseNo=3;break;
	default:
		cout<<"please entrance your choice again!";
	}

	infile>>time_domain>>out_difft>>FFT_m>>frequency;infile.ignore(100,'\n');
	FFT_N=1;
	for (i=0;i<FFT_m;i++) FFT_N *= 2;
	//air properties
	infile>>velocity_coef>>velocity_method;infile.ignore(100,'\n');
	infile>>AirPara.adiabatic_coef>>AirPara.Pav>>AirPara.sound_speed;infile.ignore(100,'\n');
	AirPara.Aav=pow(AirPara.sound_speed,2)/AirPara.adiabatic_coef/AirPara.Pav;
	//absorbing layer structure
	infile>>PML_xbd.alpha>>PML_xbd.resistivity;infile.ignore(100,'\n');
	infile>>PML_ybd.alpha>>PML_ybd.resistivity>>PML_ybd.time_step;infile.ignore(100,'\n');
	infile>>PML_zbd.alpha>>PML_zbd.resistivity;infile.ignore(100,'\n');
	for (i=0;i<3;i++){
		switch(num_PML[i]){
		case 0:
			switch(i){
			case 0:
				PML_zbd.judge=0;
				break;
			case 1:
				PML_ybd.judge=0;
				break;
			case 2:
				PML_xbd.judge=0;
				break;
			default:
				cout<<"please entrance your choice again!";
			}
		break;
		case 1:
			switch(i){
			case 0:
				PML_zbd.judge=2;
				break;
			case 1:
				PML_ybd.judge=2;
				break;
			case 2:
				PML_xbd.judge=2;
				break;
			default:
				cout<<"please entrance your choice again!";
			}
		break;
		case 10:
			switch(i){
			case 0:
				PML_zbd.judge=1;
				break;
			case 1:
				PML_ybd.judge=1;
				break;
			case 2:
				PML_xbd.judge=1;
				break;
			default:
				cout<<"please entrance your choice again!";
			}
		break;
		case 11:
			switch(i){
			case 0:
				PML_zbd.judge=3;
				break;
			case 1:
				PML_ybd.judge=3;
				break;
			case 2:
				PML_xbd.judge=3;
				break;
			default:
				cout<<"please entrance your choice again!";
			}
		break;
		default:
				cout<<"please entrance your choice again!";
		}
	}
	//porous media structure
	infile>>PorePara.q>>PorePara.porosity>>PorePara.resistivity;infile.ignore(100,'\n');
	PorePara.eff_density=PorePara.q*PorePara.q/PorePara.porosity/AirPara.Aav;
	PorePara.Kp=1.0/AirPara.adiabatic_coef/AirPara.Pav;

	//source and receiver
	infile>>source.x>>source.y>>source.z>>receiver.x>>receiver.y>>receiver.z;infile.ignore(100,'\n');

	//air space dimension
	infile>>AirStep.diff_t;infile.ignore(100,'\n');
	infile>>AirStep.diff_x>>AirStep.diff_y>>AirStep.diff_z;infile.ignore(100,'\n');
	
	// input domain of size for decomposition
	infile>>mpi_var.mpi_xarea>>mpi_var.mpi_yarea>>mpi_var.mpi_zarea;
	infile.ignore(100,'\n');
	// output.data
	infile.ignore(100,'\n');infile.ignore(100,'\n');
	infile>>restart>>restart_out>>restart_infile;

	// 1 50 ct003

	// curve setting,adds the porous mediea hill,
	/* curve_judge: judge the curve or not; curve_coefa,curve_coefb, curve_coefc
	and curve_coefd: The coefficient of equations:y=a*x^3+b*x^2+c*x+d;
	*/

	infile.ignore(100,'\n');infile.ignore(100,'\n');
	infile>>hill.curve_judge>>hill.curve_coefa>>hill.curve_coefb>>
			  hill.curve_coefc>>hill.curve_coefd>>hill.curve_coefe>>
			  hill.curve_coeff>>hill.curve_coefg>>hill.curve_coefh>>
			  hill.curve_coefi>>hill.curve_coefj>>hill.curve_coefk>>
			  hill.curve_coefl>>hill.curve_points;
	infile.ignore(100,'\n');

	infile >> hill.q >> hill.porosity >> hill.resistivity; infile.ignore(100, '\n');
	hill.eff_density=hill.q * hill.q /hill.porosity/AirPara.Aav;
	hill.Kp=1.0/AirPara.adiabatic_coef/AirPara.Pav;

	infile >> canopy.q >> canopy.porosity >> canopy.resistivity; infile.ignore(100, '\n');
	canopy.eff_density = canopy.q *canopy.q / canopy.porosity / AirPara.Aav;
	canopy.Kp = 1.0 / AirPara.adiabatic_coef / AirPara.Pav;

	infile >> trunk.q >> trunk.porosity >> trunk.resistivity; infile.ignore(100, '\n');
	trunk.eff_density = trunk.q * trunk.q / trunk.porosity / AirPara.Aav;
	trunk.Kp = 1.0 / AirPara.adiabatic_coef / AirPara.Pav;

	infile>>hill.curve_xstar>>hill.curve_xstop>>hill.curve_ystar>>hill.curve_ystop
							>>hill.curve_zstar>>hill.curve_zstop;
	infile.ignore(100,'\n');
	//porous media setting
	infile.ignore(100,'\n');infile.ignore(100,'\n');
	if(air_boundary.air_south==porous_media)
	{
		infile>>SouthStep.diff_x>>SouthStep.diff_y>>SouthStep.diff_z;
		SouthStep.diff_t=AirStep.diff_t;
		infile.ignore(100,'\n');
	}
	else
	{
		infile.ignore(100,'\n');infile.ignore(100,'\n');
	}

	infile.close();
	delete[] boundary_choice;
	delete[] num_PML;
}


airpore::~airpore()
{
}


void airpore::Cal_velocity()
{
	AirMedia->Cal_velocity();
	AirMedia->MPIexchange_velocity();
	switch(CaseNo){
	case 1:
		{
			AirMedia->UpdateBC_velocity(NorthBC);
			AirMedia->UpdateBC_velocity(EastBC);
			AirMedia->UpdateBC_velocity(WestBC);
			AirMedia->UpdateBC_velocity(FrontBC);
			AirMedia->UpdateBC_velocity(BackBC);
			if(air_boundary.air_south==rigid || air_boundary.air_south==absorbing)
			{
				AirMedia->UpdateBC_velocity(SouthBC);
			}
			else
			{
				if (mpi_var.mpi_coords[2]==0){
					south_pore->Cal_velocity();
					south_pore->MPIexchange_velocity();
					south_pore->UpdateBC_velocity(SouthBC);
					south_pore->UpdateBC_velocity(EastBC);
					south_pore->UpdateBC_velocity(WestBC);
					south_pore->UpdateBC_velocity(FrontBC);
					south_pore->UpdateBC_velocity(BackBC);
					AirMedia->UpdateBC_velocity(SouthBC,*south_pore);
				}
			}
			AirMedia->Update_PML_Qvw();
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::Cal_pressure()
{
	AirMedia->Cal_pressure();
	AirMedia->MPIexchange_pressure();
	switch(CaseNo){
	case 1:
		{
			AirMedia->UpdateBC_pressure(NorthBC);
			AirMedia->UpdateBC_pressure(EastBC);AirMedia->UpdateBC_pressure(WestBC);
			AirMedia->UpdateBC_pressure(FrontBC);AirMedia->UpdateBC_pressure(BackBC);
			if(air_boundary.air_south==rigid || air_boundary.air_south==absorbing){
				AirMedia->UpdateBC_pressure(SouthBC);
			}else{
				if (mpi_var.mpi_coords[2]==0){
					south_pore->Cal_pressure();
					south_pore->MPIexchange_pressure();
					south_pore->UpdateBC_pressure(SouthBC);
					south_pore->UpdateBC_pressure(EastBC);
					south_pore->UpdateBC_pressure(WestBC);
					south_pore->UpdateBC_pressure(FrontBC);
					south_pore->UpdateBC_pressure(BackBC);
					AirMedia->UpdateBC_pressure(SouthBC,*south_pore);
				}
			}
			AirMedia -> Update_PML_Qp();
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::get_output(int mpi_rank1,int mpi_size1)
{
	mpi_rank=mpi_rank1;
	mpi_size=mpi_size1;
	//initialize all computational domains
	MPI_Initialize();
	SetInitialCond();
	//output coordinate (combine all domains into one coordinate system)
	//output_coordi();//for plot3d in Tecplot
	//main calculation
	double **pressure_temp;
	int *temp_i,*temp_j,*temp_k,IMAX,JMAX,KMAX,NO_Rec_y;
	int i,j,index_pressure;
	Rec_count=5;
	pressure_temp=new double *[Rec_count];
	for(i=0;i<Rec_count;i++) pressure_temp[i]=new double[time_domain];
	for (i=0;i<Rec_count;i++){
		for (j=0;j<time_domain;j++){
			pressure_temp[i][j]=0.0;
		}
	}
	IMAX = AirMedia->get_whole_IMAX();
	JMAX = AirMedia->get_whole_JMAX();
	KMAX = AirMedia->get_whole_KMAX();
	if(receiver.y==0){
		NO_Rec_y=1;
	}else{
		NO_Rec_y=Rec_count;
	}
	temp_i=new int [NO_Rec_y];temp_j=new int [NO_Rec_y];
	temp_k=new int [NO_Rec_y];
	for(i=0;i<NO_Rec_y;i++){
		int ii,jj,kk;
		struct Position rec;
		if (i == 0) { rec.x = receiver.x; rec.y = -0.23; rec.z = 1.75; }			//
		else if (i == 1){ rec.x = receiver.x; rec.y = -0.46; rec.z = 0.46; }
		else if (i == 2){ rec.x = receiver.x; rec.y = -0.1; rec.z = 0.1; }
		else if (i == 3){ rec.x = receiver.x; rec.y = 1.75; rec.z = 0.23; }
		else if (i == 4){ rec.x = receiver.x; rec.y = 1.75; rec.z = 0.46; }
		else if (i == 5){ rec.x = receiver.x; rec.y = 0.6; rec.z = 0.01; }
		else if (i == 6){ rec.x = receiver.x; rec.y = 0.7; rec.z = 0.01; }
		else if (i == 7){ rec.x = receiver.x; rec.y = 0.8; rec.z = 0.01; }
		else if (i == 8){ rec.x = receiver.x; rec.y = 0.9; rec.z = 0.01; }
		else if (i == 9){ rec.x = receiver.x; rec.y = 1.0; rec.z = 0.01; }
		else if (i == 10){ rec.x = receiver.x; rec.y = 1.2; rec.z = 0.01; }
		else if (i == 11){ rec.x = receiver.x; rec.y = 1.4; rec.z = 0.01; }
		else if (i == 12){ rec.x = receiver.x; rec.y = 1.6; rec.z = 0.01; }
		else if (i == 13){ rec.x = receiver.x; rec.y = 1.8; rec.z = 0.01; }
		else if (i == 14){ rec.x = receiver.x; rec.y = 2.0; rec.z = 0.01; }
		else{			
			rec.x = receiver.x;
			rec.y = receiver.y;
			rec.z = receiver.z;
		}
		AirMedia->get_position(ii,jj,kk,rec);
		if (ii == 0 && jj == 0 && kk == 0){ ii = 1; jj = 1; kk = 1; cout << "receiver # " << i << " out of domain!" << endl; }
		if (ii > IMAX || jj > JMAX || kk > KMAX){ ii = 1; jj = 1; kk = 1; cout << "receiver # " << i << " out of domain!" << endl; }
		temp_i[i]=ii;temp_j[i]=jj;temp_k[i]=kk;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0) cout << "Receiver position all set!" << endl;
	//mkdir("./solution", 0777);
	_mkdir("solution");// for windows
	//if it exists, remove it
	if (restart == 0){
		index_pressure = 0;
		if (mpi_rank == 0){
			for (i = 0; i < NO_Rec_y; i++){
				char tempnameid[10] = "";
				char ptyfilename[100] = ("./solution/pty");

				if (i < 10 && i >= 0) strcat_s(ptyfilename, "000");
				else if (i >= 10 && i < 100) strcat_s(ptyfilename, "00");
				else if (i >= 100 && i < 1000) strcat_s(ptyfilename, "0");
				sprintf_s(tempnameid, "%d", i);
				strcat_s(ptyfilename, tempnameid);// pty001
				strcat_s(ptyfilename, ".dat");	// pty001.dat

				ofstream ofout(ptyfilename, ios::out | ios::binary);
				if (ofout.is_open()) remove(ptyfilename);
				ofout.close();
			}
			_mkdir("myfile");
			//mkdir("/panfs/pfs.local/work/zczheng/j391z772/Results/3D/20180822_3D_ANSI_study", 0777);
			//mkdir("/panfs/pfs.local/work/zczheng/j391z772/Results/3D/20180822_3D_ANSI_study/Rigid_ground", 0777);
			//mkdir("/panfs/pfs.local/scratch/zczheng/j391z772/restart_backup/20180822_3D_ANSI_study", 0777);
			//mkdir("/panfs/pfs.local/scratch/zczheng/j391z772/restart_backup/20180822_3D_ANSI_study/Rigid_ground", 0777);
			//mkdir("/panfs/pfs.local/scratch/zczheng/j391z772/restart_backup/20180822_3D_ANSI_study/Rigid_ground/myfile", 0777);
		}
	}


	int MF_count=0,temp_time=0;//MF_count to count how many frames when moving
	int time_series1,time_series2,time_series3;
	int time_series4,time_series5,time_series6;
	time_series1=473600;time_series2=120000;time_series3=149600;
	time_series4=238000;time_series5=296800;time_series6=473600;

	if (restart == 1) {
		index_pressure = Read_restart();
		MPI_Barrier(MPI_COMM_WORLD);
		restart_timer = index_pressure;
		AirMedia->Set_Timer(restart_timer);
		//restart = 0;
	}
	
	ofstream outfile[5];

	if (mpi_rank == 0){
		for (i = 0; i < NO_Rec_y; i++){
			char tempnameid1[10] = "";
			char pty_infile[100] = ("./solution/pty");

			if (i < 10 && i >= 0) strcat_s(pty_infile, "000");
			else if (i >= 10 && i < 100) strcat_s(pty_infile, "00");
			else if (i >= 100 && i < 1000) strcat_s(pty_infile, "0");
			sprintf_s(tempnameid1, "%d", i);
			strcat_s(pty_infile, tempnameid1);	// pty001
			strcat_s(pty_infile, ".dat");	// pty001.dat
			if (mpi_rank != 0)sprintf_s(pty_infile, "%d", mpi_rank);

			outfile[i].open(pty_infile, ios::app);
			outfile[i].setf(ios::scientific, ios::floatfield);
			outfile[i].precision(6);
			if (mpi_rank != 0){
				outfile[i].close();
				remove(pty_infile);
			}
		}
		cout << "Receiver files all set!" << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	int n;
	int MF_limit;
	double MF_lefttime1;
	double MF_lefttime2;
	double temp;
	do {
		//when count=0, computation will cover whole frame; next, computation will cover half of frame
		//	double MF_timeinterval;
		if(MF_count==(int)((JMAX-move_frame.JMAX)/(move_frame.lead_DJ-1)+1.e-6)||(JMAX==move_frame.JMAX)){
			MF_limit=time_domain-temp_time;
		}else if(MF_count==0){
			MF_limit=(int)((move_frame.JMAX-1)*AirStep.diff_y/AirPara.sound_speed/AirStep.diff_t*move_frame.ratio);
		}else{
			MF_limit=(int)((move_frame.trail_DJ-1)*AirStep.diff_y/AirPara.sound_speed/AirStep.diff_t);
			MF_lefttime1=0.0;
			MF_lefttime2=0.0;
			if (MF_count%10==0){
				MF_lefttime1=(move_frame.trail_DJ-1)*AirStep.diff_y/AirPara.sound_speed/AirStep.diff_t-MF_limit;
				MF_lefttime1=MF_lefttime1*10.0;
			}
			if(MF_count%100==0){
				MF_lefttime2=(move_frame.trail_DJ-1)*AirStep.diff_y/AirPara.sound_speed/AirStep.diff_t-MF_limit;
				MF_lefttime2=MF_lefttime2*100-int(MF_lefttime2*10)*10;
			}
			MF_limit=MF_limit+int(MF_lefttime1)+int(MF_lefttime2);
			if ((temp_time+MF_limit)>time_domain)MF_limit=time_domain-temp_time;
		}

		SetMovingFrame(MF_count);

		MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == 0) cout << "Moving frame set!" << ", MF_limit = " << MF_limit << endl;
		for(n=0;n<MF_limit;n++){
			//print out data
			if (restart == 1){
				n = index_pressure;
				restart = 0;
			}
			index_pressure=index_pressure+1;
			for (i = 0; i < NO_Rec_y; i++){
				temp = AirMedia->get_pressure(temp_i[i], temp_j[i], temp_k[i], MF_count);
				if (mpi_rank == 0){
					if (fabs(temp)<1.0e-100) temp = 0.0;
					outfile[i] << temp << endl;
					pressure_temp[i][index_pressure - 1] = temp;
				}
			}

			if (mpi_rank == 0) cout << "Time step = " << index_pressure << ", #3 receiver reading is " << pressure_temp[2][index_pressure - 1];

			MPI_Barrier(MPI_COMM_WORLD);
			
			//calculate pressure and velocity to advance
			clock_t begin_calp, begin_calv, end_calv, end_temp;
			begin_calv = clock();
			if ((scheme == FB_p_v) || (scheme == LeapTrap)){
				//if (mpi_rank == 0) cout << "starting cal_pressure " << endl;
				Cal_pressure();
				//if (mpi_rank == 0) cout << "starting cal_velocity " << endl;
				Cal_velocity();
				end_calv = clock();
				if (mpi_rank == 0) { cout << ", one time step using time: " << (double(end_calv - begin_calv)) / CLOCKS_PER_SEC << "s." << endl;}
			}
			else if(scheme==FB_v_p){
				Cal_velocity();
				Cal_pressure();
			}else{
				Cal_pressure();
				Cal_velocity();
				UpdateInitialCond(MF_count);
				Cal_pressure();
				Cal_velocity();
			}

			//update sourece
			//double source_frequency;
			//source_frequency = 50;
			//UpdateSource(source_frequency, index_pressure, source);

			//update data for the next time
			UpdateInitialCond(MF_count);
			if(((n+temp_time+1)%out_difft)==0){//output of contour
				mpisend_data_contour(MF_count);
				if (mpi_rank==0) get_data_contour(n+temp_time+1);
			}
			// Save restart files
			if (index_pressure%restart_out == 0) Save_restart(index_pressure);
		}
		MF_count = MF_count + 1;
		temp_time = MF_limit + temp_time;
	}while((MF_count<=(int)((JMAX-move_frame.JMAX)/(move_frame.lead_DJ-1)+1.e-6))&&(temp_time<time_domain) && move_frame.Judge==1);
	for (i = 0; i < NO_Rec_y; i++) outfile[i].close();
	//delete new variabe to vacate memory
	delete AirMedia;
	switch(CaseNo){
	case 1:
		if(air_boundary.air_south==porous_media){
			delete south_pore;
		}
		break;
	default:
		cout<<"it is wrong to delete object";
	}
	delete[] temp_i;delete[] temp_j;delete[] temp_k;
	for(i=0;i<Rec_count;i++)delete[] pressure_temp[i];
	delete[] pressure_temp;
	delete[] mpi_var.mpi_coords;
}

void airpore::UpdateSource(double f0, int N, Position)
{
	AirMedia->UpdateSource(f0, N, source);
}

void airpore::UpdateInitialCond(int MF_count)
{
	AirMedia->UpdateInitialCond(MF_count);
	switch(CaseNo){
			case 1:
				if(air_boundary.air_south==porous_media){
					if (mpi_var.mpi_coords[2]==0){
						south_pore->UpdateInitialCond(MF_count);
					}
				}
				break;
			default:
				cout<<"it is wrong to entrance boundary";
	}
}

void airpore::Save_restart(int index)
{
	char restartfile[200] = "myfile/ct", restart_temp[20];
	char restartfile0[200] = "myfile/ct", restart_temp0[20];
	char restart_pt0[20], restart_pt1[200];
	int restart_nn,i,j;
	restart_nn = index / restart_out;
	if (restart_nn < 10) {
		strcat_s(restartfile, "00"); 
		strcat_s(restartfile0, "00");
	}
	if ((restart_nn >= 10) && (restart_nn < 100)) { 
		strcat_s(restartfile, "0");
		strcat_s(restartfile0, "0"); 
	}
	
	int restart_ii, restart_id,restart_ii0;
	restart_ii = sprintf_s(restart_temp, "%d", restart_nn); 
	strcat_s(restartfile, restart_temp);	// ct003
	restart_id = sprintf_s(restart_temp, "_%d", mpi_rank);
	strcat_s(restartfile, restart_temp);	// ct003_2
	strcat_s(restartfile, ".dat");		// ct003_2.dat

	ofstream re_outfile1(restartfile, ios::out | ios::binary);	
	re_outfile1 << mpi_size << " " << mpi_rank << " " << index << endl;
	re_outfile1.close();
	save_restartfile(restartfile);
	re_outfile1.close();

	if (restart_nn > 1) {
		restart_ii0 = sprintf_s(restart_temp0, "%d", restart_nn - 1);
		strcat_s(restartfile0, restart_temp0);	// ct001
		strcat_s(restartfile0, restart_temp);		// ct001_2
		strcat_s(restartfile0, ".dat");			// ct001_2.dat

		ifstream re_inf(restartfile0, ios::in | ios::binary);
		if (!re_inf && mpi_rank == 0){ cout << "Cannot open restart output file index - 1 !!!!! :(" << endl; }

		int ret;
		ret = remove(restartfile0);
		if (mpi_rank == 0) {
			if (ret == 0) cout << "Old restart files deleted ! :)" << endl;
			else cout << "Unable to delete old files!!!!!!! :(" << endl;
		}
	}
	// save receivers data
	if (mpi_rank == 0){
		for (i = 1; i <= Rec_count; i++)
		{
			char pt_0[20] = "pty";
			char pt_1[200] = "myfile/pty";
			char pt_2[200] = "myfile/pty";

			sprintf_s(restart_temp, "%d", i);	//
			strcat_s(pt_0, restart_temp);		//pty1
			strcat_s(pt_0, ".dat");			//pty1.dat


			strcat_s(pt_1, restart_temp);
			strcat_s(pt_2, restart_temp);
			strcat_s(pt_1, "_ct");	//	.../pty1_ct
			strcat_s(pt_2, "_ct");	//	.../pty1_ct

			if (restart_nn < 10) {
				strcat_s(pt_1, "00");
			}
			if ((restart_nn >= 10) && (restart_nn < 100)) {
				strcat_s(pt_1, "0");
			}

			sprintf_s(restart_temp, "%d", restart_nn);
			strcat_s(pt_1, restart_temp);	//.../pty1_ct003
			strcat_s(pt_1, ".dat");		//.../pty1_ct003.dat

			copyFile(pt_0, pt_1);

			if (restart_nn > 2)
			{
				if (restart_nn - 2 < 10) {
					strcat_s(pt_2, "00");
				}
				if ((restart_nn - 2 >= 10) && (restart_nn < 100)) {
					strcat_s(pt_2, "0");
				}
				sprintf_s(restart_temp, "%d", restart_nn - 2);
				strcat_s(pt_2, restart_temp);	//.../pty1_ct001
				strcat_s(pt_2, ".dat");		//.../pty1_ct001.dat

				int ret;
				ret = remove(pt_2);
				if (ret == 0) cout << pt_2 << " deleted ! :)" << endl;
				else cout << "Unable to delete old files!!!!!!! :(" << endl;
			}
			cout << pt_1 << " archieved!" << endl;
		}
	}
}

void airpore::copyFile(const char* src_loc, const char* dest_loc)
{
	fstream in, out;
	in.open(src_loc, fstream::in | fstream::binary);

	if (in.is_open())
	{
		out.open(dest_loc, fstream::out);

		char tmp;
		while (in.read(&tmp, 1))
		{
			out.write(&tmp, 1);
		}
		out.close();
	}
	else  cout << "can not open file " << src_loc << endl; 
	in.close();
}

void airpore::save_restartfile(char *restartfile)
{
	AirMedia->save_restart_cal(restartfile);
	AirMedia->save_restart_air(restartfile);
	switch(CaseNo){
	case 1:
		if (air_boundary.air_south==porous_media){
			if (mpi_var.mpi_coords[2]==0){
				south_pore->save_restart_cal(restartfile);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

int airpore::Read_restart()
{
	int i, j, restart_nn;
	int Re_size, Re_rank, index;
	int int_pois, *point_pois;
	char temp_string[20];
	
	*restart_infile1 = '\0';
	strcat_s(restart_infile1,"myfile/");
	strcat_s(restart_infile1, restart_infile);
	sprintf_s(temp_string, "_%d", mpi_rank);
	strcat_s(restart_infile1, temp_string);
	strcat_s(restart_infile1, ".dat");

	cout << "input restart file name is " << restart_infile << ", restart_infile full location is " << restart_infile1 << endl;

	ifstream infile;
	infile.open(restart_infile1, ios::in | ios::binary);
	if (!infile) cout << "Cannot open restart input files in mpi # "<< mpi_rank <<"!!!!! :(" << endl; 
	else cout << "file " << restart_infile1 << " loaded ! :) " << endl;

	infile >> Re_size >> Re_rank >> index;
	infile.ignore(100, '\n');

	restart_nn = index / restart_out;

	if (Re_size != mpi_size) cout << "Wrong MPI nodes size!!!!!!!!!!! :(" << endl;

	int_pois = infile.tellg();
	point_pois = &int_pois;
	infile.close();

	if (mpi_rank == 0) cout << "Loading input files ..." << endl;

	input_restartfile(point_pois);

	if (mpi_rank == 0) cout << "All input cal data loaded !!! :)" << endl;
	
	// Load receiver readings
	if (mpi_rank == 0){
		for (i = 1; i <= Rec_count; i++)
		{
			char pt_0[20] = "pty";
			char pt_1[200] = "myfile/pty";

			sprintf_s(temp_string, "%d", i);
			strcat_s(pt_0, temp_string);
			strcat_s(pt_0, ".dat");	// pt001.dat

			strcat_s(pt_1, temp_string);
			strcat_s(pt_1, "_");
			strcat_s(pt_1, restart_infile);	//	.../pty1_ct003
			strcat_s(pt_1, ".dat");		//.../pty1_ct003.dat

			copyFile(pt_1, pt_0);
			cout << pt_0 << " loaded!" << endl;
		}
	}

	if (mpi_rank == 0) cout << "All input receivers data loaded !!! :)" << endl;

	return index;
}

void airpore::input_restartfile(int *point_pois)
{
	AirMedia->input_restart_cal(restart_infile1,point_pois);
	AirMedia->input_restart_air(restart_infile1,point_pois);
	switch(CaseNo){
	case 1:
		if (air_boundary.air_south==porous_media){
			if (mpi_var.mpi_coords[2]==0){
				south_pore->input_restart_cal(restart_infile,point_pois);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}

}

void airpore::SetMovingFrame(int MF_count)
{
	AirMedia->SetMovingFrame(MF_count);
	AirMedia->SetWindProfile(MF_count);
	AirMedia->SetQMove(MF_count);
	switch(CaseNo){
	case 1:
		if(air_boundary.air_south==porous_media){
			if (mpi_var.mpi_coords[2]==0){
				south_pore->SetMovingFrame(MF_count);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}
void airpore::MPI_Initialize()
{
	int *mpi_dim;
	//bool *mpi_periodic,mpi_reorder;
	mpi_dim=new int[3];
	//mpi_periodic=new bool[3];
	mpi_var.mpi_coords=new int[3];
	//mpi_reorder=true;
	int mpi_periodic[3];
	mpi_periodic[0] = 0; mpi_periodic[1] = 0; mpi_periodic[2] = 0;
	int mpi_reorder;
	mpi_reorder = 1;
	
	mpi_dim[0]=mpi_var.mpi_xarea;
	mpi_dim[1]=mpi_var.mpi_yarea;
	mpi_dim[2]=mpi_var.mpi_zarea;
	//mpi_periodic[0]=false;
	//mpi_periodic[1]=false;
	//mpi_periodic[2]=false;
	//mpi_var.mpi_comm3d=MPI::COMM_WORLD.Create_cart(3,mpi_dim,mpi_periodic,mpi_reorder);

	MPI_Cart_create(MPI_COMM_WORLD, 3, mpi_dim, mpi_periodic, mpi_reorder, &mpi_var.mpi_comm3d);

	mpi_var.mpi_rank = mpi_rank;
	//= mpi_var.mpi_comm3d.Get_rank();
	mpi_var.mpi_size = mpi_size;
	//mpi_var.mpi_comm3d.Shift(0,1,mpi_var.mpi_nback,mpi_var.mpi_nfront);
	//mpi_var.mpi_comm3d.Shift(1,1,mpi_var.mpi_nleft,mpi_var.mpi_nright);
	//mpi_var.mpi_comm3d.Shift(2,1,mpi_var.mpi_nbottom,mpi_var.mpi_ntop);

	MPI_Cart_shift(mpi_var.mpi_comm3d, 0, 1, &mpi_var.mpi_nback, &mpi_var.mpi_nfront);
	MPI_Cart_shift(mpi_var.mpi_comm3d, 1, 1, &mpi_var.mpi_nleft, &mpi_var.mpi_nright);
	MPI_Cart_shift(mpi_var.mpi_comm3d, 2, 1, &mpi_var.mpi_nbottom, &mpi_var.mpi_ntop);
	//mpi_var.mpi_comm3d.Get_coords(mpi_var.mpi_rank,3,mpi_var.mpi_coords);
	MPI_Cart_coords(mpi_var.mpi_comm3d, mpi_var.mpi_rank, 3, mpi_var.mpi_coords);
	delete [] mpi_dim; //delete[] mpi_periodic;
}

void airpore::SetInitialCond()
{
	int judge_porous,judge_source;
	char filename[100];
	strcpy_s(filename,"coordi_air");
	judge_porous=0;
	AirMedia= new air(scheme,AirStep,filename,move_frame,gauss_width,
		velocity_coef,velocity_method,AirPara,PML_xbd,
						PML_ybd,PML_zbd,hill,canopy,trunk,mpi_var,judge_porous,restart_timer);
	if (restart == 1) judge_source = 0;
	else judge_source = 1;
	AirMedia->Set_InitialCond(source, judge_source);

	switch(CaseNo){
	case 1:
		if(air_boundary.air_south==porous_media){
			strcpy_s(filename,"coordi_pore");
			judge_porous=1;
			south_pore=new porous(scheme,SouthStep,filename,move_frame,
								gauss_width,PorePara,mpi_var,judge_porous);
			if (mpi_var.mpi_coords[2]==0){
				south_pore->Set_InitialCond(source, judge_source);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::get_FFT_y1(int out_type)
{
	int i, j, num;
	double dt;
	num = Rec_count;
	num = 5;
	dt = AirStep.diff_t;

	for (j = 0; j < num; j++){
		char pty_file[100] = "./solution/pty";
		char p_f_file[100] = "./solution/p_f2_";
		char temp[100] = "";
		if (j <10)	{sprintf_s(temp, "000%d", j);}
		else if (j>=10 && j < 100) {sprintf_s(temp,"00%d",j);}
		else if (j>=100 && j < 1000) {sprintf_s(temp,"0%d",j);}
		else if (j > 1000) {sprintf_s(temp,"%d",j);}
		
		strcat_s(temp, ".dat");
		strcat_s(pty_file, temp);
		cout << "start to read " << pty_file << endl;
		ifstream infile(pty_file, ios::in | ios::binary);
		double *ppr, *pr, *pi;
		ppr = new double[time_domain];
		pi = new double[FFT_N];
		pr = new double[FFT_N];
		for (i = 0; i < time_domain; i++)ppr[i] = 0.0;
		for (i = 0; i < FFT_N; i++){
			pr[i] = 0.0;
			pi[i] = 0.0;
		}
		for (i = 0; i < time_domain; i++){
			infile >> ppr[i];
			infile.ignore(100, '\n');
		}
		infile.close();
		strcat_s(p_f_file, temp);
		cout << "start to write " << p_f_file << endl;
		if (receiver.y != 0){
			for (i = 0; i < time_domain; i++) { pr[i] = ppr[i]; pi[i] = 0; }
			for (i = time_domain; i < FFT_N; i++) { pr[i] = 0; pi[i] = 0; }
			FFT(1, FFT_m, pr, pi);
			FFT_output(FFT_m, pr, pi, p_f_file, dt);
		}
		delete[] pr; delete[] pi; delete[] ppr;
	}
}

//the following is the output format for the command "load data files" in Tecplot
void airpore::mpisend_data_contour(int MF_count){		
	AirMedia->mpi_send_data(MF_count);
	if(air_boundary.air_south==porous_media) south_pore->mpi_send_data(MF_count);
}
void airpore::get_data_contour(int n)
{
	char p_contour[200]="/panfs/pfs.local/work/zczheng/j391z772/Results/3D/20180822_3D_ANSI_study/Rigid_ground/pt",temp[20];
	int nn;
	nn=(n+1)/out_difft;
	if(nn<10) strcat_s(p_contour,"00");
	if((nn>=10)&&(nn<100)) strcat_s(p_contour,"0");
	int ii;
	ii=sprintf_s(temp,"%d",nn);
	strcat_s(p_contour,temp);strcat_s(p_contour,".dat");

	//output
	ofstream fp_p_t(p_contour,ios::out | ios::trunc | ios::binary);
	fp_p_t.setf(ios::scientific,ios::floatfield);
	fp_p_t.precision(6);

	switch(CaseNo){
	case 1:
		{
			if (air_boundary.air_south == rigid || air_boundary.air_south == absorbing){
				fp_p_t << "VARIABLES = \"X\", \"Y\", \"Z\", \"P\"" << endl;
				fp_p_t<<"ZONE T = \"Air Zone\",";
				fp_p_t<<*AirMedia;
			}else{

				fp_p_t<<"ZONE T = \"Air Zone\",";
				fp_p_t<<*AirMedia;
				fp_p_t<<"ZONE T = \"south pore Zone\",";
				fp_p_t<<*south_pore;
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
	fp_p_t.close();
}
