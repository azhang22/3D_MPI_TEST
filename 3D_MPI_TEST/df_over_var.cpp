#include <math.h>
#include <stdio.h>
#include "mygeneral.h"

double Geo_over_var(double geo_var1,double geo_var2,double geo_var3,
					double geo_var4,double geo_jaco)
{  	
	return (geo_var1*geo_var2-geo_var3*geo_var4)*geo_jaco ;
}

/*
double Velo_over_var(double geo_var1,double geo_var2,double geo_var3,
			double Uav_up, double Uav_down,double Vav_up,double Vav_down,
			double Wav_up, double Wav_down,double geo_x,geo_y,geo_z)
{
	return geo_var1*(Uav_up-Uav_down)/geo_x/2.0+
			geo_var2*(Vav_up-Vav_down)/geo_y/2.0+
			geo_var3*(Wav_up-Wav_down)/geo_z/2.0;
}
double Coef_velo(double geo_var1,double geo_var2,double geo_var3,
						double Uav_vel,double Vav_vel,double Wav_vel,
						double geo_xyzt)
{ return (Uav_vel*geo_var1+Vav_vel*geo_var2+Wav_vel*geo_var3)/geo_xyzt;
}
*/