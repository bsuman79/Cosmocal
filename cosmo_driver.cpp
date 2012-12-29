#include <stdio.h>
#include <stdlib.h>

#include "cosmo.h"
#include "driver.h"
#include "growth.h"

using namespace std;
int main(){

   
	printf("\n\ncalculate cosmic time, ang. diameter distance, growth factor and its log derivative and few other cosmological quantities.\n input- cosmological parameters and redshift\n");
	cout << "*****************************************************"<<endl;
	
	/*****************************/
	
	cout<< "read in parameters"<<endl;
	cin >> H0;
	cin >> Omega_M;
	cin >> wt;
	cin >> z;
	
	cosmo cosm_model(H0, Omega_M, Omega_k, wt);
	growth growth_fact(Omega_M, Omega_k, wt, cosm_model); 
	
	cout << "scale factor= "<<cosm_model.scale_fact(z)<<endl;
	cout<<"cosmic time= "<<cosm_model.cosmic_time(z)<<" Gyr"<<endl;
	cout << "ang. diameter distance= "<<  cosm_model.ang_diam(z)<<" Mpc"<<endl;
	cout << "hubble parameter at z= "<<z<<" is "<<cosm_model.hubblez(z)<<" km/s/Mpc"<<endl;
	cout << "critical density= "<<cosm_model.calc_rho_crit(z)<<" Msun/Mpc^3 "<<endl;
	cout << "Omega matter at z= "<<z<<" is "<<cosm_model.Omega_Mz(z)<<endl;
	cout << "Delta_vir at z= "<<z<<" is "<<cosm_model.Delta_vir(z)<<endl;
	growth_fact.growthfactor(z, &Da, &Dadz);
	cout << "growth factor at z= "<<z<<" is "<<Da<<" and dlog D/d log a= "<<Dadz<<endl;
	return 0;
}

 
 
