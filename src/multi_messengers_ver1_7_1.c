#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "../hdr/nr.h"
#include "../hdr/nrutil.c"
#include "../hdr/astro_const.h"
#include "../hdr/multi_messengers.h"


/* mode select */
const int opacitymode=3; //0=iron, 1=Y_e~0.3-0.5, 2=Y_e~0.1-0.2, 3=CO
const int diskwindmode=0; //0=no disk, 1=super-Edd disk, 2=disk wind
const int timebinmode=5; //0=MKKB15, 1=M+19, 2=OKM18, 3=OKM18-coarse, 4=M+20 (ALMA), 5=Galactic magnetar
const double tmax=51*86400;
const double tobs=(51+1.5*365.2422)*86400;


/* implicit parameters */
double gamma_b = 1.0e5; /* break Lorentz factor of electron injection spectrum */
double eps_mag,eps_e; /* magnetic field and electron acceleration efficiency */
double k_b_pwn; /* coefficient related to pdV work of PWN */
double Mdot_wind = 1.0e-6*M_SUN/(365.25*DAY_TO_SEC); /* mass loss rate of the progenitor star [g/s] */
double v_wind = 1.0e8; /* wind velocity of the progenito star [cm/s] */
double n_env = 100.; /* CSM gas number density [1/cc] */
double cs_csm = 1.0e7; /* CSM gas sound velocity [cm/s] */


int main ()
{
  struct _input in = load_inputfile();
  lc_sn(in.m_ns,in.r_ns,in.b_t,in.b_p,in.ome_0,in.ene_sn,in.m_ej,in.m_Ni,in.m_rpe,in.r_0,in.delta,in.d_l,in.albd_fac,in.kappa);
  
  return 0;
}


struct _input load_inputfile()
{
  double dmy[20];
  FILE *ip;
  ip = fopen("../in_para","r");
  int i=0;
  while (fscanf(ip,"%*s %lf %*[\n]",&dmy[i])!=EOF){
    i++;
  }
  fclose(ip);

  struct _input in;
  in.m_ns = dmy[0]*M_SUN;
  in.r_ns = dmy[1];
  in.b_t = dmy[2];
  in.b_p = dmy[3];
  in.ome_0 = 2.0*M_PI/(dmy[4]*1.0e-3);
  in.ene_sn = dmy[5];
  in.m_ej = dmy[6]*M_SUN;
  in.m_Ni = dmy[7]*M_SUN;
  in.m_rpe = dmy[8]*M_SUN;
  in.r_0 = dmy[9];
  in.delta = dmy[10];
  in.albd_fac = dmy[11];
  in.kappa = dmy[12];
  in.d_l = dmy[13]*PC_TO_CM;

  return in;
}


double i_ns(double m_ns, double r_ns)
{
    return 0.35*m_ns*pow(r_ns,2.0);
}


double ene_rot(double m_ns, double r_ns, double ome)
{
    return 0.50*i_ns(m_ns,r_ns)*pow(ome,2.0);
}


double ene_b(double m_ns, double r_ns, double b_t)
{
    double vol_NS = 4.0*M_PI/3.0*pow(r_ns,3.0);
    
    return pow(b_t,2.0)/8.0/M_PI*vol_NS;
}


double eps_b(double m_ns, double r_ns, double b_t)
{
    double ene_G = 3.0/4.0*G*pow(m_ns,2.0)/r_ns;
        
    return 15.0/4.0*ene_b(m_ns,r_ns,b_t)/ene_G;
}


double l_gw(double m_ns, double r_ns, double eps, double ome)
{
    /* Eq. (2) of Dall'Osso et al. (2009) */
    double sinalpha2 = 2./3.;//2./3.; for magnetars 0.; for mergers /* mean value of the inclination factor sin^2(alpha) */

    double i_ns_tmp = i_ns(m_ns,r_ns);
    
    return 2.0/5.0*G*pow(i_ns_tmp*eps,2.0)*pow(C,-5.0)*pow(ome,6.0)*sinalpha2*(1.0+15.0*sinalpha2);
}


double l_d(double r_ns, double b_p, double ome, double time)
{
    if(diskwindmode==0){
    //pulsar
    double sinalpha2 = 2./3.;//2./3.; for magnetars 0.; for mergers /* mean value of the inclination factor sin^2(alpha) */
    double mu_tmp = 0.5*b_p*pow(r_ns,3.0);
    
    double BHformation=1e20;//1e4;
    if(time<BHformation){
    return pow(mu_tmp,2.0)*pow(ome,4.0)*pow(C,-3.0)*(1.0+sinalpha2);
    }else{
    return 0;
    }
        
    }else{
    return 0.;
    }
    
}


double vol_ej(double r_ej)
{
    /* the volume of the ejecta */
    
    return 4.0*M_PI*pow(r_ej,3.0)/3.0;
}


double tau_ej(double m_ej, double delta, double r_ej, double kappa)
{
    /* the electron scattering optical depth in the ejecta */
    double tau_ej_tmp = (3.0-delta)/4.0/M_PI*kappa*m_ej/r_ej/r_ej/1.3;
    /* the factor 1.3 is put for calibrating the Arnett model */
    
    //if (tau_ej_tmp > 1.0)
        return tau_ej_tmp;
    //else
    //    return 1.0;
}


double rho_env(double r_ej)
{
  /* ambient matter density [g/cc] */
  double r_wind = sqrt(5.0*Mdot_wind*v_wind/24./M_PI/M_PRO/n_env/cs_csm/cs_csm);

  if (r_ej < r_wind) 
    return Mdot_wind/4.0/M_PI/r_ej/r_ej/v_wind;
  else 
    return 10.*M_PRO;
}


void lc_sn(double m_ns, double r_ns, double b_t, double b_p, double ome_0, double ene_sn, double m_ej, double m_Ni, double m_rpe, double r_0, double delta, double d_l, double albd_fac, double kappa)
{
    /* the condition for magnetically deformed rotation can occur (Eq. 19 of Dall'Osso+2009) */
    double ene_b_tmp = ene_b(m_ns,r_ns,b_t);
    double eps = 0.0;
    
    if (ene_b_tmp/1.0e50 < 2.1*(m_ns/1.4/M_SUN)*pow(2.0*M_PI/ome_0/1.0e-3,-2.0)*log(320.0*(m_ns/1.4/M_SUN)*pow(r_ns/1.2e6,-4)*pow(2.0*M_PI/ome_0/1.0e-3,2.0)*pow(b_p/1.0e14,-2.0)+1.0))
        eps = eps_b(m_ns,r_ns,b_t);
    else
        eps = 0.0;
    
    /* The spin down of newborn magnetar with the GW and dipole emissions */
    double k_gw = l_gw(m_ns,r_ns,eps,ome_0)/i_ns(m_ns,r_ns)/pow(ome_0,6.0);
    double k_d = l_d(r_ns,b_p,ome_0,0)/i_ns(m_ns,r_ns)/pow(ome_0,4.0);
    /* where we substitue sin^2alpha = 1.0. */

    printf("\n");
    printf("[Derived Quantities] \n");
    printf("I = %12.1e g cm^2 \n",i_ns(m_ns,r_ns));
    printf("E_rot_0 = %12.1e erg \n",ene_rot(m_ns,r_ns,ome_0));
    printf("E_B = %12.1e erg \n",ene_b_tmp);
    printf("epsilon_B = %12.1e \n",eps);
    printf("L_GW_0 = %12.1e erg/s \n",l_gw(m_ns,r_ns,eps,ome_0));
    printf("L_d_0 = %12.1e erg/s \n",l_d(r_ns,b_p,ome_0,0));

    
    int i,j,j_max;
    double initime,timedeg;
    if(timebinmode==0){
      //MKKB15
      j_max = 9;
      initime=exp(6.0*log(10.));
      timedeg=0.25;
    }else if(timebinmode==1){
      //M+19
      j_max = 25;
      initime=exp(3.0*log(10.));
      timedeg=0.25;
    }else if(timebinmode==2){
      //OKM18
      j_max = 61;
      initime=exp(3.0*log(10.));
      timedeg=0.1;
    }else if(timebinmode==3){
      j_max = 31;
      initime=exp(6.0*log(10.));
      timedeg=0.1;
    }else if(timebinmode==4){
      j_max = 41;
      timedeg=log10(tobs/tmax)/(20.-10.);
      initime=tmax/pow(10.,10.*timedeg);
    }else if(timebinmode==5){
      j_max = 81;
      initime=exp(4.0*log(10.));
      timedeg=0.1;
    }
    
    
    if(diskwindmode==2){
      gamma_b=1.; //windのminimum injection
    }
    double e_hardx = 5.0e4*EV_TO_ERG; /* in unit of [erg] */
    double e_gamma = 1.0e9*EV_TO_ERG; /* in unit of [erg] */
    
    k_b_pwn = 0.0; /* 0 or 1 (pdV work) */
    /* Eq. (17) of Murase+15 */
    if(diskwindmode==2){
      eps_mag = 0.003; //new
    }else{
      eps_mag = 0.003;
    }
        
    /* the initial condition */
    double ene_sd = 0.0;
    double ene_nb = 0.0;
    double ene_gw = 0.0;
    double ene_th = 0.5*ene_sn;
    double ene_kin = 0.5*ene_sn;
    double ene_decay = 0.0;
    double ene_Ni = 0.0;
    double ene_Co = 0.0;
    double ene_rpe = 0.0;
    
    double ene_b_pwn = 0.0; /* the magnetic energy in the pwn */
    double u_b_pwn = 0.0;
    double b_pwn = 0.0; /* in unit of [G]*/
    
    double r_ej = r_0;
    double r_w = 100*r_ns; //modified
    double vol_ej_tmp = vol_ej(r_ej);
    double tau_ej_tmp = tau_ej(m_ej,delta,r_ej,kappa);
    double v_ej = pow(2.0*ene_kin/m_ej,0.5);
    double T_ej = pow(ene_th/A_RAD/vol_ej_tmp,0.25);
    double vol_w_tmp = vol_ej(r_w);
    double v_w;
    double dm_ej = 0.0;
   
    double eps_th;
    if(diskwindmode==0){
        eps_e=0.997;//pulsar
        eps_th=eps_e;
    }else if(diskwindmode==1){
        eps_e=1.0;//BH-x
        eps_th=eps_e;
    }else if(diskwindmode==2){
        eps_e=0.1;//BH-IS
        eps_th=1.0; //0=all thermal escape, 1.0=all radiation absorbed
    }
    
    double fac_psr_dep_tmp = 1.0; 
    double fac_Ni_dep_tmp = 1.0; 
    double fac_Co_dep_tmp = 1.0;
    double spec_non_thermal_x = 1.0;
    double f_gamma_esc_x = 1.0;
    double spec_non_thermal_gamma = 1.0;
    double f_gamma_esc_gamma = 1.0;
    

    /* timescales */
    int i_max = 5000; /* resolution */
    double t_0 = r_0/v_ej;
    double t_f;
    if(timebinmode!=5){
      t_f = 3.6525e5*DAY_TO_SEC; /* in unit of [s] */ //注意
    }else{
      t_f = 3.6525e6*DAY_TO_SEC; /* in unit of [s] */ //注意
    }
    double t = t_0;
    double del_ln_t = (log(t_f)-log(t_0))/(double)(i_max-1.0);
    double t_dif_ej = tau_ej_tmp*r_ej/C;
    if(tau_ej_tmp<1){
      t_dif_ej = r_ej/C;
    }
    double t_dyn = r_ej/v_ej;
    double t_U_ab_max=t,t_B_ab_max=t,t_V_ab_max=t,t_R_ab_max=t,t_hardx_max = t,t_gamma_max=t;
    
    v_w=r_w/t_0;
    
    double ome = ome_0;
    double del_ln_ome = -t/ome*(k_gw*pow(ome,5.0)+k_d*pow(ome,3.0))*del_ln_t; /* Eq. (4) of Dall'Osso et al. (2009) */
    double f = ome/M_PI;
    
    /* the thermal SN emission */
    double l_emi_sn = ene_th/t_dif_ej;
    double ene_emi_sn = 0.0;
    double U_ab=0.0,B_ab=0.0,V_ab=0.0,R_ab=0.0,I_ab=0.0;
    double U_app=0.0,B_app=0.0,V_app=0.0,R_app=0.0,I_app=0.0;
    ab_AB_mag(l_emi_sn,T_ej,&U_ab,&B_ab,&V_ab,&R_ab,&I_ab);
    app_AB_mag(l_emi_sn,T_ej,d_l,&U_app,&B_app,&V_app,&R_app,&I_app);
    double U_ab_max=U_ab,B_ab_max=B_ab,V_ab_max=V_ab,R_ab_max=R_ab;
    double dm_15_U=0.0,dm_15_B=0.0,dm_15_V=0.0,dm_15_R=0.0;
    double V_minus_R = V_ab-R_ab;

    /* the bolometric non-thermal emission */
    double l_emi_psr = 0.0;
    double ene_emi_psr = 0.0;

    /* the emission in the hard X-ray band (10-100 keV) */
    double l_emi_hardx = 0.0;
    double l_emi_hardx_max = 0.0;

    /* the emission in the gamma-ray band (~GeV) */
    double l_emi_gamma = 0.0;
    double l_emi_gamma_max = 0.0;

    /* the bolometric Ni, Co emission */
    double l_emi_Ni = ((1.0-fac_Ni_dep_tmp)*ene_Ni+(1.0-fac_Co_dep_tmp)*ene_Co)/t_dif_ej;
    double ene_emi_Ni = 0.0;


    /* for the output */
    FILE *op1,*op2,*op3;
    char file_name[128];
    char dat[] = ".dat";
    op1 = fopen("../op/hydro.dat","w+");
    op2 = fopen("../op/lc.dat","w+");
    op3 = fopen("../op/optical.dat","w+");

    double calctime[j_max];
    for(j=0;j<j_max;j++){
        calctime[j]=initime*exp(timedeg*j*log(10.));
    }


    /* main loop for the time evolution */
    for (i=0;i<i_max;i++) {
      ome += ome*del_ln_ome;
      f = ome/M_PI; /* m = 2 bar mode */
      t = t_0*exp(del_ln_t*(double)i);
      del_ln_ome = -t/ome*(k_gw*pow(ome,5.0)+k_d*pow(ome,3.0))*del_ln_t; /* Eq. (4) of Dall'Osso et al. (2009) */
      
      ene_sd += (l_d(r_ns,b_p,ome,t)+l_disk(t))*t*del_ln_t;
      if (tau_ej_tmp*(r_w/r_ej) > C/v_w){
	ene_nb += (eps_e*l_d(r_ns,b_p,ome,t)+eps_th*l_disk(t))*t*del_ln_t;
      }else{
	ene_nb += (tau_ej_tmp*(r_w/r_ej)*(v_w/C)*eps_e*l_d(r_ns,b_p,ome,t)+eps_th*l_disk(t))*t*del_ln_t;
      }
      ene_gw += l_gw(m_ns,r_ns,eps,ome)*t*del_ln_t;
      ene_kin += ene_th/t_dyn*t*del_ln_t;
      ene_decay += (l_Ni(m_Ni,t)+l_Co(m_Ni,t)+l_rpe(m_rpe,t))*t*del_ln_t;
      
      ene_th += (fac_psr_dep_tmp*(eps_e*l_d(r_ns,b_p,ome,t)+eps_e*l_disk(t))+(eps_th-eps_e)*l_disk(t)
		 +fac_Ni_dep_tmp*l_Ni(m_Ni,t)
		 +fac_Co_dep_tmp*l_Co(m_Ni,t)
		 +l_rpe(m_rpe,t)
		 -l_emi_sn
		 -ene_th/t_dyn)*t*del_ln_t; //disk-outflow reverse shock thermal emission無視
      if (ene_th < 0.0)
	ene_th = 0.0;
      
      ene_emi_sn += l_emi_sn*t*del_ln_t;
      ene_emi_psr += (1.0-fac_psr_dep_tmp)*eps_e*(l_d(r_ns,b_p,ome,t)+l_disk(t))*t*del_ln_t;
      ene_emi_Ni += ((1.0-fac_Ni_dep_tmp)*l_Ni(m_Ni,t)+(1.0-fac_Co_dep_tmp)*l_Co(m_Ni,t))*t*del_ln_t;     
      ene_b_pwn += (eps_mag*l_d(r_ns,b_p,ome,t) - k_b_pwn*ene_b_pwn/t_dyn)*t*del_ln_t; /* Previously, the second term is neglected */

      if(diskwindmode==2){
	u_b_pwn = eps_mag*3.0*(0.02*1.99*1e33*(0.3*C*0.3*C)/2.)/4.0/M_PI/r_w/r_w/r_w; //0.003 is equivalent to 0.003 for M_w=0.02
	//u_b_pwn = eps_mag*0.3*1e-3*1.99*1e33*(0.3*C*0.3*C)/pow((r_w/0.3/C),5./3.)/4.0/M_PI/r_w/r_w/0.3/C; //new??
      }else{
	u_b_pwn = 3.0*ene_b_pwn/4.0/M_PI/r_w/r_w/r_w; //eps_magいらない
      }
      b_pwn = pow(u_b_pwn*8.0*M_PI,0.5); //be careful: epsilon_B=epsilon_B^pre/8/PI for epsilon_B^pre used in Murase et al. 2018
             
      v_ej = pow(2.0*ene_kin/(m_ej+dm_ej),0.5); //including the deceleartion effect
      if(sqrt(7.0/6.0/(3.0-delta)*ene_nb/m_ej*pow(r_ej/r_w,3.0-delta))>pow((2*l_disk(t)*r_w/(3-delta)/m_ej),1./3.)){
	v_w=(sqrt(7.0/6.0/(3.0-delta)*ene_nb/m_ej*pow(r_ej/r_w,3.0-delta))+r_w/t);//eps_magむし
      }else{
	v_w=pow((2*l_disk(t)*r_w/(3-delta)/m_ej),1./3.)+r_w/t;//eps_magむし
      }
      vol_ej_tmp = vol_ej(r_ej);
      tau_ej_tmp = tau_ej(m_ej,delta,r_ej,kappa);
      t_dif_ej = tau_ej_tmp*r_ej/C;
      if(tau_ej_tmp<1){
	t_dif_ej = r_ej/C;
      }
      t_dyn = r_ej/v_ej;
      T_ej = pow(ene_th/A_RAD/vol_ej_tmp,0.25);
      
      
      /* UV/optiacl/IR */
      l_emi_sn = ene_th/t_dif_ej;
      ab_AB_mag(l_emi_sn,T_ej,&U_ab,&B_ab,&V_ab,&R_ab,&I_ab);
      app_AB_mag(l_emi_sn,T_ej,d_l,&U_app,&B_app,&V_app,&R_app,&I_app);
      if (R_ab < R_ab_max){
	R_ab_max = R_ab;
	t_R_ab_max = t;
      }
      if (V_ab < V_ab_max){
	V_ab_max = V_ab;
	t_V_ab_max = t;
      }
      if (U_ab < U_ab_max){
	U_ab_max = U_ab;
	t_U_ab_max = t;
      }
      if (t-t_R_ab_max<15.0*DAY_TO_SEC){
	dm_15_R = R_ab-R_ab_max;
      }
      if (t-t_V_ab_max<15.0*DAY_TO_SEC){
	dm_15_V = V_ab-V_ab_max;
      }
      V_minus_R = V_ab-R_ab;
      

      /* Gamma-ray and X-ray */
      fac_psr_dep_tmp = fac_psr_dep(m_ej,delta,r_ej,albd_fac,b_pwn,gamma_b,T_ej);
      fac_Ni_dep_tmp = fac_Ni_dep(m_ej,delta,r_ej,albd_fac);
      fac_Co_dep_tmp = fac_Co_dep(m_ej,delta,r_ej,albd_fac);

      spec_non_thermal_x = spec_non_thermal(e_hardx,b_pwn,gamma_b,T_ej);
      f_gamma_esc_x = f_gamma_esc(e_hardx,m_ej,delta,r_ej,albd_fac);
      spec_non_thermal_gamma = spec_non_thermal(e_gamma,b_pwn,gamma_b,T_ej);
      f_gamma_esc_gamma = f_gamma_esc(e_gamma,m_ej,delta,r_ej,albd_fac);
            
      l_emi_psr = (1.0-fac_psr_dep_tmp)*eps_e*(l_d(r_ns,b_p,ome,t)+l_disk(t));
      l_emi_hardx = eps_e*(l_d(r_ns,b_p,ome,t)+l_disk(t))*spec_non_thermal_x*f_gamma_esc_x*e_hardx;
      l_emi_gamma = eps_e*(l_d(r_ns,b_p,ome,t)+l_disk(t))*spec_non_thermal_gamma*f_gamma_esc_gamma*e_gamma;
      l_emi_Ni = ((1.0-fac_Ni_dep_tmp)*ene_Ni+(1.0-fac_Co_dep_tmp)*ene_Co)/t_dif_ej;

      if (l_emi_hardx > l_emi_hardx_max){
	l_emi_hardx_max = l_emi_hardx;
	t_hardx_max = t;
      }
      
      if (l_emi_gamma > l_emi_gamma_max){
	l_emi_gamma_max = l_emi_gamma;
	t_gamma_max = t;
      }            
      

      /* output */
      for(j=0;j<j_max;j++){
	if((t<=calctime[j])&&(calctime[j])<t*exp(del_ln_t)){
	  fprintf(op1,"%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e \n",
		  t,ome,r_ej,v_ej,r_w,v_w,tau_ej_tmp,T_ej,b_pwn,ene_th,ene_b_pwn);
	  fprintf(op2,"%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e \n",
		  t,l_emi_sn,l_gw(m_ns,r_ns,eps,ome),l_d(r_ns,b_p,ome,t),l_disk(t),l_Ni(m_Ni,t)+l_Co(m_Ni,t),l_rpe(m_rpe,t),l_emi_hardx,l_emi_gamma,l_emi_psr,l_emi_Ni,
		  spec_non_thermal_x,f_gamma_esc_x,spec_non_thermal_gamma,f_gamma_esc_gamma);
	  fprintf(op3,"%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e \n",
		  t,U_ab,B_ab,V_ab,R_ab,I_ab,U_app,B_app,V_app,R_app,I_app,V_minus_R);
	}
      }
      
      
      /* evolve radii and so on */
      r_ej += v_ej*t*del_ln_t;
      r_w += v_w*t*del_ln_t;
      if (r_w > r_ej){
	r_w = r_ej;
      }
      vol_w_tmp = vol_ej(r_w);
      dm_ej += 4.0*M_PI*r_ej*r_ej*v_ej*t*del_ln_t*rho_env(r_ej);
      
    }
    fclose(op1);
    fclose(op2);
    fclose(op3);


    /* logging */
    printf("\n");
    printf("P_f = %12.1e s \n",2.0*M_PI/ome);
    printf("E_rot_f = %12.1e erg \n",ene_rot(m_ns,r_ns,ome));
    printf("E_sd_f+E_disk_f = %12.1e erg \n",ene_sd);
    printf("E_GW_f = %12.1e erg \n",ene_gw);
    printf("E_emi_f = %12.1e erg \n",ene_emi_sn+ene_emi_psr);
    printf("E_int_f = %12.1e erg \n",ene_th);
    printf("E_kin_f = %12.1e erg \n",ene_kin);
    
    printf("\n");
    printf("R_ab_max = %12.3e \n",R_ab_max);
    printf("V_ab_max = %12.3e \n",V_ab_max);
    printf("t_R_ab_max = %12.3e \n",t_R_ab_max/DAY_TO_SEC);
    printf("t_V_ab_max = %12.3e \n",t_V_ab_max/DAY_TO_SEC);
    printf("dm_15_R = %12.3e \n",dm_15_R);
    printf("dm_15_V = %12.3e \n",dm_15_V);
    printf("L_hardx_max = %12.3e \n",l_emi_hardx_max);
    printf("L_gamma_max = %12.3e \n",l_emi_gamma_max);
    printf("t_hardx_max = %12.3e \n",t_hardx_max/DAY_TO_SEC);
    printf("t_gamma_max = %12.3e \n",t_gamma_max/DAY_TO_SEC);    

    printf("\n");
    printf("[Conservation Check] \n");
    printf("(E_GW+E_sd+E_disk+E_rot)/E_rot_0 = %12.3e \n",(ene_gw+ene_sd+ene_rot(m_ns,r_ns,ome))/ene_rot(m_ns,r_ns,ome_0));
    printf("(E_emi+E_int+E_kin)/(E_sd+E_disk+E_SN+E_Ni) = %12.3e \n",(ene_emi_sn+ene_emi_psr+ene_emi_Ni+ene_th+ene_kin)/(ene_sd+ene_sn+ene_decay));
    printf("\n");
    
}

double s_f_num(double f)
{
    /* The offical adv.LIGO sensitivity curve */
    int i,k;
    int n_max = 10; /* data size */
    FILE *op;
    double f_data[n_max],s_f_data[n_max];
    
    
    op = fopen("ZERO_DET_high_P.txt","r"); /* read the fitting data file */
    
    for (i=0;i<n_max;i++) {
        fscanf(op,"%le %le\n",&f_data[i],&s_f_data[i]);
    }
    
    fclose(op);
    
    i=0;
    if(f < f_data[0] || f > f_data[n_max-1]){
        return INFINITY;
    } else
        while(i<n_max){
            if(f <= f_data[i+1] && f_data[i] <= f){
                return exp(2.0*(log(s_f_data[i])+(log(s_f_data[i+1])-log(s_f_data[i]))/(log(f_data[i+1])-log(f_data[i]))*(log(f)-log(f_data[i]))));
            }else i++;
        }
    
    return 0.0;
}

/* Conversion from Bolometric luminosity [erg/s] and temperature [K] to absolute AB magnitude */
void ab_AB_mag(double L_bol, double T_ph, double *U_ab, double *B_ab, double *V_ab, double *R_ab, double *I_ab)
{
    double d_L = 10.0*PC_TO_CM; /* calibrated at 10 pc */
    
    double U_ab_tmp;
    double lam_c_U = 0.3652e-4; /* central wavelengh of U band */
    double F_0_U = 3.630e-20; /* normalization flux of U band */
    double x_U = H*C/lam_c_U/K_B/T_ph;
    
    double B_ab_tmp;
    double lam_c_B = 0.4448e-4; /* central wavelengh of B band */
    double F_0_B = 3.630e-20; /* normalization flux of B band */
    double x_B = H*C/lam_c_B/K_B/T_ph;
    
    double V_ab_tmp;
    double lam_c_V = 0.5505e-4; /* central wavelengh of V band */
    double F_0_V = 3.630e-20; /* normalization flux of V band */
    double x_V = H*C/lam_c_V/K_B/T_ph;
    
    double R_ab_tmp;
    double lam_c_R = 0.6588e-4; /* central wavelengh of R band */
    double F_0_R = 3.630e-20; /* normalization flux of R band */
    //double F_0_R = 3.080e-20; /* normalization flux of R band */
    double x_R = H*C/lam_c_R/K_B/T_ph;
    
    double I_ab_tmp;
    double lam_c_I = 0.8060e-4; /* central wavelengh of I band */
    double F_0_I = 3.630e-20; /* normalization flux of I band */
    double x_I = H*C/lam_c_I/K_B/T_ph;
    
    U_ab_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_U,3.0)/pow(C,3.0)/(exp(x_U)-1.0)/(d_L*d_L)/F_0_U);
    B_ab_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_B,3.0)/pow(C,3.0)/(exp(x_B)-1.0)/(d_L*d_L)/F_0_B);
    V_ab_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_V,3.0)/pow(C,3.0)/(exp(x_V)-1.0)/(d_L*d_L)/F_0_V);
    R_ab_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_R,3.0)/pow(C,3.0)/(exp(x_R)-1.0)/(d_L*d_L)/F_0_R);
    I_ab_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_I,3.0)/pow(C,3.0)/(exp(x_I)-1.0)/(d_L*d_L)/F_0_I);
    //R_ab_tmp = -2.5*log10(0.75*L_bol/(C/lam_c_R)/(4.0*M_PI*d_L*d_L)/F_0_R);
    /* the factor 0.75 corresponds to the fraction in optical band around the peak */
    /* See Fig. 8 of Stritzinger 2009. */
    
    *U_ab = U_ab_tmp;
    *B_ab = B_ab_tmp;
    *V_ab = V_ab_tmp;
    *R_ab = R_ab_tmp;
    *I_ab = I_ab_tmp;
}

/* Conversion from Bolometric luminosity [erg/s] to apparent magnitude */
void app_AB_mag(double L_bol, double T_ph, double d_L, double *U_app, double *B_app, double *V_app, double *R_app, double *I_app)
{
    double U_app_tmp;
    double lam_c_U = 0.3652e-4; /* central wavelengh of U band */
    double F_0_U = 3.630e-20; /* normalization flux of U band */
    double x_U = H*C/lam_c_U/K_B/T_ph;
    
    double B_app_tmp;
    double lam_c_B = 0.4448e-4; /* central wavelengh of B band */
    double F_0_B = 3.630e-20; /* normalization flux of B band */
    double x_B = H*C/lam_c_B/K_B/T_ph;
    
    double V_app_tmp;
    double lam_c_V = 0.5505e-4; /* central wavelengh of V band */
    double F_0_V = 3.630e-20; /* normalization flux of V band */
    double x_V = H*C/lam_c_V/K_B/T_ph;
    
    double R_app_tmp;
    double lam_c_R = 0.6588e-4; /* central wavelengh of R band */
    double F_0_R = 3.630e-20; /* normalization flux of R band */
    double x_R = H*C/lam_c_R/K_B/T_ph;
    
    double I_app_tmp;
    double lam_c_I =  2.179e-4;//K-band //I-band=0.8060e-4; /* central wavelengh of R band */
    double F_0_I = 3.630e-20; /* normalization flux of R band */
    double x_I = H*C/lam_c_I/K_B/T_ph;
    
    U_app_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_U,3.0)/pow(C,3.0)/(exp(x_U)-1.0)/(d_L*d_L)/F_0_U);
    B_app_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_B,3.0)/pow(C,3.0)/(exp(x_B)-1.0)/(d_L*d_L)/F_0_B);
    V_app_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_V,3.0)/pow(C,3.0)/(exp(x_V)-1.0)/(d_L*d_L)/F_0_V);
    R_app_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_R,3.0)/pow(C,3.0)/(exp(x_R)-1.0)/(d_L*d_L)/F_0_R);
    I_app_tmp = -2.5*log10(L_bol/(A_RAD*pow(T_ph,4.0))*2.0*H*pow(C/lam_c_I,3.0)/pow(C,3.0)/(exp(x_I)-1.0)/(d_L*d_L)/F_0_I);
    
    *U_app = U_app_tmp;
    *B_app = B_app_tmp;
    *V_app = V_app_tmp;
    *R_app = R_app_tmp;
    *I_app = I_app_tmp;
}

/* absolute bolometric magnitude */
double ab_bol_mag(double L_bol)
{
    double ab_bol_mag_cari = 4.8;
    double L_bol_cari = 4e33;

    return ab_bol_mag_cari-2.5*log10(L_bol/L_bol_cari);
}

/* Energy injection rate by Ni decay */
double l_Ni(double m_Ni, double t)
{
    double tau_Ni = 8.8*DAY_TO_SEC;
    double epsi_Ni = 3.9e10; /* in unit of [erg/s/g] */
    
    return m_Ni*epsi_Ni*exp(-t/tau_Ni);
}

/* Energy injection rate by Co decay */
double l_Co(double m_Ni, double t)
{
    double tau_Ni = 8.8*DAY_TO_SEC;
    double tau_Co = 111.3*DAY_TO_SEC;
    double epsi_Co = 6.8e9; /* in unit of [erg/s/g] */
    
    return m_Ni*epsi_Co*(exp(-t/tau_Co)-exp(-t/tau_Ni));
    
}

/* Energy injection rate by r-process */
double l_rpe(double m_rpe, double t)
{
    double meancharge;
    double meanA;
   
    if(opacitymode==0){
    meancharge=26;
    meanA=55.845;
    }else if(opacitymode==1){
    meancharge=30;
    meanA=(55.845+78.971)/2.;
    }else if(opacitymode==2){
    meancharge=(54+79)/2.;
    meanA=(131.293+196.967)/2.;
    }else if(opacitymode==3){
    meancharge=7.0;
    meanA=14.0;
    }
    
    double A0,A1,A2,A3;
    A0 = 15.8048255330964;
    A1 = -1.20355385508098;
    A2 = -0.0223153602207339;
    A3 = 0.00123960875933877;
 
    return 5.*m_rpe*pow(10.,A0+A1*log10(t)+A2*log10(t)*log10(t)+A3*log10(t)*log10(t)*log10(t))/pow((meancharge/70.),1/3.)/(meanA/200); //2-5=possible enhancement in gammas
    
}

/* energy depostision fraction from gamma rays to the thermal bath by Compton (Eq. 46 of Murase+15) */
double gamma_ene_depo_frac_Compton(double e_gamma)
{
    double x = e_gamma/M_ELEC/C/C;
    
    if (x > 1.0e-3)
        return 3.0/4.0*SIGMA_T*(2.0*pow(1.0+x,2.0)/x/x/(1.0+2.0*x)
                            -(1.0+3.0*x)/pow(1.0+2.0*x,2.0)
                            -(1.0+x)*(2.0*x*x-2.0*x-1.0)/x/x/pow(1.0+2.0*x,2.0)
                            -4*x*x/3.0/pow(1.0+2.0*x,3.0)
                            -((1.0+x)/x/x/x-1.0/2.0/x+1.0/2.0/x/x/x)*log(1.0+2.0*x));
    else
        return SIGMA_T*x;
    
}

/* The total KN cross section in unit of cm^2 */
double sigma_kn(double e_gamma)
{
    double x = e_gamma/M_ELEC/C/C;

    if (x > 1.0e-3)
        return 3.0/4.0*SIGMA_T*((1.0+x)/x/x/x*(2.0*x*(1.0+x)/(1.0+2.0*x)-log(1.0+2.0*x))+1.0/2.0/x*log(1.0+2.0*x)-(1.0+3.0*x)/pow(1.0+2.0*x,2.0));
    else
        return SIGMA_T;
}


/* The total BH cross section in unit of cm^2 */
double sigma_BH_p(double e_gamma)
{
    /* Eq. (A1) and (A2) of Chodorowski+92 */
    double x = e_gamma/M_ELEC/C/C;
    double log2x = log(2.0*x);
    double eta = (x-2.0)/(x+2.0);
    double alpha = 1.0/137.0;
    double zeta3 = 1.2020569;
    
    if (x > 4.0)
        return 3.0*alpha/8.0/M_PI*SIGMA_T*(28.0/9.0*log2x-218.0/27.0
                                           +pow(2.0/x,2.0)*(6.0*log2x-7.0/2.0+2.0/3.0*pow(log2x,3.0)-pow(log2x,2.0)-1.0/3.0*M_PI*M_PI*log2x+2.0*zeta3+M_PI*M_PI/6.0)
                                           -pow(2.0/x,4.0)*(3.0/16.0*log2x+1.0/8.0)
                                           -pow(2.0/x,6.0)*(29.0/9.0/256.0*log2x-77.0/27.0/512.0));
    else if (x > 2.0)
        return 1.0/4.0*alpha*SIGMA_T*pow((x-2.0)/x,3.0)*(1.0+eta/2.0+23.0/40.0*eta*eta+37.0/120.0*pow(eta,3.0)+61.0/192*pow(eta,4.0));
    else
        return 0.0;
}


/* inelastisity of gamma rays in Compton */
double gamma_inelas_Compton(double e_gamma)
{
    double x = e_gamma/M_ELEC/C/C;
    
    return gamma_ene_depo_frac_Compton(e_gamma)/sigma_kn(e_gamma);
}

/* inelastisity of gamma rays in Compton */
double gamma_inelas_BH(double e_gamma)
{
    /* Eq. (48) of Murase+15 */
    double x = e_gamma/M_ELEC/C/C;
        
    if (x > 2.0)
        return (x-2.0)/x;
    else
        return 0.0;
}

/* opacity of boud-free emission */
double kappa_bf(double e_gamma, double Z_eff)
{
    double zeta = 1.0;//0.5 for OMK18, 1.0 for Murase+18; /* neutral fraction */
    
    /* See http://physics.nist.gov/PhysRefData/XrayMassCoef/tab3.html */
    //return 5.0*zeta*pow(e_gamma/EV_TO_ERG/1.0e4,-3.0)*pow(Z_eff/7.0,3.0);
    
    double ironopacity,seleopacity,xenonopacity,goldopacity;
    
    double A0,A1,A2,A3,A4;
    //iron
    if(log10(e_gamma/EV_TO_ERG/1.0e6)<-2.14801){
        A0 = -3.95919414261072;
        A1 = -2.64892215754265;
    }else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-0.695){
        A0 = -3.37291030805215;
        A1 = -2.75161208271434;
    }else{
        A0 = -1.59069232728045;
        A1 = -0.206813265289848;
    }
    ironopacity=zeta*pow(10.,A0+A1*log10(e_gamma/EV_TO_ERG/1.0e6));
    //selenium
    if(log10(e_gamma/EV_TO_ERG/1.0e6)<-2.84291){
        A0 = -3.78835348616654;
        A1 = -2.38423803432305;
    }else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-1.89764){
        A0 = -3.68604734441612;
        A1 = -2.66063041649055;
    }else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-0.55){
        A0 = -2.93083141712927;
        A1 = -2.60630263958148;
    }else{
        A0 = -1.58243193342158;
        A1 = -0.102384218895718;
    }
    seleopacity=zeta*pow(10.,A0+A1*log10(e_gamma/EV_TO_ERG/1.0e6));
    //xenon
    if(log10(e_gamma/EV_TO_ERG/1.0e6)<-2.32037){
        A0 = -3.07458863553159;
        A1 = -2.35739410975398;
    }else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-1.46141){
        A0 = -3.17731357386225;
        A1 = -2.68342346938979;
    }else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-0.282){
        A0 = -2.17345283895274;
        A1 = -2.26742402391864;
    }else{
        A0 = -1.55866825608716;
        A1 = -0.0467127630289143;
    }
    xenonopacity=zeta*pow(10.,A0+A1*log10(e_gamma/EV_TO_ERG/1.0e6));
    //gold
    if(log10(e_gamma/EV_TO_ERG/1.0e6)<-2.65645){
        A0 = -2.5113444206149;
        A1 = -2.06076550942316;
    }else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-1.92377){
        A0 = -2.45512766933321;
        A1 = -2.27191147638504;
    }else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-1.09299){
        A0 = -2.32582461907071;
        A1 = -2.39063081204358;
    }else if(log10(e_gamma/EV_TO_ERG/1.0e6)<0.02){
        A0 = -1.55199838865465;
        A1 = -1.82076527957878;
    }else{
        A0 = -1.58507209319691;
        A1 = 0.0628004018301846;
    }
    goldopacity=zeta*pow(10.,A0+A1*log10(e_gamma/EV_TO_ERG/1.0e6));
    
    if(opacitymode==0){
    return ironopacity;
    }else if(opacitymode==1){
    return 0.5*ironopacity+0.5*seleopacity;
    }else if(opacitymode==2){
    return 0.5*xenonopacity+0.5*goldopacity;
    }else if(opacitymode==3){
    return 5.0*zeta*pow(e_gamma/EV_TO_ERG/1.0e4,-3.0)*pow(Z_eff/7.0,3.0); //approximate formula used in Murase et al. 2015
    }
    
}


/* escape fraction of gamma rays interms of energy */
double f_gamma_dep(double e_gamma, double m_ej, double delta, double r_ej, double albd_fac)
{
    double mu_e; /* electron mean molecular weight */
    //double Z_eff = 7.0; /* effective nuclear weight */
    /* this corresponds to C:O = 1:1 */
    double Z_eff;
    if(opacitymode==0){
    Z_eff = 24.21; /* effective nuclear weight */
    mu_e = 2.148;
    }else if(opacitymode==1){
    Z_eff = 26.74;
    mu_e = 2.2353;
    }else if(opacitymode==2){
    Z_eff = 53.90;
    mu_e= 2.4622;
    }else if(opacitymode==3){
    Z_eff = 2.0;
    mu_e= 7.0;
    }
    
    double tau_Compton = (3.0-delta)/4.0/M_PI*m_ej*sigma_kn(e_gamma)/mu_e/M_PRO/r_ej/r_ej;
    double tau_BH = (3.0-delta)/4.0/mu_e/M_PI*m_ej*(1.0+Z_eff)*sigma_BH_p(e_gamma)/M_PRO/r_ej/r_ej;
    double tau_bf = (1.0-albd_fac)*(3.0-delta)/4.0/M_PI*m_ej*kappa_bf(e_gamma,Z_eff)/r_ej/r_ej;
    
    double power_Compton=0.0;
    if (tau_Compton > 1.0)
        power_Compton = tau_Compton*tau_Compton;
    else
        power_Compton = tau_Compton;
    
    double f_gamma_dep_tmp = (1.0-pow(1.0-gamma_inelas_Compton(e_gamma),power_Compton))+(1.0-exp(-(tau_BH+tau_bf)));
        
    if (f_gamma_dep_tmp > 1.0)
        return 1.0;
    else
        return f_gamma_dep_tmp;
}

double f_gamma_esc(double e_gamma, double m_ej, double delta, double r_ej, double albd_fac)
{
    double mu_e; /* electron mean molecular weight */
    //double Z_eff = 7.0; /* effective nuclear weight */
    /* this corresponds to C:O = 1:1 */
    double Z_eff;
    if(opacitymode==0){
    Z_eff = 24.21; /* effective nuclear weight */
    mu_e = 2.148;
    }else if(opacitymode==1){
    Z_eff = 26.74;
    mu_e = 2.2353;
    }else if(opacitymode==2){
    Z_eff = 53.90;
    mu_e= 2.4622;
    }else if(opacitymode==3){
    Z_eff = 7.0;
    mu_e = 2.0;
    }
    
    double tau_Compton = (3.0-delta)/4.0/M_PI*m_ej*sigma_kn(e_gamma)/mu_e/M_PRO/r_ej/r_ej;
    double tau_BH = (3.0-delta)/4.0/mu_e/M_PI*m_ej*(1.0+Z_eff)*sigma_BH_p(e_gamma)/M_PRO/r_ej/r_ej;
    //double tau_bf = 2.0*(1.0-albd_fac)*(3.0-delta)/4.0/M_PI*m_ej*kappa_bf(e_gamma,Z_eff)/r_ej/r_ej;
    double tau_bf = (1.0-albd_fac)*(3.0-delta)/4.0/M_PI*m_ej*kappa_bf(e_gamma,Z_eff)/r_ej/r_ej;
    
    double tau_abs = (1.0+gamma_inelas_Compton(e_gamma))*(tau_BH+tau_bf);
    double tau_eff = sqrt((tau_abs+tau_Compton)*tau_abs);
    
    double power_Compton=0.0;
    if (tau_Compton > 1.0)
        power_Compton = tau_Compton*tau_Compton;
    else
        power_Compton = tau_Compton;
    
    
    //return exp(-tau_eff);
    return exp(-(tau_BH+tau_bf))*(exp(-(tau_Compton))+(1.0-exp(-(tau_Compton)))*pow(1.0-gamma_inelas_Compton(e_gamma),power_Compton));
    //return (exp(-(tau_Compton))+(1.0-exp(-(tau_Compton)))*pow(1.0-gamma_inelas_Compton(e_gamma),power_Compton));

}

/* escacpe fraction of Ni-decay energy */
double fac_Ni_dep(double m_ej, double delta, double r_ej, double albd_fac)
{    
    /* see table 1 of Nadyozhin 1994 */
    double e_Ni_1 = 1.56e6*EV_TO_ERG;
    double e_Ni_2 = 0.812e6*EV_TO_ERG;
    double e_Ni_3 = 0.75e6*EV_TO_ERG;
    double e_Ni_4 = 0.48e6*EV_TO_ERG;
    double e_Ni_5 = 0.27e6*EV_TO_ERG;
    double e_Ni_6 = 0.158e6*EV_TO_ERG;
    
    double prob_Ni_1 = 0.14;
    double prob_Ni_2 = 0.86;
    double prob_Ni_3 = 0.495;
    double prob_Ni_4 = 0.366;
    double prob_Ni_5 = 0.365;
    double prob_Ni_6 = 1.0;
    
    double e_Ni_tot = e_Ni_1*prob_Ni_1+e_Ni_2*prob_Ni_2+e_Ni_3*prob_Ni_3+e_Ni_4*prob_Ni_4+e_Ni_5*prob_Ni_5+e_Ni_6*prob_Ni_6;
    
    
    return (prob_Ni_1*e_Ni_1*f_gamma_dep(e_Ni_1,m_ej,delta,r_ej,albd_fac)
            +prob_Ni_2*e_Ni_2*f_gamma_dep(e_Ni_2,m_ej,delta,r_ej,albd_fac)
            +prob_Ni_3*e_Ni_3*f_gamma_dep(e_Ni_3,m_ej,delta,r_ej,albd_fac)
            +prob_Ni_4*e_Ni_4*f_gamma_dep(e_Ni_4,m_ej,delta,r_ej,albd_fac)
            +prob_Ni_5*e_Ni_5*f_gamma_dep(e_Ni_5,m_ej,delta,r_ej,albd_fac)
            +prob_Ni_6*e_Ni_6*f_gamma_dep(e_Ni_6,m_ej,delta,r_ej,albd_fac))/e_Ni_tot;
    
    /*
    double kappa_gamma = 0.03;
    double tau_eff = (3.0-delta)*kappa_gamma*m_ej/4.0/M_PI/(r_ej*r_ej);
    
    return 1.0-exp(-tau_eff);
    */
}

/* escape fraction of Co-decay energy */
double fac_Co_dep(double m_ej, double delta, double r_ej, double albd_fac)
{    
    /* see table 1 of Nadyozhin 1994 */
    double e_Co_1 = 3.24e6*EV_TO_ERG;
    double e_Co_2 = 2.60e6*EV_TO_ERG;
    double e_Co_3 = 2.03e6*EV_TO_ERG;
    double e_Co_4 = 1.77e6*EV_TO_ERG;
    double e_Co_5 = 1.36e6*EV_TO_ERG;
    double e_Co_6 = 1.24e6*EV_TO_ERG;
    double e_Co_7 = 1.04e6*EV_TO_ERG;
    double e_Co_8 = 0.990e6*EV_TO_ERG;
    double e_Co_9 = 0.847e6*EV_TO_ERG;
    double e_Co_10 = 0.511e6*EV_TO_ERG;
    double e_Co_eleposi = 0.632e6*EV_TO_ERG;
    
    double prob_Co_1 = 0.125;
    double prob_Co_2 = 0.17;
    double prob_Co_3 = 0.12;
    double prob_Co_4 = 0.16;
    double prob_Co_5 = 0.043;
    double prob_Co_6 = 0.68;
    double prob_Co_7 = 0.14;
    double prob_Co_8 = 0.028;
    double prob_Co_9 = 1.0;
    double prob_Co_10 = 0.38;
    double prob_Co_eleposi = 0.18;
    
    double e_Co_tot = e_Co_1*prob_Co_1+e_Co_2*prob_Co_2+e_Co_3*prob_Co_3+e_Co_4*prob_Co_4+e_Co_5*prob_Co_5+e_Co_6*prob_Co_6+e_Co_7*prob_Co_7+e_Co_8*prob_Co_8+e_Co_9*prob_Co_9+e_Co_10*prob_Co_10+e_Co_eleposi*prob_Co_eleposi;
    
    
    return (prob_Co_1*e_Co_1*f_gamma_dep(e_Co_1,m_ej,delta,r_ej,albd_fac)
            +prob_Co_2*e_Co_2*f_gamma_dep(e_Co_2,m_ej,delta,r_ej,albd_fac)
            +prob_Co_3*e_Co_3*f_gamma_dep(e_Co_3,m_ej,delta,r_ej,albd_fac)
            +prob_Co_4*e_Co_4*f_gamma_dep(e_Co_4,m_ej,delta,r_ej,albd_fac)
            +prob_Co_5*e_Co_5*f_gamma_dep(e_Co_5,m_ej,delta,r_ej,albd_fac)
            +prob_Co_6*e_Co_6*f_gamma_dep(e_Co_6,m_ej,delta,r_ej,albd_fac)
            +prob_Co_7*e_Co_7*f_gamma_dep(e_Co_7,m_ej,delta,r_ej,albd_fac)
            +prob_Co_8*e_Co_8*f_gamma_dep(e_Co_8,m_ej,delta,r_ej,albd_fac)
            +prob_Co_9*e_Co_9*f_gamma_dep(e_Co_9,m_ej,delta,r_ej,albd_fac)
            +prob_Co_10*e_Co_10*f_gamma_dep(e_Co_10,m_ej,delta,r_ej,albd_fac)
            +prob_Co_eleposi*e_Co_eleposi
            )/e_Co_tot;
    
    /*
    double kappa_gamma = 0.03;
    double tau_eff = (3.0-delta)*kappa_gamma*m_ej/4.0/M_PI/(r_ej*r_ej);
    
    return 1.0-exp(-tau_eff);//+e_Co_eleposi/e_Co_tot;
    */
}

/* escaoe fraction of pulsar emission */
double fac_psr_dep(double m_ej, double delta, double r_ej, double albd_fac, double b_pwn, double gamma_b, double T_ej)
{
    double e_gamma_min;
    double e_gamma_max_tmp;
    double e_gamma_gamma_ani_tmp;
    double e_gamma_syn_b_tmp;
    if(diskwindmode==1){
    e_gamma_min = 10.0*EV_TO_ERG;
    e_gamma_max_tmp = 1e5;
    e_gamma_gamma_ani_tmp = e_gamma_gamma_ani(T_ej);
    e_gamma_syn_b_tmp = 1e3;
    }else{
    e_gamma_min = 1.0*EV_TO_ERG; /* eV */
    e_gamma_max_tmp = e_gamma_max(b_pwn);
    e_gamma_gamma_ani_tmp = e_gamma_gamma_ani(T_ej);
    e_gamma_syn_b_tmp = e_gamma_syn_b(b_pwn,gamma_b);
    }
    
    if (e_gamma_gamma_ani_tmp < e_gamma_max_tmp)
        e_gamma_max_tmp = e_gamma_gamma_ani_tmp;
    
    int i;
    int i_max = 1000; /* resolution */
    
    double e_tmp = e_gamma_min;
    double del_ln_e = 0.0; 
    double frac_psr_dep_tmp = 0.0;
    
    if (e_gamma_max_tmp > e_gamma_syn_b_tmp){
        del_ln_e = (log(e_gamma_max_tmp)-log(e_gamma_min))/(double)(i_max+1);
        for (i=0;i<=i_max;i++) {
            e_tmp = e_gamma_min*exp(del_ln_e*(double)i);
            if (i == 0 || i==i_max){
                frac_psr_dep_tmp += (1.0/2.0)*f_gamma_dep(e_tmp,m_ej,delta,r_ej,albd_fac)*spec_non_thermal(e_tmp,b_pwn,gamma_b,T_ej)*e_tmp*del_ln_e;
            } else if (i % 2 == 0) {
                frac_psr_dep_tmp += (2.0/3.0)*f_gamma_dep(e_tmp,m_ej,delta,r_ej,albd_fac)*spec_non_thermal(e_tmp,b_pwn,gamma_b,T_ej)*e_tmp*del_ln_e;
            } else {
                frac_psr_dep_tmp += (4.0/3.0)*f_gamma_dep(e_tmp,m_ej,delta,r_ej,albd_fac)*spec_non_thermal(e_tmp,b_pwn,gamma_b,T_ej)*e_tmp*del_ln_e;
            }
        }
    } else {
        e_gamma_max_tmp = e_gamma_syn_b_tmp;
        del_ln_e = (log(e_gamma_max_tmp)-log(e_gamma_min))/(double)(i_max+1);
        for (i=0;i<=i_max;i++) {
            e_tmp = e_gamma_min*exp(del_ln_e*(double)i);
            if (i == 0 || i==i_max){
                frac_psr_dep_tmp += (1.0/3.0)*f_gamma_dep(e_tmp,m_ej,delta,r_ej,albd_fac)*spec_non_thermal(e_tmp,b_pwn,gamma_b,T_ej)*e_tmp*del_ln_e;
            } else if (i % 2 == 0) {
                frac_psr_dep_tmp += (2.0/3.0)*f_gamma_dep(e_tmp,m_ej,delta,r_ej,albd_fac)*spec_non_thermal(e_tmp,b_pwn,gamma_b,T_ej)*e_tmp*del_ln_e;
            } else {
                frac_psr_dep_tmp += (4.0/3.0)*f_gamma_dep(e_tmp,m_ej,delta,r_ej,albd_fac)*spec_non_thermal(e_tmp,b_pwn,gamma_b,T_ej)*e_tmp*del_ln_e;
            }
        }

    }

    return frac_psr_dep_tmp;
}

double l_disk(double time)
{
    if(diskwindmode==0){
    //no fall back
    return 0;
    }else if(diskwindmode==1){
    //BH model - disk
    double ulxlum=1.0*1e40;//15.*2*2.8*1.26*1e38;
    //double acclum=0.1*(2./3.)*(0.01/2.)*1.99*1e33*(C*C)/pow((time/2.),5./3.);
    double acclum=0.1*1e-3*1.99*1e33*(C*C)/pow((time/1.),5./3.);
    if(ulxlum<acclum){
        return ulxlum;
    }else{
        return acclum;
    }
    }else if(diskwindmode==2){
    //BH model - internal dissipation
    //return 0.1*(2./3.)*(0.01/2)*1.99*1e33*(0.3*C*0.3*C)/pow((time/2.),5./3.);
    return 0.3*1e-3*1.99*1e33*(0.3*C*0.3*C)/pow((time/1.),5./3.); //new
    }
}


/* maximum energy of electrons */
double gamma_e_max(double b_pwn)
{
    /* Eq. (21) of Murase+15, but neglecting Y */
    double eta = 1;
    
    if(diskwindmode==1){
    return 1e5*EV_TO_ERG; //BHdisk
    }else{
    return sqrt(6.0*M_PI*ELEC/eta/SIGMA_T/b_pwn); //modified
    }

}

/* possible maxium energy of photons */
double e_gamma_max(double b_pwn)
{
    /* Eq. (22) of Murase+15 in unit of [erg] */
    return gamma_e_max(b_pwn)*M_ELEC*C*C;
}

/* maximum energy of photons limited by gamma-gamma */
double e_gamma_gamma_ani(double T_ej)
{
    return pow(M_ELEC*C*C,2.0)/2.0/K_B/T_ej;
}

double e_gamma_syn_b(double b_pwn, double gamma_b)
{
    /* Eq. (24) of Murase+15 in unit of [erg] */
    
    return 3.0/2.0*H/(2.0*M_PI)*gamma_b*gamma_b*ELEC*b_pwn/M_ELEC/C;
}

double spec_non_thermal(double e_gamma, double b_pwn, double gamma_b, double T_ej)
{
  /* psr non-thermal emission injection spectrum "E*dF/dE/(eps_e*L_d) [erg^-1]" */
  /* We assume a broken power law with the low and high energy spectral indices are -p_1 and -2 */
  /* This is motivated by more detailed calculation by Murase+15 */
  
  if(diskwindmode==1){//1=super-Edd disk
    double e_gamma_min = 10.0*EV_TO_ERG;
    double e_gamma_syn_b_tmp = 1e3*EV_TO_ERG;
    double e_gamma_max_tmp = 1e5*EV_TO_ERG;
    double e_gamma_gamma_ani_tmp = e_gamma_gamma_ani(T_ej);
    double p_1 = 0.0; //modified
    double norm_fac = 0.0;
    
    if (e_gamma_gamma_ani_tmp < e_gamma_max_tmp)
      e_gamma_max_tmp = e_gamma_gamma_ani_tmp;
    
    if (e_gamma_max_tmp > e_gamma_syn_b_tmp){
      norm_fac = 1.0/(1.0/(2.0-p_1)*(1.0-pow(e_gamma_min/e_gamma_syn_b_tmp,2.0-p_1))+log(e_gamma_max_tmp/e_gamma_syn_b_tmp))/e_gamma_syn_b_tmp;
      if (e_gamma < e_gamma_syn_b_tmp && e_gamma >= e_gamma_min)
	return norm_fac*pow(e_gamma/e_gamma_syn_b_tmp,-p_1+1.0);
      else if (e_gamma > e_gamma_syn_b_tmp && e_gamma <= e_gamma_max_tmp)
	return norm_fac*pow(e_gamma/e_gamma_syn_b_tmp,-1);
      else
	return 0.0;
    } else {
      norm_fac = 1.0/(1.0/(2.0-p_1)*(1.0-pow(e_gamma_min/e_gamma_syn_b_tmp,2.0-p_1)))/e_gamma_syn_b_tmp;
      if (e_gamma < e_gamma_syn_b_tmp && e_gamma >= e_gamma_min)
	return norm_fac*pow(e_gamma/e_gamma_syn_b_tmp,-p_1+1.0);
      else
	return 0.0;
    }
  }
  else{//0=no disk, 2=disk wind
    double p_1 = 1.5; //photon index -- modified
    double e_gamma_min;
    if(diskwindmode==0){
      e_gamma_min= 1.0*EV_TO_ERG; //0=no disk for pulsar
    }else if(diskwindmode==2){
      e_gamma_min = 1.0*1e-6*EV_TO_ERG; //2=disk wind
    }
    
    double e_gamma_max_tmp = e_gamma_max(b_pwn);
    double e_gamma_gamma_ani_tmp = e_gamma_gamma_ani(T_ej);
    double e_gamma_syn_b_tmp = e_gamma_syn_b(b_pwn,gamma_b);
    double norm_fac = 0.0;
    
    if (e_gamma_gamma_ani_tmp < e_gamma_max_tmp)
      e_gamma_max_tmp = e_gamma_gamma_ani_tmp;
    
    if(e_gamma_max_tmp > e_gamma_syn_b_tmp && e_gamma_min < e_gamma_syn_b_tmp){
      norm_fac = 1.0/((1.0/(2.0-p_1))*(1.0-pow(e_gamma_min/e_gamma_syn_b_tmp,2.0-p_1))+log(e_gamma_max_tmp/e_gamma_syn_b_tmp))/e_gamma_syn_b_tmp;
      if(e_gamma < e_gamma_syn_b_tmp && e_gamma >= e_gamma_min)
        return norm_fac*pow(e_gamma/e_gamma_syn_b_tmp,-p_1+1.0);
      else if(e_gamma > e_gamma_syn_b_tmp && e_gamma <= e_gamma_max_tmp)
        return norm_fac*pow(e_gamma/e_gamma_syn_b_tmp,-1);
      else
        return 0.0;
    }
    else if(e_gamma_min > e_gamma_syn_b_tmp){
      norm_fac = 1.0/log(e_gamma_max_tmp/e_gamma_syn_b_tmp)/e_gamma_syn_b_tmp;
      if (e_gamma < e_gamma_max_tmp && e_gamma > e_gamma_min)
	return norm_fac*pow(e_gamma/e_gamma_syn_b_tmp,-1);
      else
	return 0.0;
    }
    else{
      norm_fac = 1.0/((1.0/(2.0-p_1))*(1.0-pow(e_gamma_min/e_gamma_syn_b_tmp,2.0-p_1)))/e_gamma_syn_b_tmp;
      if (e_gamma < e_gamma_syn_b_tmp && e_gamma >= e_gamma_min)
        return norm_fac*pow(e_gamma/e_gamma_syn_b_tmp,-p_1+1.0);
      else
        return 0.0;
    }
  }
}
