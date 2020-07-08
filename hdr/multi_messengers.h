struct _input load_inputfile();

struct _input{
  double m_ns;
  double r_ns;
  double b_t;
  double b_p;
  double ome_0;
  double ene_sn;
  double m_ej;
  double m_Ni;
  double m_rpe;
  double r_0;
  double delta;
  double albd_fac;
  double kappa;
  double d_l;
};


double i_ns(double m_ns, double r_ns);
double ene_rot(double m_ns, double r_ns, double ome);
double ene_b(double m_ns, double r_ns, double b_t);
double eps_b(double m_ns, double r_ns, double b_t);
double l_gw(double m_ns, double r_ns, double eps, double ome);
double l_d(double r_ns, double b_p, double ome, double time);
double l_disk(double time);
double vol_ej(double r_ej);
double tau_ej(double m_ej, double delta, double r_ej, double kappa);
double rho_env(double r_ej);

void lc_sn(double m_ns, double r_ns, double b_t, double b_p, double ome_0, double ene_sn, double m_ej, double m_Ni, double m_rpe, double r_0, double delta, double d_l, double albd_fac, double kappa,
           char s_P0[], char s_Bp[], char s_Bt[], char s_mej[]);

double s_f_num(double f);

void ab_AB_mag(double L_bol, double T_ph, double *U_ab, double *B_ab, double *V_ab, double *R_ab, double *I_ab);
void app_AB_mag(double L_bol, double T_ph, double d_L, double *U_app, double *B_app, double *V_app, double *R_app, double *I_app);
double ab_bol_mag(double L_bol);

double l_Ni(double m_Ni, double t);
double l_Co(double m_Ni, double t);
double l_rpe(double m_rpe, double t);
double gamma_ene_depo_frac_Compton(double e_gamma);
double sigma_kn(double e_gamma);
double sigma_BH_p(double e_gamma);
double gamma_inelas_Compton(double e_gamma);
double gamma_inelas_BH(double e_gamma);
double kappa_bf(double e_gamma, double Z_eff);
double f_gamma_dep(double e_gamma, double m_ej, double delta, double r_ej, double albd_fac);
double f_gamma_esc(double e_gamma, double m_ej, double delta, double r_ej, double albd_fac);
double fac_Ni_dep(double m_ej, double delta, double r_ej, double albd_fac);
double fac_Co_dep(double m_ej, double delta, double r_ej, double albd_fac);
double fac_psr_dep(double m_ej, double delta, double r_ej, double albd_fac, double b_pwn, double gamma_b, double T_ej);
double e_gamma_max(double b_pwn);
double e_gamma_gamma_ani(double T_ej);
double e_gamma_syn_b(double b_pwn, double gamma_b);
double spec_non_thermal(double e_gamma, double b_pwn, double gamma_b, double T_ej);
