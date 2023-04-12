#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
 

#define NEQ 21
typedef double real;

real IFNUMBER_0(real t, real stim_start, real stim_end, real stim_period, real stim_duration, real stim_amplitude) {
    if(((t>=stim_start)&&(t<=stim_end))&&(((t-stim_start)-(floor(((t-stim_start)/stim_period))*stim_period))<=stim_duration)) {
        return stim_amplitude;
    } else {
        return 0.000000e+00;
    }

}

real IFNUMBER_1(real V) {
    if(V==(-4.713000e+01)) {
        return 3.200000e+00;
    } else {
        return ((3.200000e-01*(V+4.713000e+01))/(1.000000e+00-exp(((-1.000000e-01)*(V+4.713000e+01)))));
    }

}

real IFNUMBER_2(real V) {
    if(V<(-4.000000e+01)) {
        return (1.350000e-01*exp(((V+8.000000e+01)/(-6.800000e+00))));
    } else {
        return 0.000000e+00;
    }

}

real IFNUMBER_3(real V) {
    if(V<(-4.000000e+01)) {
        return ((3.560000e+00*exp((7.900000e-02*V)))+(3.100000e+05*exp((3.500000e-01*V))));
    } else {
        return (1.000000e+00/(1.300000e-01*(1.000000e+00+exp(((V+1.066000e+01)/(-1.110000e+01))))));
    }

}

real IFNUMBER_4(real V) {
    if(V<(-4.000000e+01)) {
        return (((((-1.271400e+05)*exp((2.444000e-01*V)))-(3.474000e-05*exp(((-4.391000e-02)*V))))*(V+3.778000e+01))/(1.000000e+00+exp((3.110000e-01*(V+7.923000e+01)))));
    } else {
        return 0.000000e+00;
    }

}

real IFNUMBER_5(real V) {
    if(V<(-4.000000e+01)) {
        return ((1.212000e-01*exp(((-1.052000e-02)*V)))/(1.000000e+00+exp(((-1.378000e-01)*(V+4.014000e+01)))));
    } else {
        return ((3.000000e-01*exp(((-2.535000e-07)*V)))/(1.000000e+00+exp(((-1.000000e-01)*(V+3.200000e+01)))));
    }

}

real IFNUMBER_6(real V) {
    if(fabs((V+1.410000e+01))<1.000000e-10) {
        return 1.500000e-03;
    } else {
        return ((3.000000e-04*(V+1.410000e+01))/(1.000000e+00-exp(((V+1.410000e+01)/(-5.000000e+00)))));
    }

}

real IFNUMBER_7(real V) {
    if(fabs((V-3.332800e+00))<1.000000e-10) {
        return 3.783612e-04;
    } else {
        return ((7.389800e-05*(V-3.332800e+00))/(exp(((V-3.332800e+00)/5.123700e+00))-1.000000e+00));
    }

}

real IFNUMBER_8(real V) {
    if(fabs((V-1.990000e+01))<1.000000e-10) {
        return 6.800000e-04;
    } else {
        return ((4.000000e-05*(V-1.990000e+01))/(1.000000e+00-exp(((V-1.990000e+01)/(-1.700000e+01)))));
    }

}

real IFNUMBER_9(real V) {
    if(fabs((V-1.990000e+01))<1.000000e-10) {
        return 3.150000e-04;
    } else {
        return ((3.500000e-05*(V-1.990000e+01))/(exp(((V-1.990000e+01)/9.000000e+00))-1.000000e+00));
    }

}

real IFNUMBER_10(real V) {
    if(fabs((V+1.000000e+01))<1.000000e-10) {
        return (4.579000e+00/(1.000000e+00+exp(((V+1.000000e+01)/(-6.240000e+00)))));
    } else {
        return ((1.000000e+00-exp(((V+1.000000e+01)/(-6.240000e+00))))/((3.500000e-02*(V+1.000000e+01))*(1.000000e+00+exp(((V+1.000000e+01)/(-6.240000e+00))))));
    }

}

real IFNUMBER_11(real V) {
    if(fabs((V-7.900000e+00))<1.000000e-10) {
        return ((6.000000e+00*2.000000e-01)/1.300000e+00);
    } else {
        return ((6.000000e+00*(1.000000e+00-exp(((-(V-7.900000e+00))/5.000000e+00))))/(((1.000000e+00+(3.000000e-01*exp(((-(V-7.900000e+00))/5.000000e+00))))*1.000000e+00)*(V-7.900000e+00)));
    }

}

void set_initial_conditions(real *x0, real *values) { 

    x0[0] = values[0]; //V
    x0[1] = values[1]; //m
    x0[2] = values[2]; //h
    x0[3] = values[3]; //j
    x0[4] = values[4]; //oa
    x0[5] = values[5]; //oi
    x0[6] = values[6]; //ua
    x0[7] = values[7]; //ui
    x0[8] = values[8]; //xr
    x0[9] = values[9]; //xs
    x0[10] = values[10]; //d
    x0[11] = values[11]; //f
    x0[12] = values[12]; //f_Ca
    x0[13] = values[13]; //u
    x0[14] = values[14]; //v
    x0[15] = values[15]; //w
    x0[16] = values[16]; //Na_i
    x0[17] = values[17]; //K_i
    x0[18] = values[18]; //Ca_i
    x0[19] = values[19]; //Ca_up
    x0[20] = values[20]; //Ca_rel

}

static int solve_model(real time, real *sv, real *rDY) {

    //State variables
    const real V =  sv[0];
    const real m =  sv[1];
    const real h =  sv[2];
    const real j =  sv[3];
    const real oa =  sv[4];
    const real oi =  sv[5];
    const real ua =  sv[6];
    const real ui =  sv[7];
    const real xr =  sv[8];
    const real xs =  sv[9];
    const real d =  sv[10];
    const real f =  sv[11];
    const real f_Ca =  sv[12];
    const real u =  sv[13];
    const real v =  sv[14];
    const real w =  sv[15];
    const real Na_i =  sv[16];
    const real K_i =  sv[17];
    const real Ca_i =  sv[18];
    const real Ca_up =  sv[19];
    const real Ca_rel =  sv[20];

    //Parameters
    real stim_amplitude = (-2.000000e+03);
    real stim_start = 5.000000e+00;
    real stim_end = 1.000000e+05;
    real stim_period = 1.000000e+03;
    real stim_duration = 2.000000e+00;
    real Cm = 1.000000e+02;
    real g_Na = 7.800000e+00;
    real R = 8.314300e+00;
    real T = 3.100000e+02;
    real F = 9.648670e+01;
    real Na_o = 1.400000e+02;
    real K_o = 5.400000e+00;
    real g_K1 = 9.000000e-02;
    real g_to = 1.652000e-01;
    real K_Q10 = 3.000000e+00;
    real g_Kr = 2.941176e-02;
    real g_Ks = 1.294118e-01;
    real g_Ca_L = 1.237500e-01;
    real i_NaK_max = 5.993387e-01;
    real Km_Na_i = 1.000000e+01;
    real Km_K_o = 1.500000e+00;
    real Ca_o = 1.800000e+00;
    real g_B_Na = 6.744375e-04;
    real g_B_Ca = 1.131000e-03;
    real g_B_K = 0.000000e+00;
    real I_NaCa_max = 1.600000e+03;
    real gamma = 3.500000e-01;
    real K_mNa = 8.750000e+01;
    real K_mCa = 1.380000e+00;
    real K_sat = 1.000000e-01;
    real i_CaP_max = 2.750000e-01;
    real K_rel = 3.000000e+01;
    real tau_tr = 1.800000e+02;
    real I_up_max = 5.000000e-03;
    real K_up = 9.200000e-04;
    real Ca_up_max = 1.500000e+01;
    real CMDN_max = 5.000000e-02;
    real Km_CMDN = 2.380000e-03;
    real TRPN_max = 7.000000e-02;
    real Km_TRPN = 5.000000e-04;
    real CSQN_max = 1.000000e+01;
    real Km_CSQN = 8.000000e-01;
    real V_cell = 2.010000e+04;
    real i_st = 0.000000e+00;
    i_st = IFNUMBER_0(time, stim_start, stim_end, stim_period, stim_duration, stim_amplitude);
    real E_Na = (((R*T)/F)*log((Na_o/Na_i)));
    real alpha_m = 0.000000e+00;
    alpha_m = IFNUMBER_1(V);
    real beta_m = (8.000000e-02*exp(((-V)/1.100000e+01)));
    real alpha_h = 0.000000e+00;
    alpha_h = IFNUMBER_2(V);
    real beta_h = 0.000000e+00;
    beta_h = IFNUMBER_3(V);
    real alpha_j = 0.000000e+00;
    alpha_j = IFNUMBER_4(V);
    real beta_j = 0.000000e+00;
    beta_j = IFNUMBER_5(V);
    real E_K = (((R*T)/F)*log((K_o/K_i)));
    real alpha_oa = (6.500000e-01*pow((exp(((V-(-1.000000e+01))/(-8.500000e+00)))+exp((((V-(-1.000000e+01))-4.000000e+01)/(-5.900000e+01)))), (-1.000000e+00)));
    real beta_oa = (6.500000e-01*pow((2.500000e+00+exp((((V-(-1.000000e+01))+7.200000e+01)/1.700000e+01))), (-1.000000e+00)));
    real oa_infinity = pow((1.000000e+00+exp((((V-(-1.000000e+01))+1.047000e+01)/(-1.754000e+01)))), (-1.000000e+00));
    real alpha_oi = pow((1.853000e+01+(1.000000e+00*exp((((V-(-1.000000e+01))+1.037000e+02)/1.095000e+01)))), (-1.000000e+00));
    real beta_oi = pow((3.556000e+01+(1.000000e+00*exp((((V-(-1.000000e+01))-8.740000e+00)/(-7.440000e+00))))), (-1.000000e+00));
    real oi_infinity = pow((1.000000e+00+exp((((V-(-1.000000e+01))+3.310000e+01)/5.300000e+00))), (-1.000000e+00));
    real g_Kur = (5.000000e-03+(5.000000e-02/(1.000000e+00+exp(((V-1.500000e+01)/(-1.300000e+01))))));
    real alpha_ua = (6.500000e-01*pow((exp(((V-(-1.000000e+01))/(-8.500000e+00)))+exp((((V-(-1.000000e+01))-4.000000e+01)/(-5.900000e+01)))), (-1.000000e+00)));
    real beta_ua = (6.500000e-01*pow((2.500000e+00+exp((((V-(-1.000000e+01))+7.200000e+01)/1.700000e+01))), (-1.000000e+00)));
    real ua_infinity = pow((1.000000e+00+exp((((V-(-1.000000e+01))+2.030000e+01)/(-9.600000e+00)))), (-1.000000e+00));
    real alpha_ui = pow((2.100000e+01+(1.000000e+00*exp((((V-(-1.000000e+01))-1.950000e+02)/(-2.800000e+01))))), (-1.000000e+00));
    real beta_ui = (1.000000e+00/exp((((V-(-1.000000e+01))-1.680000e+02)/(-1.600000e+01))));
    real ui_infinity = pow((1.000000e+00+exp((((V-(-1.000000e+01))-1.094500e+02)/2.748000e+01))), (-1.000000e+00));
    real alpha_xr = 0.000000e+00;
    alpha_xr = IFNUMBER_6(V);
    real beta_xr = 0.000000e+00;
    beta_xr = IFNUMBER_7(V);
    real xr_infinity = pow((1.000000e+00+exp(((V+1.410000e+01)/(-6.500000e+00)))), (-1.000000e+00));
    real alpha_xs = 0.000000e+00;
    alpha_xs = IFNUMBER_8(V);
    real beta_xs = 0.000000e+00;
    beta_xs = IFNUMBER_9(V);
    real xs_infinity = pow((1.000000e+00+exp(((V-1.990000e+01)/(-1.270000e+01)))), (-5.000000e-01));
    real i_Ca_L = (((((Cm*g_Ca_L)*d)*f)*f_Ca)*(V-6.500000e+01));
    real d_infinity = pow((1.000000e+00+exp(((V+1.000000e+01)/(-8.000000e+00)))), (-1.000000e+00));
    real tau_d = 0.000000e+00;
    tau_d = IFNUMBER_10(V);
    real f_infinity = (exp(((-(V+2.800000e+01))/6.900000e+00))/(1.000000e+00+exp(((-(V+2.800000e+01))/6.900000e+00))));
    real tau_f = (9.000000e+00*pow(((1.970000e-02*exp(((-pow(3.370000e-02, 2.000000e+00))*pow((V+1.000000e+01), 2.000000e+00))))+2.000000e-02), (-1.000000e+00)));
    real f_Ca_infinity = pow((1.000000e+00+(Ca_i/3.500000e-04)), (-1.000000e+00));
    real tau_f_Ca = 2.000000e+00;
    real sigma = ((1.000000e+00/7.000000e+00)*(exp((Na_o/6.730000e+01))-1.000000e+00));
    real E_Ca = (((R*T)/(2.000000e+00*F))*log((Ca_o/Ca_i)));
    real i_NaCa = (((Cm*I_NaCa_max)*(((exp((((gamma*F)*V)/(R*T)))*pow(Na_i, 3.000000e+00))*Ca_o)-((exp(((((gamma-1.000000e+00)*F)*V)/(R*T)))*pow(Na_o, 3.000000e+00))*Ca_i)))/(((pow(K_mNa, 3.000000e+00)+pow(Na_o, 3.000000e+00))*(K_mCa+Ca_o))*(1.000000e+00+(K_sat*exp(((((gamma-1.000000e+00)*V)*F)/(R*T)))))));
    real i_CaP = (((Cm*i_CaP_max)*Ca_i)/(5.000000e-04+Ca_i));
    real i_rel = ((((K_rel*pow(u, 2.000000e+00))*v)*w)*(Ca_rel-Ca_i));
    real tau_u = 8.000000e+00;
    real tau_w = 0.000000e+00;
    tau_w = IFNUMBER_11(V);
    real w_infinity = (1.000000e+00-pow((1.000000e+00+exp(((-(V-4.000000e+01))/1.700000e+01))), (-1.000000e+00)));
    real i_tr = ((Ca_up-Ca_rel)/tau_tr);
    real i_up = (I_up_max/(1.000000e+00+(K_up/Ca_i)));
    real i_up_leak = ((I_up_max*Ca_up)/Ca_up_max);
    real Ca_CMDN = ((CMDN_max*Ca_i)/(Ca_i+Km_CMDN));
    real Ca_TRPN = ((TRPN_max*Ca_i)/(Ca_i+Km_TRPN));
    real Ca_CSQN = ((CSQN_max*Ca_rel)/(Ca_rel+Km_CSQN));
    real V_i = (V_cell*6.800000e-01);
    real V_rel = (4.800000e-03*V_cell);
    real V_up = (5.520000e-02*V_cell);
    real B2 = ((1.000000e+00+((TRPN_max*Km_TRPN)/pow((Ca_i+Km_TRPN), 2.000000e+00)))+((CMDN_max*Km_CMDN)/pow((Ca_i+Km_CMDN), 2.000000e+00)));
    real m_inf = (alpha_m/(alpha_m+beta_m));
    real tau_m = (1.000000e+00/(alpha_m+beta_m));
    real h_inf = (alpha_h/(alpha_h+beta_h));
    real tau_h = (1.000000e+00/(alpha_h+beta_h));
    real j_inf = (alpha_j/(alpha_j+beta_j));
    real tau_j = (1.000000e+00/(alpha_j+beta_j));
    real i_K1 = (((Cm*g_K1)*(V-E_K))/(1.000000e+00+exp((7.000000e-02*(V+8.000000e+01)))));
    real i_to = ((((Cm*g_to)*pow(oa, 3.000000e+00))*oi)*(V-E_K));
    real tau_oa = (pow((alpha_oa+beta_oa), (-1.000000e+00))/K_Q10);
    real tau_oi = (pow((alpha_oi+beta_oi), (-1.000000e+00))/K_Q10);
    real i_Kur = ((((Cm*g_Kur)*pow(ua, 3.000000e+00))*ui)*(V-E_K));
    real tau_ua = (pow((alpha_ua+beta_ua), (-1.000000e+00))/K_Q10);
    real tau_ui = (pow((alpha_ui+beta_ui), (-1.000000e+00))/K_Q10);
    real i_Kr = ((((Cm*g_Kr)*xr)*(V-E_K))/(1.000000e+00+exp(((V+1.500000e+01)/2.240000e+01))));
    real tau_xr = pow((alpha_xr+beta_xr), (-1.000000e+00));
    real i_Ks = (((Cm*g_Ks)*pow(xs, 2.000000e+00))*(V-E_K));
    real tau_xs = (5.000000e-01*pow((alpha_xs+beta_xs), (-1.000000e+00)));
    real f_NaK = pow(((1.000000e+00+(1.245000e-01*exp(((((-1.000000e-01)*F)*V)/(R*T)))))+((3.650000e-02*sigma)*exp((((-F)*V)/(R*T))))), (-1.000000e+00));
    real i_B_Na = ((Cm*g_B_Na)*(V-E_Na));
    real i_B_Ca = ((Cm*g_B_Ca)*(V-E_Ca));
    real i_B_K = ((Cm*g_B_K)*(V-E_K));
    real i_Na = (((((Cm*g_Na)*pow(m, 3.000000e+00))*h)*j)*(V-E_Na));
    real i_NaK = ((((((Cm*i_NaK_max)*f_NaK)*1.000000e+00)/(1.000000e+00+pow((Km_Na_i/Na_i), 1.500000e+00)))*K_o)/(K_o+Km_K_o));
    real Fn = (1.000000e+03*(((1.000000e-15*V_rel)*i_rel)-((1.000000e-15/(2.000000e+00*F))*((5.000000e-01*i_Ca_L)-(2.000000e-01*i_NaCa)))));
    real B1 = ((((2.000000e+00*i_NaCa)-((i_CaP+i_Ca_L)+i_B_Ca))/((2.000000e+00*V_i)*F))+(((V_up*(i_up_leak-i_up))+(i_rel*V_rel))/V_i));
    real u_infinity = pow((1.000000e+00+exp(((-(Fn-3.417500e-13))/1.367000e-15))), (-1.000000e+00));
    real tau_v = (1.910000e+00+(2.090000e+00*pow((1.000000e+00+exp(((-(Fn-3.417500e-13))/1.367000e-15))), (-1.000000e+00))));
    real v_infinity = (1.000000e+00-pow((1.000000e+00+exp(((-(Fn-6.835000e-14))/1.367000e-15))), (-1.000000e+00)));
    rDY[0] = ((-((((((((((((i_Na+i_K1)+i_to)+i_Kur)+i_Kr)+i_Ks)+i_B_Na)+i_B_Ca)+i_NaK)+i_CaP)+i_NaCa)+i_Ca_L)+i_st))/Cm);
    rDY[1] = ((m_inf-m)/tau_m);
    rDY[2] = ((h_inf-h)/tau_h);
    rDY[3] = ((j_inf-j)/tau_j);
    rDY[4] = ((oa_infinity-oa)/tau_oa);
    rDY[5] = ((oi_infinity-oi)/tau_oi);
    rDY[6] = ((ua_infinity-ua)/tau_ua);
    rDY[7] = ((ui_infinity-ui)/tau_ui);
    rDY[8] = ((xr_infinity-xr)/tau_xr);
    rDY[9] = ((xs_infinity-xs)/tau_xs);
    rDY[10] = ((d_infinity-d)/tau_d);
    rDY[11] = ((f_infinity-f)/tau_f);
    rDY[12] = ((f_Ca_infinity-f_Ca)/tau_f_Ca);
    rDY[13] = ((u_infinity-u)/tau_u);
    rDY[14] = ((v_infinity-v)/tau_v);
    rDY[15] = ((w_infinity-w)/tau_w);
    rDY[16] = ((((-3.000000e+00)*i_NaK)-(((3.000000e+00*i_NaCa)+i_B_Na)+i_Na))/(V_i*F));
    rDY[17] = (((2.000000e+00*i_NaK)-(((((i_K1+i_to)+i_Kur)+i_Kr)+i_Ks)+i_B_K))/(V_i*F));
    rDY[18] = (B1/B2);
    rDY[19] = (i_up-(i_up_leak+((i_tr*V_rel)/V_up)));
    rDY[20] = ((i_tr-i_rel)*pow((1.000000e+00+((CSQN_max*Km_CSQN)/pow((Ca_rel+Km_CSQN), 2.000000e+00))), (-1.000000e+00)));

	return 0;  

}

void solve_ode(real *sv, float final_time, char *file_name) {

    real rDY[NEQ];

    real reltol = 1e-5;
    real abstol = 1e-5;
    real _tolerances_[NEQ];
    real _aux_tol = 0.0;
    //initializes the variables
    real dt = 0.000001;
    real time_new = 0.0;
    real previous_dt = dt;

    real edos_old_aux_[NEQ];
    real edos_new_euler_[NEQ];
    real *_k1__ = (real*) malloc(sizeof(real)*NEQ);
    real *_k2__ = (real*) malloc(sizeof(real)*NEQ);
    real *_k_aux__;

    const real _beta_safety_ = 0.8;

    const real __tiny_ = pow(abstol, 2.0f);

    if(time_new + dt > final_time) {
       dt = final_time - time_new;
    }

    solve_model(time_new, sv, rDY);
    time_new += dt;

    for(int i = 0; i < NEQ; i++){
        _k1__[i] = rDY[i];
    }

    FILE *f = fopen(file_name, "w");
    fprintf(f, "#t, V, m, h, j, oa, oi, ua, ui, xr, xs, d, f, f_Ca, u, v, w, Na_i, K_i, Ca_i, Ca_up, Ca_rel\n");

    real min[NEQ];
    real max[NEQ];

    for(int i = 0; i < NEQ; i++) {
       min[i] = DBL_MAX;
       max[i] = DBL_MIN;
    }

    while(1) {

        for(int i = 0; i < NEQ; i++) {
            //stores the old variables in a vector
            edos_old_aux_[i] = sv[i];
            //computes euler method
            edos_new_euler_[i] = _k1__[i] * dt + edos_old_aux_[i];
            //steps ahead to compute the rk2 method
            sv[i] = edos_new_euler_[i];
        }

        time_new += dt;
        solve_model(time_new, sv, rDY);
        time_new -= dt;//step back

        double greatestError = 0.0, auxError = 0.0;
        for(int i = 0; i < NEQ; i++) {
            // stores the new evaluation
            _k2__[i] = rDY[i];
            _aux_tol = fabs(edos_new_euler_[i]) * reltol;
            _tolerances_[i] = (abstol > _aux_tol) ? abstol : _aux_tol;

            // finds the greatest error between  the steps
            auxError = fabs(((dt / 2.0) * (_k1__[i] - _k2__[i])) / _tolerances_[i]);

            greatestError = (auxError > greatestError) ? auxError : greatestError;
        }
        ///adapt the time step
        greatestError += __tiny_;
        previous_dt = dt;
        ///adapt the time step
        dt = _beta_safety_ * dt * sqrt(1.0f/greatestError);

        if (time_new + dt > final_time) {
            dt = final_time - time_new;
        }

        //it doesn't accept the solution
        if ((greatestError >= 1.0f) && dt > 0.00000001) {
            //restore the old values to do it again
            for(int i = 0;  i < NEQ; i++) {
                sv[i] = edos_old_aux_[i];
            }
            //throw the results away and compute again
        } else{//it accepts the solutions

            if (time_new + dt > final_time) {
                dt = final_time - time_new;
            }

            _k_aux__ = _k2__;
            _k2__	= _k1__;
            _k1__	= _k_aux__;

            //it steps the method ahead, with euler solution
            for(int i = 0; i < NEQ; i++){
                sv[i] = edos_new_euler_[i];
            }
            fprintf(f, "%lf ", time_new);
            for(int i = 0; i < NEQ; i++) {
                fprintf(f, "%lf ", sv[i]);
                if(sv[i] < min[i]) min[i] = sv[i];
                if(sv[i] > max[i]) max[i] = sv[i];
            }

            fprintf(f, "\n");

            if(time_new + previous_dt >= final_time) {
                if(final_time == time_new) {
                    break;
                } else if(time_new < final_time) {
                    dt = previous_dt = final_time - time_new;
                    time_new += previous_dt;
                    break;
                }
            } else {
                time_new += previous_dt;
            }

        }
    }

    char *min_max = malloc(strlen(file_name) + 9);
    sprintf(min_max, "%s_min_max", file_name);
    FILE* min_max_file = fopen(min_max, "w");
    for(int i = 0; i < NEQ; i++) {
        fprintf(min_max_file, "%lf;%lf\n", min[i], max[i]);
    }
    fclose(min_max_file);
    free(min_max);
    
    free(_k1__);
    free(_k2__);
}
int main(int argc, char **argv) {

	real *x0 = (real*) malloc(sizeof(real)*NEQ);

    real values[21];
    values[0] = (-8.118000e+01); //V
    values[1] = 2.908000e-03; //m
    values[2] = 9.649000e-01; //h
    values[3] = 9.775000e-01; //j
    values[4] = 3.043000e-02; //oa
    values[5] = 9.992000e-01; //oi
    values[6] = 4.966000e-03; //ua
    values[7] = 9.986000e-01; //ui
    values[8] = 3.296000e-05; //xr
    values[9] = 1.869000e-02; //xs
    values[10] = 1.367000e-04; //d
    values[11] = 9.996000e-01; //f
    values[12] = 7.755000e-01; //f_Ca
    values[13] = 2.350000e-112; //u
    values[14] = 1.000000e+00; //v
    values[15] = 9.992000e-01; //w
    values[16] = 1.117000e+01; //Na_i
    values[17] = 1.390000e+02; //K_i
    values[18] = 1.013000e-04; //Ca_i
    values[19] = 1.488000e+00; //Ca_up
    values[20] = 1.488000e+00; //Ca_rel
	set_initial_conditions(x0, values);

	solve_ode(x0, strtod(argv[1], NULL), argv[2]);


	return (0);
}