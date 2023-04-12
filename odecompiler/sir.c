#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
 

#define NEQ 3
typedef double real;
real n = 1.000000e+03;
real init_i = 3.000000e+00;

void set_initial_conditions(real *x0, real *values) { 

    x0[0] = values[0]; //S
    x0[1] = values[1]; //I
    x0[2] = values[2]; //R

}

static int solve_model(real time, real *sv, real *rDY) {

    //State variables
    const real S =  sv[0];
    const real I =  sv[1];
    const real R =  sv[2];

    //Parameters
    real beta = (1.000000e+00/n);
    real gamma = 4.000000e-02;
    rDY[0] = (((-beta)*S)*I);
    rDY[1] = (((beta*S)*I)-(gamma*I));
    rDY[2] = (gamma*I);

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
    fprintf(f, "#t, S, I, R\n");

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

    real values[3];
    values[0] = (n-init_i); //S
    values[1] = init_i; //I
    values[2] = 0.000000e+00; //R
	set_initial_conditions(x0, values);

	solve_ode(x0, strtod(argv[1], NULL), argv[2]);


	return (0);
}