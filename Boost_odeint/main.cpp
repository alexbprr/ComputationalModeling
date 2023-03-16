#include "OdeSolver.h"

int main(){
    ODE ode(odesystem); 
    std::vector<double> times;
    int tfinal = 10;
    double dt = 0.1;
    for (double t = 0; t <= tfinal; t += dt)
        times.push_back(t);
    solve(&ode,tfinal,dt,times,{"Bacteria", "Neutrophil"},{10,5},"output.csv");
    plot(&ode,times,1,"Bacteria","output.csv");
    return 0;
}