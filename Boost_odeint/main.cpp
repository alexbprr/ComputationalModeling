#include "OdeSolver.h"

int main(){
    std::vector<double> times;
    int tfinal = 10;
    double dt = 0.1;

    ODE ode(odesystem);     
    
    for (double t = 0; t <= tfinal; t += dt)
        times.push_back(t);
    
    std::vector<std::string> populations = {"Bacteria", "Neutrophil"};
    std::string fname = "output.csv";
    
    solve(&ode,tfinal,dt,times,populations,{10,5},fname);

    for (int i = 0; i < populations.size(); i++)
        plot(times,i+1,populations.at(i),fname);

    return 0;
}