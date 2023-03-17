#include "OdeSolver.h"

int main(){
    int tfinal = 10;
    double dt = 0.1;

    ODE ode(odesystem);         
    
    std::vector<std::string> populations = {"Bacteria", "Neutrophil"};
    std::string fname = "output.csv";
    
    solve(&ode,tfinal,dt,populations,{10,5},fname);

    std::vector<double> times;
    for (double t = 0; t <= tfinal; t += dt)
        times.push_back(t);
        
    for (int i = 0; i < populations.size(); i++)
        plot(times,i+1,populations.at(i),fname);

    return 0;
}