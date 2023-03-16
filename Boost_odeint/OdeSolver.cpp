#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xcsv.hpp>
#include <xtensor/xrandom.hpp>

#include "OdeSolver.h"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

void odesystem(const std::vector<double> &u , std::vector<double> &dudt , const double /* t */) {
    double rb = 0.2, 
    lnb = 0.05, 
    sn = 0.1, 
    mn = 0.1;
    double B = u[0];
    double N = u[1];    
    dudt[0] = rb*B/(1 + B) - lnb*N*B;
    dudt[1] = sn*B*N - mn*N; 
}

bool isReferenceTime(std::vector<double> times, double ct){
    int j = 0;
    for (double t: times) {
        if (abs(ct - t) <= 1.0E-06) {
            //Remove t from times 
            times.erase(times.begin()+j);
            return true; 
        }
        j++;
    }
    return false; 
}

void createFile(std::string name, std::vector<std::string> varNames){
    std::ofstream fp;
    fp.open(name);
    fp.close();
}

void save(std::string filename, double t, std::vector<double> values){
    std::ofstream fp;    
    std::stringstream ss;     
    fp.open(filename, ios::app);
    fp << t; 
    for (double v: values)
        fp << "," << v;        
    fp << endl; 
    fp.close();    
}  

void solve(ODE *ode, double tfinal, double dt, std::vector<double> times, std::vector<std::string> varNames, 
                std::vector<double> u0, std::string fname){
    double t = 0;
    ode->u = u0;
    createFile(fname,varNames);
    save(fname,0,ode->u);
    
    runge_kutta_cash_karp54<std::vector<double>> stepper;
    auto c_stepper = make_controlled(1.E-08, 1.E-08, stepper);
    
    for (t = dt; t <= tfinal; t += dt) {
        
        c_stepper.stepper().do_step(ode->odeSystem, ode->u, t, dt);
            
        if (isReferenceTime(times,t))
            save(fname,t,ode->u);    
    }
} 

std::vector<std::vector<double>> readCSV_MultidimensionalArray(std::string fname){
    std::ifstream f(fname);
    std::string line, val;                  /* string for line & value */
    std::vector<std::vector<double>> array;    /* vector of vector<int>  */

    while (std::getline (f, line)) {        /* read each line */
        std::vector<double> v;                 /* row vector v */
        std::stringstream s (line);         /* stringstream line */
        while (getline (s, val, ','))       /* get each value (',' delimited) */
            v.push_back (std::stod (val));  /* add to row vector */
        array.push_back (v);                /* add row vector to array */
    }

    return array;
}

void plot(ODE *ode, std::vector<double> times, int id, std::string vname, std::string fname){    
    std::vector<std::vector<double>> data = readCSV_MultidimensionalArray(fname);
    std::vector<double> col;
    for (auto& row: data){    
        col.push_back(row.at(id));    
    }
    assert(times.size() == col.size());
    
    plt::figure_size(1200, 780);
    plt::xlabel("days");
    plt::plot(times, col, {{"label", vname}});
    plt::legend();
    plt::save("./plots.png");
}
