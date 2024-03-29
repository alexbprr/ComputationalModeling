#include "OdeSolver.h"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
std::vector<std::string> colors = {"IndianRed", "Aqua", "Aquamarine", "Azure", "Beige", "Bisque", "Black", "BlanchedAlmond", "Blue"};

void odesystem(const std::vector<double> &u , std::vector<double> &dudt , const double /* t */) {
    double beta = 0.006;
    double v = 0.1; 
    double S = u[0];
    double I = u[1];  
    double R = u[2];  
    dudt[0] = -beta*S*I;
    dudt[1] = beta*S*I - v*I;
    dudt[2] = v*I;
}  

void solve(ODE *ode, double tfinal, double dt, std::vector<std::string> varNames, 
                std::vector<double> u0, std::string fname){
    double t = 0;
    ode->u = u0;
    createFile(fname,varNames);
    save(fname,0,ode->u);
    
    runge_kutta_cash_karp54<std::vector<double>> stepper;
    auto c_stepper = make_controlled(1.E-08, 1.E-08, stepper);
    
    for (t = dt; t <= tfinal; t += dt) {
        
        c_stepper.stepper().do_step(ode->odeSystem, ode->u, t, dt);            
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

void plot(std::vector<double> times, int id, std::string vname, std::string fname){    
    std::vector<std::vector<double>> data = readCSV_MultidimensionalArray(fname);
    std::vector<double> col;
    for (auto& row: data){    
        col.push_back(row.at(id));    
    }
    assert(times.size() == col.size());
    
    plt::figure_size(1200, 780);
    plt::xlabel("days");
    plt::plot(times, col, {{"label", vname}, {"color",colors[(id-1)%colors.size()]}, {"linewidth", "2"}});
    plt::legend();
    
    stringstream ss;
    ss << vname << "_plot.png";
    plt::save(ss.str());
    ss.str("");
}
