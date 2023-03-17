#ifndef H_ODESOLVER
#define H_ODESOLVER

#include <cmath>
#include <ctime>
#include <ostream>
#include <sstream>
#include <unordered_map> 
#include <fstream> 
#include <numeric>
#include <boost/numeric/odeint.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost::algorithm;
using namespace boost::numeric::odeint;

typedef struct TODE {
    std::vector<double> u; //Variables's vector     
    std::function<void (const std::vector<double> &u , std::vector<double> &dudt , const double /* t */)> odeSystem;
    TODE(){}
    TODE(std::function<void (const std::vector<double> &u , std::vector<double> &dudt , const double /* t */)> s){
        odeSystem = s;
    }
} ODE; 

void odesystem(const std::vector<double> &u , std::vector<double> &dudt , const double /* t */) ;

void solve(ODE *ode, double tfinal, double dt, std::vector<std::string> varNames, 
                std::vector<double> u0, std::string fname); 

void createFile(std::string name, std::vector<std::string> varNames); 

void save(std::string filename, double t, std::vector<double> values); 

void plot(std::vector<double> times, int id, std::string vname, std::string fname);

#endif 