#include <cmath>
#include <ctime>
#include <map> 
#include <cstdlib> 
#include <fstream> 
#include <numeric>
#include <iostream> 
#include <boost/numeric/odeint.hpp>
#include "de/DifferentialEvolution.h"
//#include "matplotlibcpp.h"
//namespace plt = matplotlibcpp;

using namespace std;
using namespace de;
using namespace boost::numeric::odeint;

runge_kutta_cash_karp54<std::vector<double>> stepper; 
double N0;
double r, K, a, m;

void odesystem( const std::vector<double> &u , std::vector<double> &dudt , const double /* t */ ) {
    double H = u[0];
    double P = u[1];        
    dudt[0] = r*H - a*H*P;
    dudt[1] = a*H*P - m*P;
}

vector<double> advance(double t, double dt, std::vector<double> u){ //to do: pass the parameters 
    stepper.do_step(odesystem, u, t, dt);    
    return u;
}

std::vector<std::vector<double>> readCSV_to_MultidimensionalArray(std::string fname){
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

class TestDE : public de::IOptimizable{
private:    
    double tfinal; 
    double dt;
    std::vector<std::vector<double> > data;

public:    
    TestDE(double tf, double deltat): tfinal(tf), dt(deltat){}
    void init(std::string fname){
        data = readCSV_to_MultidimensionalArray(fname);        
    }
    double EvaluateCost(std::vector<double> inputs) const override {        
        int s = data[0].size()-1;
        std::vector<double> u;
        u.reserve(s);
        u.resize(s);

        u[0] = 10; //condicao inicial da presa
        u[1] = 2; //condicao inical do predador
        r = inputs[0]; 
        a = inputs[1];
        m = inputs[2];

        double errorH = 0, errorP = 0, sumexactH = 0, sumexactP = 0; 
        int i = 0;
        for (double t = 0; t <= tfinal; t += dt){

            if (abs(t - data[i][0]) < 0.01){
                double H = data[i][1];
                double P = data[i][2];
                errorH += (u[0] - H)*(u[0] - H); //Soma dos quadrados dos erros entre N obtido pela EDO e o dado de N
                sumexactH += H*H; //"Valor exato da população"
                errorP += (u[1] - P)*(u[1] - P); //Soma dos quadrados dos erros entre N obtido pela EDO e o dado de N
                sumexactP += P*P;

                i++; 
            }
            //std::cout << data.size() << std::endl;
            if (i >= data.size()) 
                break;

            u = advance(t,dt,u);
        }
        
        errorH = sqrt(errorH/sumexactH); //Erro norma 2 */  
        errorP = sqrt(errorP/sumexactP);
        return errorH + errorP;
    }

    unsigned int NumberOfParameters() const override{
        return 3;
    }

    std::vector<Constraints> GetConstraints() const override{
        std::vector<Constraints> constr(NumberOfParameters());
        constr[0] = Constraints(0, 1, true);
        constr[1] = Constraints(0, 1, true);
        constr[2] = Constraints(0, 1, true);
        return constr;
    }

    static bool terminationCondition(const DifferentialEvolution& de){
        if (de.GetBestCost() <= 0.01)
            return true;
        return false; 
    }
};

int main(){
    TestDE deInstance(50,0.01);
    deInstance.init("data/data2.csv");
    int populationSize = 100, maxIterations = 30; 
    de::DifferentialEvolution de(deInstance, populationSize, 
                std::time(nullptr), true, TestDE::terminationCondition);
    std::pair<std::vector<double>,std::vector<double>> costs = 
          de.Optimize(maxIterations, true);

    //{{"label", vname}, {"color",colors[(id-1)%colors.size()]}, {"linewidth", "2"}}
    /*plt::figure_size(1200, 780);
    plt::plot(costs.first, costs.second, {{"label", "Cost"}});
    plt::xlabel("iteration");
    plt::title("DE cost evolution");
    plt::legend();
    plt::save("./de_cost.png");*/
    return 0;
}
