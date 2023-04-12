#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <sstream>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace std;
std::vector<std::string> colors = {"IndianRed", "Aqua", "Aquamarine", "Azure", "Beige", "Bisque", "Black", "BlanchedAlmond", "Blue", "Green", "Yellow", "Purple", "Grey"};

void createFile(std::string name){
    std::ofstream fp;
    fp.open(name);
    fp.close();
}

template <typename T>
void save(std::string filename, double t, std::vector<T> values){
    std::ofstream fp;    
    std::stringstream ss;     
    fp.open(filename, std::ios::app);
    fp << t; 
    for (T v: values)
        fp << "," << v;        
    fp << endl; 
    fp.close();    
}

template <typename T>
std::vector<std::vector<T>> readCSV_MultidimensionalArray(std::string fname){
    std::ifstream f(fname);
    std::string line, val;                  
    std::vector<std::vector<T>> array;    

    while (std::getline (f, line)) {      
        std::vector<T> v;                 
        std::stringstream s (line);        
        while (getline (s, val, ','))      
            v.push_back (std::stod (val)); 
        array.push_back (v);               
    }

    return array;
}

void plot(std::vector<double> times, int id, std::string vname, std::string fname){    
    std::vector<std::vector<double>> data = readCSV_MultidimensionalArray<double>(fname);
    std::vector<double> col;
    for (auto& row: data){    
        col.push_back(row.at(id));    
    }
    
    plt::figure_size(1200, 780);
    plt::xlabel("days");
    plt::plot(times, col, {{"label", vname}, {"color",colors[(id-1)%colors.size()]}, {"linewidth", "2"}});
    plt::legend();
    
    std::stringstream ss;
    ss << vname << "_plot.png";
    plt::save(ss.str());
    ss.str("");
}

// Function to calculate propensity values for each reaction
std::vector<double> calcPropensity(double A, double B, double C, double k1, double k2) {
    std::vector<double> a(2);
    a[0] = k1 * A * B;
    a[1] = k2 * C;
    return a;
}

// Function to update the system state based on a chosen reaction
void updateState(int& A, int& B, int& C, int reaction) {
    switch(reaction) {
        case 0: // Reaction 1
            C++;
            A--;
            B--;
            break;
        case 1: // Reaction 2
            A++;
            B++;
            C--;
            break;
    }
}

// Function to choose the next reaction to occur based on the current propensity values
int chooseReaction(double a0, const std::vector<double>& a, std::mt19937& gen) {
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double r1 = dis(gen);
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        sum += a[i];
        if (r1 < sum/a0) {
            return i;
        }
    }
    return -1; // No reaction chosen
}

void gillespie(std::string fname){
    // Initialize system state and parameters
    int A = 100;
    int B = 50;
    int C = 0;
    double k1 = 0.01;
    double k2 = 0.03;
    //int num_steps = 50;

    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());    
    std::vector<double> times;

    createFile(fname);
    
    // Main loop
    for (double t = 0; t < 50;) {
        std::uniform_real_distribution<> dis(0.0, 1.0);
        // Calculate propensity values for each reaction
        std::vector<double> a = calcPropensity(A, B, C, k1, k2);
        
        double a0 = std::accumulate(a.begin(), a.end(), 0.0);
        
        // Choose the next reaction to occur based on the current propensity values
        int reaction = chooseReaction(a0, a, gen);
        
        // If no reaction chosen, break out of loop
        if (reaction == -1) {
            std::cout << "No reaction chosen!" << std::endl;
            break;
        }

        // Update the system state based on the chosen reaction
        updateState(A, B, C, reaction);
        
        // Calculate the time elapsed for this step
        double r2 = dis(gen);
        double tau = (1.0/a0) * std::log(1.0/r2);

        // Update the time
        t += tau;
        
        std::vector<int> values = {A, B, C};
        times.push_back(t);
        save(fname, t, values);
    }

    plot(times,1,"A",fname);
    plot(times,2,"B",fname);
    plot(times,3,"C",fname);
}

// Main function to run the Gillespie algorithm
int main() {
    gillespie("gillespie_output.csv");
    return 0;
}