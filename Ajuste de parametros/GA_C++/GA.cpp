#include <stdlib.h>
#include <cstdio>
#include <memory.h> 
#include <time.h>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <boost/numeric/odeint.hpp>

template <typename T>
std::pair<T*, T*> selectIndividuals(T** population, size_t populationSize, float* fitnesses, float upperBound) {
    // Get a number between 0 and upper bound
    std::pair<T*, T*> chosen = std::make_pair(population[rand() % populationSize], population[rand() % populationSize]);
    float probA = static_cast<float>(rand()) / RAND_MAX * upperBound;
    float probB = static_cast<float>(rand()) / RAND_MAX * upperBound;
    size_t i = 0;
    // Rotate the roulette
    while ((probA > 0 || probB > 0) && i < populationSize) {
        if (probA > 0)
        {
            probA -= fitnesses[i];
            if (probA <= 0) chosen.first = population[i];
        }
        if (probB > 0)
        {
            probB -= fitnesses[i];
            if (probB <= 0) chosen.second = population[i];
        }
        i++;
    }
    return chosen;
}

template <typename T>
std::pair<T*, T*> crossover(std::pair<T*, T*> parents, size_t chromosomeSize, float crossoverRate, float mutationChance, void (*mutate)(T*, size_t))
{
    size_t sizeA = static_cast<size_t>(chromosomeSize * crossoverRate);
    size_t sizeB = chromosomeSize - sizeA;
    T* childA = new T[chromosomeSize];
    T* childB = new T[chromosomeSize];
    // Copy first part of chromosomes
    memcpy(childA, parents.first, sizeA * sizeof(T));
    memcpy(childB, parents.second, sizeA * sizeof(T));
    // Copy second part of chromosomes
    memcpy(childA + sizeA, parents.second + sizeA, sizeB * sizeof(T));
    memcpy(childB + sizeA, parents.first + sizeA, sizeB * sizeof(T));
    // Apply mutation
    for (size_t i = 0; i < chromosomeSize; i++)
    {
        if (static_cast<float>(rand()) / RAND_MAX <= mutationChance) mutate(childA, i);
        if (static_cast<float>(rand()) / RAND_MAX <= mutationChance) mutate(childB, i);
    }
    return std::make_pair(childA, childB);
}

template <typename T>
float* getFitnesses(T** population, size_t populationSize, float (*getFitness)(T*), bool maximization)
{
    float* fitnesses = new float[populationSize];
    for (size_t i = 0; i < populationSize; i++){
        if (getFitness != nullptr)
            fitnesses[i] = getFitness(population[i]);
        else 
            fitnesses[i] = 0; 
    }
    // Invert values on minimization problems
    if (maximization == false) 
    {
        float maxFitness = std::numeric_limits<float>::min();
        for (size_t i = 0; i < populationSize; i++)
            if (fitnesses[i] > maxFitness)
                maxFitness = fitnesses[i];

        for (size_t i = 0; i < populationSize; i++)        
            fitnesses[i] = maxFitness - fitnesses[i];
    }
    return fitnesses;
}

template <typename T>
void clear(T** population, size_t populationSize)
{
    for (size_t i = 0; i < populationSize; i++)
    {
        delete population[i];
    }
    delete[] population;
}

template <typename T>
T** geneticAlgorithm(size_t chromosomeSize, size_t populationSize, size_t maxNumGenerations, float crossoverRate, float mutationChance,
  T** (*generateRandomPopulation)(size_t, size_t), float (*getFitness)(T*), void (*mutate)(T*, size_t), bool maximization)
{
    T** population = generateRandomPopulation(chromosomeSize, populationSize);
    for (size_t i = 0; i < maxNumGenerations; i++)
    {
        float* fitnesses = getFitnesses(population, populationSize, getFitness, maximization);
        // Get upper bound
        float upperBound = 0;
        for (size_t i = 0; i < populationSize; i++)
        {
            upperBound += fitnesses[i];
        }
        T** newPopulation = new T*[populationSize];
        for (size_t j = 0; j < populationSize; j++)
        {
            std::pair<T*, T*> parents = selectIndividuals(population, populationSize, fitnesses, upperBound);
            std::pair<T*, T*> children = crossover(parents, chromosomeSize, crossoverRate, mutationChance, mutate);
            newPopulation[j] = children.first;
            //TO DO: verificar indice
            j++;
            newPopulation[j] = children.second;
        }    
        delete fitnesses;
        clear(population, populationSize);
        population = newPopulation;
    }
    return population;
}

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

std::vector<double> advance(double t, double dt, std::vector<double> u){ //to do: pass the parameters 
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

// -------------------------------------------- Root fiding problem 

float minValue;
float maxValue;
size_t populationSize;
size_t maxNumGenerations;
float mutationChance;
float crossoverRate;

//Return random value for first parameter
float getRandomValue()
{    
    float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    return r * (maxValue - minValue) + minValue;
}

float** generateRandomPopulation(size_t chromosomeSize, size_t populationSize)
{
    float** population = new float*[populationSize];
    for (size_t i = 0; i < populationSize; i++) {
        population[i] = new float[chromosomeSize];
        population[i][0] = getRandomValue();
    }
    return population;
}

void mutate(float* chromosome, size_t index){
    chromosome[index] = getRandomValue();
}

float getFitness(float* chromosome){
    std::vector<std::vector<double> > data = readCSV_to_MultidimensionalArray("data.csv"); 
    int s = data[0].size()-1;
    std::vector<double> u;
    u.reserve(s);
    u.resize(s);

    u[0] = 10; //condicao inicial da presa
    u[1] = 2; //condicao inical do predador
    r = chromosome[0]; 
    a = chromosome[1];
    m = chromosome[2];

    double errorH = 0, errorP = 0, sumexactH = 0, sumexactP = 0; 
    int i = 0;
    double tfinal = 50, dt = 0.01;
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

int main(){    
    std::cout << "Enter GA parameters: " << std::endl;
    std::cin >> minValue >> maxValue >> populationSize >> maxNumGenerations >> mutationChance >> crossoverRate;
    int nParams = 3;

    srand(time(0));
    float** solutions = geneticAlgorithm(nParams, populationSize, maxNumGenerations, crossoverRate, mutationChance,
                &generateRandomPopulation, &getFitness, &mutate, false);    

    std::map<float, float*> orderedSolutions;
    for (size_t i = 0; i < populationSize; i++){
        orderedSolutions[getFitness(solutions[i])] = solutions[i];
    }
    
    int count = 0;
    for (std::map<float, float*>::iterator it = orderedSolutions.begin(); it != orderedSolutions.end(); ++it){
        std::cout << "#" << count++ << " r = " << it->second[0] << ", a = " << it->second[1] 
        << ", m = " << it->second[2] << " (fitness = " << it->first << ")" << std::endl;
    }

    //std::cout << "Best solution (min value) = " << orderedSolutions.begin()->second[0] << std::endl; 
    //std::cout << "Best solution (max value) = " << orderedSolutions.rbegin()->second[0] << std::endl; 

    clear(solutions, populationSize);
    return 0;
}