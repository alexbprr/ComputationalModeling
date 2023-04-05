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
        for (size_t j = 0; j < populationSize; )
        {
            std::pair<T*, T*> parents = selectIndividuals(population, populationSize, fitnesses, upperBound);
            std::pair<T*, T*> children = crossover(parents, chromosomeSize, crossoverRate, mutationChance, mutate);
            newPopulation[j] = children.first;
            j++;
            newPopulation[j] = children.second;
            j++;
        }    
        delete fitnesses;
        clear(population, populationSize);
        population = newPopulation;
    }
    return population;
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
    float x = chromosome[0];
    return std::pow(x,2);
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
        std::cout << "#" << count++ << " x = " << it->second[0] << " (fitness = " << it->first << ")" << std::endl;
    }
    std::cout << "Best solution (min value) = " << orderedSolutions.begin()->second[0] << std::endl; 
    std::cout << "Best solution (max value) = " << orderedSolutions.rbegin()->second[0] << std::endl; 

    clear(solutions, populationSize);
    return 0;
}