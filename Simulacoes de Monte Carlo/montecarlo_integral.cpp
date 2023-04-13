#include <ctime>
#include <iostream>
#include <random>
using namespace std;

int main(){
	/*
    	Integral de x ^ 2 por monte carlo no intervalo de 
    	domínio com x de 0 a 2 = 8 / 3.
	max(x ^ 2) em(0, 2) = 4.
	Definindo retângulo com área = (2 - 0)x(4 - 0) = 8 base x altura
	Definindo f(x) = x ^ 2
	*/

	//ponto de máximo x = 2 e valor máximo H = 2 ^ 2 = 4
	int  H = 4;

	//área do retângulo
	int a = 0;
	int b = 2;
	int A = (b - a) * H;

	//declarando o gerador de números pseudo-aleatórios
	mt19937 mt(time(NULL));

	//definindo contador
	int count = 0;

	int N = 1000000;
	for (unsigned int i = 0; i < N; ++i) {
		//gera x uniforme entre 0 e 2

		uniform_real_distribution<float> linear_x(0.f, 2.f);
		double x = linear_x(mt);


		//gera y uniforme entre 0 e 4
		uniform_real_distribution<float> linear_y(0.f, 4.f);
		double y = linear_y(mt);

		//calculando z = f(x) = x ^ 2 para cada x gerado
		double f = x * x;


		//verificando se a coordenada está embaixo de f(x) e somando a contagem de pontos
		if (y <= f) {
			++count;
		}
	}

	double prop = (double)count / N;
	double integral = A * prop;
	cout << "O valor da integral de f(x)=x^2 = " << integral << endl;
    
    return 0;
}