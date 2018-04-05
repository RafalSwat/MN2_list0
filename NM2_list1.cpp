#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>

void zad1()
{
	std::ofstream file_1;
	file_1.open("data_1.csv", std::ios::out);

	int M = 200; // dodatnie przyrosty w czasie
	double h = 0.1; // krok w przestrzeni
	double k = 0.005125; // krok w czsie
	double n = 1 / h - 1; // liczba kroków w przestrzeni
	double s = k / (h * h); 

	double w[10] = {};	

	//ustawianie wartosci temp dla t = 0. ( po położeniu od 0 do 1 co h)
	for (int i = 0; i < n + 1; i++)
	{
		w[i] = sin(M_PI * i * h);
	}

	double t = 0;

	for (int j = 1; j < M; j++)
	{
		t = j*k; /// inkremenacja czasu (po czasie od 0 do 1 co k, *powinno byc 195.. ~ 200 czyli M)
		double v[10] = {}; // tab. przechowująca wartości temp. dla t != 0
		for (int i = 1; i < n; i++)
		{
			v[i] = s * w[i - 1] + (1 - 2 * s) * w[i] + s * w[i + 1];
		}
		double solution[10] = {};
		for (int i = 0; i < n + 1; i++)
		{
			w[i] = v[i]; // przepisanie wartosci temp.
			file_1 << w[i] << ',';
			solution[i] = pow(M_E, -(M_PI*M_PI)* t) * sin(M_PI * i * h); // no to cos nie wyszlo
			std::cout << w[i] << "-->" << solution[i] << std::endl;

		}
		file_1 << std::endl;
	}
}
void tri(double n,double a[10], double c[10], double d[10], double v[10])// funkcja rozwiązująca ukłąd równań (zad2)
{
	n = 20;
	double b[10];
	for (int k = 0; k < 10; k++)
	{
		b[k] = v[k];
	}
	for (int i = 2; i < n; i++)
	{
		d[i] = d[i] - (a[i - 1] / d[i - 1]) * c[i - 1];
		b[i] = b[i] - (a[i - 1] / d[i - 1]) * b[i - 1];
	}
}
void zad2()
{
	std::fstream file_2;
	file_2.open("data_1.csv", std::ios::out);
	
	int M = 50; // dodatnie przyrosty w czasie
	double h = 0.05; // krok w przestrzeni (x należy do (0,1), 1/n=h)
	double k = 0.05; /// krok w czsie ( przy M=50, k powinno być 0.02 bo 1/M, jesli istotnie ma być to kwadrat jednostkowy)
	double n = 20; // liczba kroków w przestrzeni
	double s = k / (h * h);

	double v[10] = {};
	//ustawianie wartosci temp dla t = 0. ( po położeniu od 0 do 1 co h)
	for (int i = 1; i < n - 1; i++)
	{
		v[i] = (h * i + h * i * h * i) * pow(M_PI, h * i); 
	}
	double t = 0;
	//ustawiawianie wartości współczynników c, a(tri - macierz-układ równań)
	double a[10] = {};
	double c[10] = {};
	for (int i = 1; i < n - 1; i++)
	{
		a[i] = -s;
		c[i] = -s;
	}
	//ustawianie wartości współczynika d(diagonal) (call tri - rozwiąrz układ równań)
	double d[10] = {};
	for (int j = 1; j < M; j++)
	{
		for (int i = 1; i < n; i++)
		{
			d[i] = 1 + 2 * s;

		}
	}

	file_2.close();
}


int main()
{

	zad1();
	//zad2();
	system("pause");
	return 0;
}
