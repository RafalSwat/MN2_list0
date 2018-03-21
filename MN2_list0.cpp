#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>

void show(double tab[5][6])
{
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			std::cout << tab[i][j] << '\t';
		}
		std::cout << std::endl;
	}
}

int el_podst(double tab[5][6]) // wyszukuje najwięszy element w macierzy (element podstawowy)
{
	double max = tab[0][0];
	int indicator = 0;
	for (int i = 0; i < 5; i++)
	{
		if (tab[i][0] > max)
		{
			max = tab[i][0];
			indicator = i;
		}
	}
	return indicator;
}

bool zad1()
{
	const double eps = 1e-12; // stała - przybliżenie zera

	const size_t N = 5;
	double A[N][N] = { {0,0,2,1,2},{0,1,0,2,-1},{1,2,0,-2,0},{0,0,0,-1,1},{0,1,-1,1,-1} }; // przykładowe dane
	double b[N] = { 1,1,-4,-2,-1 }; // przykładowe dane
	double x[N];
	double Ab[N][N + 1];
	//tworzenie rozszerzonej macierzy
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N + 1; j++)
			if (j < N)
			{
				Ab[i][j] = A[i][j];
			}
			else Ab[i][j] = b[i];
	}
	//kopiowanie tablicy
	double Ab_copy[N][N + 1] = {};
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			Ab_copy[i][j] = Ab[i][j];
			x[j] = Ab[i][N]; /// nie wiem o co tu chodzi ale działa, metoda prób i błędów
		}
	}
	//printowanie macierzy
	show(Ab_copy);
	std::cout << std::endl;
	//ustawianie elementu podstawowego
	int indicator = el_podst(Ab_copy);
	for (int j = 0; j < N + 1; j++)
	{
		std::swap(Ab_copy[0][j], Ab_copy[indicator][j]);
		std::swap(x[0], x[indicator]);
	}
	std::cout << std::endl;
	//printowanie nowej macierzy
	show(Ab_copy);
	std::cout << std::endl;

	double m = 0; // mnożnik
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N ; j++)
		{
			if (fabs(Ab_copy[i][i]) < eps) return false; //fabs zwraca wartosc absolutna
			m = Ab_copy[j][i] / Ab_copy[i][i];
			for (int k = i + 1; k <= N; k++)
			{
				Ab_copy[j][k] = Ab_copy[j][k] - m * Ab_copy[i][k];
			}
		}
	}
	//etap wyliczania zmiennych
	double s = 0; // zlicza sumę iloczynów
	for (int i = N-1; i >= 0; i--)
	{
		s = Ab_copy[i][N]; // do zmiennej s trafia wektor b (ostatnia kolumna macierzy)
		for (int j = N; j >= i + 1; j--)
		{
			s = s - Ab_copy[i][j] * x[j];
		}
		if (fabs(Ab_copy[i][i]) < eps) return false;
		x[i] = s / Ab_copy[i][i];
		std::cout << 'x' << i << "= " << x[i] << std::endl;
	}
	return true;
}


int main()
{
	zad1();

	system("pause");
	return 0;
}
