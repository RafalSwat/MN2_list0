#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

void metoda_roznic_skonczonych(int &n, int &M, int &f);
void metoda_Galerkina();


double g(double x, double y);
double f(double x, double y);
double u1(double x, double y);
double u2(double x, double y);
double u3(double x, double y);
double u4(double x, double y);
int el_podst(double ** tab); // wyszukuje najwięszy element w macierzy (element podstawowy)
bool eliminacja_gaussa();

/*double g(double x, double y)
{
return 4*x*y*(x*x - y*y);
}*/
double g(double x, double y)
{
	return pow(10, -4)*sin(3 * M_PI*x)*sin(3 * M_PI*y);
}
double f(double x, double y)
{
	return x * x + y * y;
}
double u1(double x, double y)
{
	return 1;
}
double u2(double x, double y)
{
	return x * x - y * y;
}
double u3(double x, double y)
{
	return x * x*x*x + 6 * x*x*y*y + y * y*y*y;
}
double u4(double x, double y)
{
	return std::pow(x, 6) - 15 * std::pow(x, 4)*y*y + 15 * x*x*std::pow(y, 4) - std::pow(y, 6);
}
int el_podst(double ** tab) // wyszukuje najwięszy element w macierzy (element podstawowy)
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
bool eliminacja_gaussa()
{
	const double eps = 1e-12; // stała - przybliżenie zera
	const size_t N = 5;

	double A[N][N] = { { 1,0,0,0,0 },{ 1,1,1,1,1 },{ 1,3,9,27,81 },{ 1,5,25,125,625 },{ 1,6,36,216,1296 } };// dane zad3
	double b[N] = { 1,1,-4,-2,-1 }; // dane zad3
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
		for (int j = i + 1; j < N; j++)
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
	for (int i = N - 1; i >= 0; i--)
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

void metoda_roznic_skonczonych(int &n, int &M, int &f)
{
	std::ofstream file;
	file.open("data_zad1.csv", std::ios::out);
	double h = 0.001;
	double k = 0.001;
	double **matrix = new double *[n + 1];                // definition of base to hold results

	for (int i = 0; i < n + 1; i++) matrix[i] = new double[n + 1];
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			if (i == 0) matrix[i][j] = g(0, j*k);
			else if (j == 0) matrix[i][j] = g(i*h, 0);
			else if (i == n) matrix[i][j] = g(1, j*k);
			else if (j == n) matrix[i][j] = g(i*h, 1);
			else matrix[i][j] = 0;
		}
	}

	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			std::cout << matrix[i][j] << '\t';
		}
		std::cout << std::endl;
	}

	for (int k = 1; k < M; k++)
	{
		for (int j = 1; j < n; j++)
		{
			for (int i = 1; i < n; i++)
			{
				matrix[i][j] = (matrix[i - 1][j] + matrix[i + 1][j] + matrix[i][j - 1] + matrix[i][j + 1]) / 4;
			}
		}
	}

	for (int i = 0; i < n + 1; i += f)
	{
		for (int j = 0; j < n + 1; j += f)
		{
			std::cout << matrix[i][j] << '\t';
			//file << matrix[i][j] << ',';
		}
		std::cout << std::endl;
		//file << '\n';
	}

	for (int i = 0; i < n + 1; i++) delete matrix[i];
	delete[] matrix;
}

void metoda_Galerkina()
{
	double x1 = 0, x2 = 1, x3 = 1, x4 = 1;
	double y1 = 2, y2 = 0, y3 = 1, y4 = 2;

	size_t size = 4;

	double ** u = new double*[size];
	double * solution = new double[size];
	double * c = new double[size];

	for (size_t i = 0; i < size; i++)
	{
		u[i] = new double[size];
		for (size_t j = 0; j < size; j++)
		{
			if (j == 0) u[i][j] = 1;
			else if (j == 1 && i == 0) u[i][j] = u2(x1, y1);
			else if (j == 1 && i == 1) u[i][j] = u2(x2, y2);
			else if (j == 1 && i == 2) u[i][j] = u2(x3, y3);
			else if (j == 1 && i == 3) u[i][j] = u2(x4, y4);
			else if (j == 2 && i == 0) u[i][j] = u3(x1, y1);
			else if (j == 2 && i == 1) u[i][j] = u3(x2, y2);
			else if (j == 2 && i == 2) u[i][j] = u3(x3, y3);
			else if (j == 2 && i == 3) u[i][j] = u3(x4, y4);
			else if (j == 3 && i == 0) u[i][j] = u4(x1, y1);
			else if (j == 3 && i == 1) u[i][j] = u4(x2, y2);
			else if (j == 3 && i == 2) u[i][j] = u4(x3, y3);
			else if (j == 3 && i == 3) u[i][j] = u4(x4, y4);
		}
	}
	// Przypisanie wartości wekotrowi z f (tylko 4 wart, nie opłaca się robić poprzez funkcje)
	solution[0] = 4;
	solution[1] = 1;
	solution[2] = 2;
	solution[3] = 5;

	for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			std::cout << u[i][j] << '\t';
		}
		std::cout << " C" << i << '\t' << solution[i];
		std::cout << std::endl;
	}




	for (int i = 0; i < size; i++) delete[] u[i];
	delete[] u;
}





int main()
{
	//metoda_roznic_skonczonych();
	metoda_Galerkina();

	system("pause");
	return 0;
}
