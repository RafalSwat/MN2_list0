#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <float.h>

bool eliminacja_gaussa();
int el_podst(double tab[5][6]);
void metoda_newtona2();
void metoda_siecznych();
double metoda_brenta(double a, double b, double t);
void Runge_Kutta();

void show(double tab[5][6]);
double f1(double &x);
double df1(double &x);
double f(double x);
double gv(double teta, double t, double v);

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
bool eliminacja_gaussa()
{
	const double eps = 1e-12; // stała - przybliżenie zera

	const size_t N = 5;
	///double A[N][N] = { {0,0,2,1,2},{0,1,0,2,-1},{1,2,0,-2,0},{0,0,0,-1,1},{0,1,-1,1,-1} }; // dane zad2
	///double b[N] = { 1,1,-4,-2,-1 }; // dane zad2

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
double f1(double &x)
{
	return tan(x) - x;
}
double df1(double &x)
{
	return (1 / (pow(cos(x), 2))) - 1;
}
void metoda_newtona2()
{
	double next = 4.5;
	double prev;
	double eps = 1e-8;
	do
	{
		prev = next;
		next = prev - f1(prev) / df1(prev);
		std::cout << prev << '\t' << next << '\t' << fabs(next - prev) << std::endl;
	} while (fabs(next - prev) > eps);
}
void metoda_siecznych()
{
	double x0;
	double x1 = 4.6;
	double x2 = 4.4;
	double fx0;
	double fx1 = f1(x1);
	double fx2 = f1(x2);
	double eps = 1e-8;
	int i = 64;

	while ((i > 0) && (abs(x2 - x1) > eps))
	{
		if (fabs(fx1 - fx2) < eps)
		{
			std::cout << "zle pkt startowe" << std::endl;
			break;
		}
		x0 = x1 - fx1*((x2 - x1) / (fx1 - fx2));
		fx0 = f1(x0);
		if (fx0 < eps)
		{
			std::cout << '!' << x0 << std::endl;
			break;
		}
		x2 = x1;
		fx2 = fx1;
		x1 = x0;
		fx1 = fx0;
		i = i - 1;
		if (i > 0) std::cout << "przekroczony limit obiegow" << std::endl;
		std::cout << x0 << std::endl;
	}
}
double f(double x)
{
	return tan(x) - x;
}
double metoda_brenta(double a, double b, double t)
{
	double x, d, e, m, p, q, r, tol, t2, u, v, w, fu, fv, fw, fx;
	double eps = sqrt(DBL_EPSILON);
	const double c = (3 - sqrt(5)) / 2;
	v = w = x = a + c*(b - a);
	e = 0;
	fv = fw = fx = f(x);
	//main loop	
	while (true)
	{
		m = 0.5*(a + b);
		tol = eps*fabs(x) + t;
		t2 = 2 * tol;
		//Check stopping criterion
		if (fabs(x - m) <= t2 - 0.5*(b - a)) break;
		p = q = r = 0;
		if (fabs(e) > tol)
		{//fit parabola
			r = (x - w)*(fx - fv);
			q = (x - v)*(fx - fw);
			p = (x - v)*q - (x - w)*r;
			q = 2 * (q - r);
			if (q > 0)p = -p; else q = -q;
			r = e; e = d;
		}
		if (fabs(p) < fabs(0.5*q*r) && p > q*(a - x) && p < q*(b - x))
		{
			// a "parabolic interpolation" step
			d = p / q;
			u = x + d;
			// f must not be evaluated too close to a or b
			if (u - a < t2 || b - u < t2)
			{
				d = x < m ? tol : -tol;
			}
		}
		else
		{  // a "golden section" step
			e = (x < m ? b : a) - x;
			d = c*e;
		}
		// f must not be evaluated too close to x		
		u = x + (fabs(d) >= tol ? d : (d > 0 ? tol : -tol));
		fu = f(u);
		// update a,b,v,w and x
		if (fu <= fx)
		{
			if (u < x)b = x; else a = x;
			v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
		}
		else
		{
			if (u < x) a = u; else b = u;
			if (fu <= fw || w == x)
			{
				v = w; fv = fw; w = u; fw = fu;
			}
			else if (fu <= fv || v == x || v == w)
			{
				v = u; fv = fu;
			}
		}
	}
	return x;
}
double gv(double teta, double t, double v) // func. przyspieszenia
{
	double A = 0.5;
	double Q = 2;
	return (A*cos(t)) - v / Q - sin(teta);

}
void Runge_Kutta()
{
	std::ofstream file;
	file.open("./data.csv", std::ios::out);

	double kv1, kv2, kv3, kv4;
	double kx1, kx2, kx3, kx4;

	std::vector<double> v; v.push_back(0);
	std::vector<double> x; x.push_back(0.1);
	std::vector<double> t; t.push_back(0);

	double h = 0.015;
	size_t n = 10000;

	for (size_t i = 0; i < n; i++)
	{
		kx1 = v[i];
		kv1 = gv(x[i], t[i], v[i]);

		kx2 = v[i] + kx1 / 2;
		kv2 = gv(x[i] + kv1 / 2, t[i] + h / 2, v[i]);

		kx3 = v[i] + kx2 / 2;
		kv3 = gv(x[i] + kv2 / 2, t[i] + h / 2, v[i]);

		kx4 = v[i] + kx3;
		kv4 = gv(x[i] + kv3, t[i] + h, v[i]);

		v.push_back(v[i] + (h / 6) * (kv1 + 2 * kv2 + 2 * kv3 + kv4));
		x.push_back(x[i] + (h / 6) * (kx1 + 2 * kx2 + 2 * kx3 + kx4));

		t.push_back(t[i] + h);

		std::cout << ", v = " << v[i] << ", x = " << x[i] << ", t = " << t[i] << std::endl;
		file << x[i] << ' ' << t[i] << std::endl;
	}
	file.close();
}

int main()
{
	//eliminacja_gaussa();
	//metoda_siecznych();
	//double res = metoda_brenta(4.2, 4.6, 1e-16);
	//std::cout << res << std::endl;
	system("pause");
	return 0;
}
