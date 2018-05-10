#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "sol4.h"

void exercise1()
{
	std::ofstream file;
	file.open("./NM2_list4_ex1.csv", std::ios::out);
	double u[50][100];
	double x = 1;
	double t = 0;
	double h = 0.02;
	double k = 0.01;
	double alpha = 2;
	double s = double((alpha*k) / (2 * h));

	for (size_t i = 0; i < 50; i++)
	{
		u[i][0] = std::pow(1 - x, 2);
		x += h;
	}
	for (size_t i = 1; i < 49; i++)
	{
		for (size_t j = 1; j < 99; j++)
		{
			///u[i][j] = (2 * std::pow(s, 2) + s) * u[i + 1][j] + (1 - 4 * std::pow(s, 2) * u[i][j] + (2 * std::pow(s, 2) - s) * u[i - 1][j]);
			u[i][j] = u[i][j - 1] - (k / (2 * h)) * alpha * (u[i + 1][j - 1] - u[i - 1][j - 1]) + (std::pow(k, 2) / (2 * std::pow(h, 2))) * alpha * alpha
				* (u[i + 1][j - 1] - 2 * u[i][j - 1] + u[i - 1][j - 1]);

		}
	}
	for (size_t i = 0; i < 50; i++)
	{
		for (size_t j = 0; j < 100; j++)
		{
			//std::cout << u[i][j] << '\t' ;
			file << u[i][j] << '\t';
		}
		//std::cout << std::endl;
		file << std::endl;
	}
}
void exercise2()
{
	std::ofstream file;
	file.open("./NM2_list4_ex2.csv", std::ios::out);
	
	double m = 6;
	double h = 0.5;
	double k = 3;
	double n = 1;
	
	size_t size = 100;

	double * v = new double[size];
	double * w = new double[size];
	
	for (int i = 0; i < size; i++)
	{
		w[i] = 0;
		v[i] = 0;
	}

	v[0] = 0;
	v[1] = -cos(h) / 8;
	v[2] = 0;

	for (int i = 1; i < m; i++)
	{
		h = h / 2;
		n = 2 * n + 1;
		for (int j = 0; j < (n + 1) / 2; j++)
		{
			w[2 * j] = v[j];
		}
		for (int j = 1; j < (n + 1) / 2; j++)
		{
			w[2 * j - 1] = (v[j - 1] + v[j]) / 2;
		}
		for (int j = 0; j < n + 1; j++)
		{
			v[j] = w[j];
		}
		for (int p = 1; p < k; p++)
		{
			for (int j = 1; j < n; j++)
			{
				v[j] = (v[j - 1] + v[j + 1] - h*h*(cos(j*h))) / 2;
			}
		}
		for (int i = 0; i < sizeof(v); i++)
		{
			std::cout << v[i] << std::endl;
			file << v[i] << std::endl;
		}
		std::cout << "h = " << h << std::endl;
	}
	
	delete [] v;
	delete [] w;
}
