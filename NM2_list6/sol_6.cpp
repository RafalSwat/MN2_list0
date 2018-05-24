#include <iostream>
#include <vector>
#include <cmath>
#include "sol_6.h"

double f(double x)
{
	return 2 * std::sin(x) - (x*x) / 10;
}
double df(double x)
{
	double dx = 10e-8;
	return (f(x + dx) - f(x)) / dx;
}
double d2f(double x)
{
	double dx = 10e-8;
	return ((f(x + dx) - 2 * f(x) + f(x - dx)) / (dx*dx));
}
void golden_ratio()
{
	double x_opt;
	double x_low = 0;
	double x_high = 4;
	double eps_current;
	double eps = 10e-6;
	unsigned int max_iter = 100;

	double R = (sqrt(5) - 1) / 2;
	double xl = x_low;
	double xu = x_high;
	double d = R * (xu - xl);
	double x1 = xl + d;
	double x2 = xu - d;
	double fx;
	double f1 = f(x1);
	double f2 = f(x2);

	if (f1 > f2)
	{
		x_opt = x1;
		fx = f1;
	}
	else
	{
		x_opt = x2;
		fx = f2;
	}
	for (int i = 1; i < max_iter; i++)
	{
		d = R * d;
		if (f1 > f2)
		{
			xl = x2;
			x2 = x1;
			x1 = xl + d;
			f2 = f1;
			f1 = f(x1);
			x_opt = x1;
			fx = f1;
		}
		else
		{
			xu = x1;
			x1 = x2;
			x2 = xu - d;
			f1 = f2;
			f2 = f(x2);
			x_opt = x2;
			fx = f2;
		}
		if (x_opt != 0)
		{
			eps_current = (1 - R) * std::abs((xu - xl) / x_opt);
		}
		if (eps_current <= eps)
		{
			std::cout << "numbers of iterations:         " << i << std::endl;
			std::cout << "max for x:                     " << x_opt << std::endl;
			std::cout << "value of function on extremum: " << fx << std::endl;
			break;
		}
	}
}
void interpolation()
{
	unsigned int max_iter = 1000;
	unsigned int iteration = 0;
	double eps = 10e-6;
	double x0 = 0;
	double x1 = 2;
	double x2 = 4;
	double x3 = 999;
	double previous_max = 999;
	double current_max = 0;

	while (abs(current_max - previous_max) > eps && iteration < max_iter)
	{
		x3 = (f(x0)*(x1*x1 - x2 * x2) + f(x1)*(x2*x2 - x0 * x0) + f(x2) * (x0*x0 - x1 * x1)) / (2 * f(x0) * (x1 - x2) + 2 * f(x1) * (x2 - x0) + 2 * f(x2)*(x0 - x1));
		x0 = x1;
		x1 = x2;
		x2 = x3;
		iteration++;
		if (f(x3) > current_max)
		{
			previous_max = current_max;
			current_max = f(x3);
		}
	}
	std::cout << "numbers of iterations:         " << iteration << std::endl;
	std::cout << "max for x:                     " << x3 << std::endl;
	std::cout << "value of function on extremum: " << current_max << std::endl;
}
void newton_method()
{
	int M = 100;
	int iteration = 0;
	double x0 = 2;
	double x1 = 0;
	double eps = 10e-8;
	double delta = 10e-4;
	double v = df(x0);

	for (int k = 1; k < M; k++)
	{
		if (abs(v) < eps)
		{
			break;
		}
		else
		{
			x1 = x0 - (v / d2f(x0));
			v = df(x1);

			if (abs(x1 - x0) < delta || (abs(v) < eps))
			{
				break;
			}
			x0 = x1;
			iteration = k;
		}
	}
	std::cout << "numbers of iterations:         " << iteration << std::endl;
	std::cout << "max for x:                     " << x1 << std::endl;
	std::cout << "value of function on extremum: " << f(x1) << std::endl;
}
