#include <iostream>
#include <vector>
#include "sol_6.h"

int main()
{
	std::cout << std::endl << "----------------- GOLDEN RATIO ------------------" << std::endl;
	golden_ratio();
	std::cout << std::endl << "----------------- INTERPOLATION -----------------" << std::endl;
	interpolation();
	std::cout << std::endl << "----------------- NEWTON METHOD -----------------" << std::endl;
	newton_method();
	std::cout << std::endl << "----------------- BRENT METHOD ------------------" << std::endl;
	brent_method();


	system("pause");
	return 0;
}
