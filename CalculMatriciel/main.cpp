#include <iostream>
#include "LinearAlgebra.h"

float f(float x)
{
	return x;
}

int main(int argc, char** argv)
{
	float I = Linear::calculus::trapeze(f, 0, 1);
	std::cout << I << std::endl;

	std::cin.ignore().get();
	return 0;
}