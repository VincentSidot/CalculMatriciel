#pragma once
#include <iostream>
#include <functional>
#include "Matrice.h"

const float tol = 1e-12;

namespace Linear
{
	template<class T>
	bool isSquare(const Matrice<T> &mat)
	{
		return mat.col() == mat.ligne();
	}

	void operationLigne(Matrice<float> &m1, float lambda1,size_t i1 , float lambda2, size_t i2) // M[i1] <- l1*M[i1]+l2*M[i2]
	{
		if (i1 > m1.ligne() || i2 > m1.ligne())
		{
			return;
		}
		for (size_t j = 1; j <= m1.col(); j++)
		{
			m1(i1, j) = lambda1 * m1(i1, j) + lambda2 * m1(i2, j);
		}
	}

	float abs(float x)
	{
		if (x < 0)
			return -x;
		else
			return x;
	}


	Matrice<float> pivotGauss(Matrice<float> m)
	{
		if (isSquare(m))
		{
			size_t n = m.ligne();
			for (size_t j = 1; j <= n; j++)
			{
				for (size_t i = j; i < n; i++)
				{
					if (abs(m(i, j)) > tol)
					{
						for (size_t ni = i + 1; ni <= n; ni++)
						{
							operationLigne(m, 1, ni, -m(ni, j)/ m(i, j), i);
						}
					}
				}
			}
			return m;
		}
	}
	float det(Matrice<float> m)
	{
		if (isSquare(m))
		{
			m = pivotGauss(m);
			float rep = 1;
			for (size_t i = 1; i <= m.col(); i++)
			{
				rep *= m(i, i);
			}
			return rep;
		}
	}
	bool isInversible(const Matrice<float> &m)
	{
		return det(m) != 0;
	}


	Matrice<float> inverse(Matrice<float> m)
	{
		if (isSquare(m) && isInversible(m))
		{
			Matrice<float> rep = identity<float>(m.col());
			size_t n = m.ligne();
			for (size_t j = 1; j <= n; j++)
			{
				for (size_t i = j; i < n; i++)
				{
					if (abs(m(i, j)) > tol)
					{
						for (size_t ni = i + 1; ni <= n; ni++)
						{
							float f = -m(ni, j) / m(i, j);
							operationLigne(m, 1, ni, f, i);
							operationLigne(rep, 1, ni, f, i);
						}
					}
				}
			}
			for (size_t i = 1; i <= n; i++)
			{
				float f = m(i, i);
				for (size_t j = 1; j <= n; j++)
				{
					m(i, j) /=  f;
					rep(i, j) /= f;
				}
			}
			for (size_t j = n; j > 0; j--)
			{
				for (size_t i = n; i > 0; i--)
				{
					if (abs(m(i, j)) > tol)
					{
						for (size_t ni = i - 1; ni > 0; ni--)
						{
							float f = -m(ni, j) / m(i, j);
							operationLigne(m, 1, ni, f, i);
							operationLigne(rep, 1, ni, f, i);
						}
					}
				}
			}
			return rep;
		}
	}

	bool IsEchelone(const Matrice<float> &m1)
	{
		if (!isSquare(m1))
			return false;
		for (size_t j = 1; j <= m1.col(); j++)
		{
			for (size_t i = j + 1; i <= m1.col(); i++)
			{
				if (m1(i, j) != 0)
					return false;
			}
		}
		return true;
	}

	void ffix(Matrice<float> &m, float ttol = 1e-5)
	{
		for (size_t i = 1; i <= m.ligne(); i++)
		{
			for (size_t j = 1; j <= m.col(); j++)
			{
				if (abs(m(i, j)) < ttol)
				{
					m(i, j) = 0;
				}
			}
		}
	}
	Matrice<float> fix(Matrice<float> m, float ttol = 1e-5)
	{
		ffix(m, ttol);
		return m;
	}


	namespace calculus
	{
		float derivative(std::function<float(float)> f, float a, float h) 
		{
			//tex:
			// $$ f'(a) = \frac{f(a+h)-f(a-h)}{2h}+o(h^2) $$

			return (f(a + h) - f(a - h)) / 2 / h;

		}

		float second_derivative(std::function<float(float)> f, float a, float h)
		{
			//tex:
			//$ f''(a) = \frac{f(a+h)+f(a-h)-2f(a)}{h^2} + o(h^2)$

			return (f(a + h) + f(a - h) - 2 * f(a)) / h / h;

		}

		float trapeze(std::function<float(float)> f, float a, float b, size_t n = 10)
		{
			//tex:
			// $ \int_a^b{f(x)dx} = \frac{b-a}{n}(\frac{f(a)+f(b)}{2} + \sum_{k=1}^{n-1}{f(a+k\frac{b-a}{n})} + o(\frac{(b-a)^3}{12n^2}\sup_{[a,b]}|f''|)$
			float rep = f(a) + f(b);
			rep /= 2;
			for (size_t k = 1; k < n; k++)
			{
				rep += f(a + k * (b - a) / n);
			}
			return rep * (b - a) / n;
		}

	};

};