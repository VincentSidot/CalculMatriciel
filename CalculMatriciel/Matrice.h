#pragma once
#include <sstream>
#include <iostream>
#include <string>

template<class T>
class Matrice
{
public:
	
	Matrice(size_t col, size_t ligne) : m_col{ col }, m_ligne{ ligne }
	{

		m_data = new T[m_col*m_ligne]{ 0 };
	}
	Matrice(T** mat, size_t col, size_t ligne) : m_col{ col }, m_ligne{ ligne }
	{
		m_data = new T[m_col*m_ligne]{ 0 };
		for (size_t i = 0; i < m_ligne; i++)
		{
			//tex:
			//Index : $$ M(i,j) --> M[i + j*m_ligne] $$
			for (size_t j = 0; i < m_col; j++)
			{
				m_data[i + j * m_ligne] = mat[i][j];
			}
		}
	}

	Matrice(T* mat, size_t col, size_t ligne) : m_col{ col }, m_ligne{ ligne }
	{
		m_data = new T[m_col*m_ligne]{ 0 };
		for (size_t i = 0; i < m_col*m_ligne; i++)
		{
			m_data[i] = mat[i];
		}
	}

	Matrice(const Matrice<T> &othr) : m_col{ othr.col() }, m_ligne{ othr.ligne() }
	{
		m_data = new T[m_col*m_ligne]{ 0 };
		for (size_t i = 0; i < m_ligne; i++)
		{
			//tex:
			//Index : $$ M(i,j) --> M[i + j*\text{m_ligne}] $$
			for (size_t j = 0; j < m_col; j++)
			{
				m_data[i + j * m_ligne] = othr(i + 1, j + 1);
			}
		}
	}

	~Matrice()
	{
		delete[] m_data;
	}

	size_t ligne() const
	{
		return m_ligne;
	}
	size_t col() const
	{
		return m_col;
	}
	T* data() const
	{
		return m_data;
	}


	T& operator()(size_t i, size_t j)
	{
		if (i > 0 && i <= m_ligne && j > 0 && j <= m_col)
			return m_data[i - 1 + m_ligne * (j - 1)];
		else
			throw "Bad arguments";
	}
	T operator()(size_t i, size_t j) const
	{
		if (i > 0 && i <= m_ligne && j > 0 && j <= m_col)
			return m_data[i - 1 + m_ligne * (j - 1)];
		else
			throw "Bad arguments";
	}
	Matrice<T>& operator=(const Matrice<T> mat)
	{
		this->~Matrice();
		m_col = mat.col();
		m_ligne = mat.ligne();
		m_data = new T[m_col*m_ligne];
		for (size_t i = 0; i < m_col*m_ligne; i++)
		{
			m_data[i] = mat.data()[i];
		}
		return *this;
	}
	bool operator==(const Matrice<T> &mat)
	{
		if (m_col != mat.col() || m_ligne != mat.ligne())
			return false;
		for (size_t i = 0; i < m_col*m_ligne; i++)
		{
			if (m_data[i] != mat.data[i])
				return false;
		}
		return true;
	}
	void operator+=(const Matrice<T> &mat)
	{
		if (mat.col() != m_col || mat.ligne() != m_ligne)
		{
			throw "Bad size";
		}
		for (size_t i = 0; i < m_ligne*m_col; i++)
		{
			m_data[i] += mat.data()[i];
		}
	}
	void operator*=(const T &scalar)
	{
		for (size_t i = 0; i < m_ligne*m_col; i++)
		{
			m_data[i] *= scalar;
		}
	}
	
	std::string disp() const
	{
		std::stringstream ss;
		for (size_t i = 1; i <= this->m_ligne; i++)
		{
			ss << "| ";
			for (size_t j = 1; j <= this->m_col; j++)
			{
				ss << this->operator()(i, j) << " ";
			}
			ss << "|" << std::endl;
		}
		return ss.str();
	}

protected:
	
	T* m_data;
	size_t m_col;
	size_t m_ligne;


};

template<class T>
Matrice<T> operator*(const T & scalar, Matrice<T> m2)
{
	m2 *= scalar;
	return m2;
}
template<class T>
Matrice<T> operator*(Matrice<T> m1, const Matrice<T> &m2)
{
	if (m1.col() == m2.ligne())
	{
		Matrice<T> mat(m2.col(), m1.ligne());
		for (size_t i = 1; i <= mat.ligne(); i++)
		{
			for (size_t j = 1; j <= mat.col(); j++)
			{
				mat(i, j) = T(0);
				for (size_t k = 1; k <= m1.col(); k++)
				{
					mat(i, j) += m1(i,k)*m2(k,j);
				}
			}
		}
		return mat;
	}
	else
	{
		throw "Bad size";
	}
}
template<class T>
Matrice<T> operator*(Matrice<T> m1, const T &scalar)
{
	m1 *= scalar;
	return m1;
}
template<class T>
Matrice<T> operator*(const T &scalar, const Matrice<T> &m2)
{
	return operator*(m2, scalar);
}
template<class T>
Matrice<T> operator+(Matrice<T> m1, const Matrice<T> &m2)
{
	m1 += m2;
	return m1;
}

template<class T>
Matrice<T> transpose(const Matrice<T>& m)
{
	Matrice<T> mat(m.ligne(), m.col());
	for (size_t i = 1; i <= mat.ligne(); i++)
	{
		for (size_t j = 1; j <= mat.col(); j++)
		{
			mat(i, j) = m(j, i);
		}
	}
	return mat;
}


template<class T>
std::ostream& operator<<(std::ostream &os, const Matrice<T> &m1)
{
	os << m1.disp();
	return os;
}


// Special Matrice

template<class T>
Matrice<T> identity(size_t n)
{
	Matrice<T> mat(n, n);
	for (size_t i = 1; i <= n; i++)
	{
		mat(i, i) = 1;
	}
	return mat;
}

template<class T>
Matrice<T> anti_identity(size_t n)
{
	Matrice<T> mat(n, n);
	for (size_t i = 1; i <= n; i++)
	{
		mat(n + 1 - i, i) = 1;
	}
	return mat;
}

template<class T>
Matrice<T> M(size_t i,size_t j,size_t ligne,size_t col)
{
	Matrice<T> mat(ligne, col);
	mat(i, j) = 1;
	return mat;
}

template<class T>
Matrice<T> M(size_t i, size_t j, size_t n)
{
	Matrice<T> mat(n,n);
	mat(i, j) = 1;
	return mat;
}