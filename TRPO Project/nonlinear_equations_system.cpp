#include "nonlinear_equations_system.h"

#include <cmath>


Equation::Equation()
{

}

Equation::Equation(vector<double> coefs, vector<double> pows) :
	coefficients(coefs),
	powers(pows)
{

}

Equation::~Equation()
{
	coefficients.erase(coefficients.begin(), coefficients.end());
	powers.erase(powers.begin(), powers.end());
}

void Equation::setCoefficients(vector<double> coefs)
{
	this->coefficients = coefs;
}

void Equation::setPowers(vector<double> pows)
{
	this->powers = pows;
}

void Equation::printEquation()
{
	for (size_t i = 0; i < coefficients.size() - 1; i++)
	{
		if (coefficients[i])
		{
			if (powers[i] > 1)
			{
				cout << coefficients[i] << "x" << i + 1 << "^" << powers[i] << "  ";
			}
			else
			{
				cout << coefficients[i] << "x" << i + 1 << "  ";
			}
		}
	}

	if (coefficients[coefficients.size() - 1])
	{
		cout << " = (" << coefficients[coefficients.size() - 1]*(-1) << ")";
	}
	else cout << " = 0";

	cout << endl;
}

double Equation::derivativeAdjective(size_t adjectiveIndex, double point)
{
	double h = 0.000000001; //step to calculate derivative

							//approximately calculate the first derivative 
							//left way
	return (functionToDerevative(coefficients[adjectiveIndex], powers[adjectiveIndex], point) -
		functionToDerevative(coefficients[adjectiveIndex], powers[adjectiveIndex], (point - h))) / h;
}

double Equation::functionToDerevative(double coef, double power, double point)
{
	return coef * std::pow(point, power);
}

EquationSystem::EquationSystem(vector<vector<double>> coefs, vector<vector<double>> pows) :
	size(coefs.size()),
	equations(new Equation[size])
{
	for (size_t i = 0; i < size; i++)
	{
		equations[i].setCoefficients(coefs[i]);
		equations[i].setPowers(pows[i]);
	}
}

EquationSystem::~EquationSystem()
{
	delete[] equations;
}


void EquationSystem::printSystem()
{
	for (size_t i = 0; i < size; i++)
	{
		equations[i].printEquation();
	}

	cout << endl;
}

vector<vector<double>> EquationSystem::matrixJacobi(vector<double> point)
{
	vector<vector<double>> matr;
	for (size_t i = 0; i < size; i++)
	{
		matr.push_back(vector<double>());
		size_t tempSize = equations[i].getCoefficientsSize();

		for (size_t j = 0; j < tempSize; j++)
		{
			matr[matr.size()-1].push_back(equations[i].derivativeAdjective(j, point[i]));
		}
	}
	return matr;
}


int Equation::getCoefficientsSize()
{
	return coefficients.size() - 1;
}
