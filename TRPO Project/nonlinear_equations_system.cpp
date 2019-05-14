#include "nonlinear_equations_system.h"

#include <cmath>


Equation::Equation() :
	doubleCoefIndex(0)
{

}

Equation::Equation(vector<double> coefs, vector<double> pows) :
	coefficients(coefs),
	powers(pows),
	doubleCoefIndex(0)
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

void Equation::setDoubleCoefIndex(size_t value)
{
	this->doubleCoefIndex = value;
}

size_t Equation::getDoubleCoefIndex()
{
	return doubleCoefIndex;
}

void Equation::printEquation()
{
	for (size_t i = 0; i < coefficients.size() - 1; i++)
	{
		if (coefficients[i])
		{
			if (coefficients[i] != 1)
			{
				cout << coefficients[i];
			}

			cout << "x" << i + 1;

			if (powers[i] > 1)
			{
				cout << "^" << powers[i];

				if (doubleCoefIndex && (i + 1) == doubleCoefIndex)  //Only for checkingSystem() method
					cout << " + x" << i + 1;
			}
			
			cout << " + ";
		}
	}

	if (coefficients[coefficients.size() - 1])
	{
		cout << " = (" << coefficients[coefficients.size() - 1] * (-1) << ")";
	}
	else cout << " = 0";

	cout << endl;
}

double Equation::derivativeAdjective(size_t adjectiveIndex, double point)
{
	double h = 0.000000001; //step to calculate derivative

							//approximately calculate the first derivative 
							//left way

	if (adjectiveIndex + 1 == doubleCoefIndex)  //Only for checkingSystem() method
	{
		return (functionToDerevativeDouble(coefficients[adjectiveIndex], point) -
			functionToDerevativeDouble(coefficients[adjectiveIndex], (point - h))) / h;
	}

	return (functionToDerevative(coefficients[adjectiveIndex], powers[adjectiveIndex], point) -
		functionToDerevative(coefficients[adjectiveIndex], powers[adjectiveIndex], (point - h))) / h;
}

double Equation::functionToDerevative(double coef, double power, double point)
{
	return coef * std::pow(point, power);
}

double Equation::functionToDerevativeDouble(double coef, double point)
{
	return (point + coef * std::pow(point, 2));
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

void EquationSystem::setDoubleCoefOfEquation(vector<size_t> values)
{
	if (size != values.size())
	{
		return;
	}

	for (size_t i = 0; i < size; i++)
	{
		equations[i].setDoubleCoefIndex(values[i]);
	}
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
			matr[matr.size()-1].push_back(equations[i].derivativeAdjective(j, point[j]));
		}
	}
	return matr;
}


int Equation::getCoefficientsSize()
{
	return coefficients.size() - 1;
}



EquationSystem& checkingSystem(int size)
{
	vector<size_t> doubleCoefIndexes;
	vector<vector<double>> coefficients;
	vector<vector<double>> powers;

	double sum = 0;

	for (size_t i = 0; i < size; i++)
	{
		sum = 0;

		coefficients.push_back(vector<double>());
		powers.push_back(vector<double>());

		for (size_t j = 0; j < size; j++)
		{
			sum += -0.5 * (3 * size + 1) - 2 * (1 + 2 * (double (i + 1) / size) 
				+ std::pow(i + 1, 2) / std::pow(size, 2));

			if (i == j)
			{
				doubleCoefIndexes.push_back(i + 1);
				coefficients[i].push_back(2 * size);
				powers[i].push_back(2);
			}
			else
			{
				coefficients[i].push_back(1);
				powers[i].push_back(1);
			}
		}

		coefficients[i].push_back(sum);
	}

	EquationSystem* newSystem = new EquationSystem(coefficients, powers);
	newSystem->setDoubleCoefOfEquation(doubleCoefIndexes);

	return *newSystem;
}