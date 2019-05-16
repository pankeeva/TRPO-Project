#include "nonlinear_equations_system.h"

#include <cmath>
#include <omp.h>


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
				//central derivative

	if (adjectiveIndex + 1 == doubleCoefIndex)  //Only for checkingSystem() method
	{
		return (functionToDerevativeDouble(coefficients[adjectiveIndex], (point + h)) -
			functionToDerevativeDouble(coefficients[adjectiveIndex], (point - h))) / (2 * h);
	}
	double a = coefficients[adjectiveIndex];
	double b = powers[adjectiveIndex];
	double c = point + h;
	return (functionToDerevative(coefficients[adjectiveIndex], powers[adjectiveIndex], (point + h)) -
		functionToDerevative(coefficients[adjectiveIndex], powers[adjectiveIndex], (point - h))) / (2 * h);
}

double Equation::functionToDerevative(double coef, double power, double point)
{
	return coef * std::pow(point, power);
}

double Equation::functionToDerevativeDouble(double coef, double point)
{
	return (point + coef * std::pow(point, 2));
}

EquationSystem::EquationSystem(vector<vector<double>> coefs, vector<vector<double>> pows)
{
	this->size = coefs.size();
	this->equations.clear();
	for (size_t i = 0; i < size; i++)
	{
		Equation e;
		e.setCoefficients(coefs[i]);
		e.setPowers(pows[i]);
		this->equations.push_back(e);
	}
}

vector<Equation> EquationSystem::getEquations()
{
	return this->equations;
}

EquationSystem::~EquationSystem()
{
	this->equations.clear();
}

void EquationSystem::setDoubleCoefOfEquation(vector<size_t> values)
{
	if (size != values.size())
	{
		return;
	}
#pragma omp parallel for
	for (int i = 0; i < size; i++)
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
			matr[matr.size() - 1].push_back(equations[i].derivativeAdjective(j, point[j]));
		}
	}
	return matr;
}


int Equation::getCoefficientsSize()
{
	return coefficients.size() - 1;
}

int EquationSystem::getEquationsCount()
{
	return this->size;
}

int Equation::getCoefficientsCount()
{
	return this->coefficients.size();
}

double EquationSystem::calculateInPoint(int i)
{
	return -0.5 * (3 * size + 1) - 2 * (1 + 2 * (double(i + 1) / size)
		+ std::pow(i + 1, 2) / std::pow(size, 2));
}

EquationSystem& checkingSystem(int size)
{
	vector<size_t> doubleCoefIndexes;
	vector<vector<double>> coefficients;
	vector<vector<double>> powers;

	//double sum = 0;

	for (size_t i = 0; i < size; i++)
	{
		//sum = 0;

		coefficients.push_back(vector<double>());
		powers.push_back(vector<double>());

		for (size_t j = 0; j < size; j++)
		{
			/*sum += -0.5 * (3 * size + 1) - 2 * (1 + 2 * (double(i + 1) / size)
				+ std::pow(i + 1, 2.0) / std::pow(size, 2.0));
*/
			if (i == j)
			{
				/*doubleCoefIndexes.push_back(i + 1);
				coefficients[i].push_back(2 * size);
				powers[i].push_back(2);*/

				doubleCoefIndexes.push_back(1);
				coefficients[i].push_back(1);
				powers[i].push_back(1);
			}
			else
			{
				/*coefficients[i].push_back(1);
				powers[i].push_back(1);*/

				coefficients[i].push_back(0);
				powers[i].push_back(1);
			}
		}

		//coefficients[i].push_back(sum);

		coefficients[i].push_back(1);
	}

	EquationSystem* newSystem = new EquationSystem(coefficients, powers);
	newSystem->setDoubleCoefOfEquation(doubleCoefIndexes);

	return *newSystem;
}

vector<double> Equation::getCoefficients()
{
	return this->coefficients;
}
vector<double> Equation::getPowers()
{
	return this->powers;
}