#pragma once

#include <vector>
#include <iostream>

using namespace std;


class Equation {
public:
	//after calling this constructor, you need to call 
	//setCoefficients() and setPowers() methods
	Equation();
	Equation(vector<double> coefs, vector<double> pows);
	~Equation();

	void setCoefficients(vector<double> coefs);
	void setPowers(vector<double> pows);
	void setDoubleCoefIndex(size_t value);
	size_t getDoubleCoefIndex();
	void printEquation();
	double derivativeAdjective(size_t adjectiveIndex, double point);
	double functionToDerevative(double coef, double power, double point);
	double functionToDerevativeDouble(double coef, double point);    //Only for checkingSystem() method
	int getCoefficientsSize();
private:
	vector<double> coefficients; //the last element is a free member of the equation
	vector<double> powers;
	size_t doubleCoefIndex; //ex.: (x1 + x1^2). Only for checkingSystem() method
};


class EquationSystem {
public:
	EquationSystem(vector<vector<double>> coefs, vector<vector<double>> pows);
	~EquationSystem();

	void setDoubleCoefOfEquation(vector<size_t> values);
	void printSystem();
	vector<vector<double>> matrixJacobi(vector<double> point);
private:
	int size; //the number of equations in the system
	Equation* equations;
};



EquationSystem& checkingSystem(int size);
