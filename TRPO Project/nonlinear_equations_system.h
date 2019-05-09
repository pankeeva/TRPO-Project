#pragma once

#include <vector>
#include <iostream>

using namespace std;


class Equation {
public:
	//after calling this constructor, you need to call 
	//setCoefficients() and setPowers() methods
	Equation(); 
	Equation(vector<double> coefs,  vector<double> pows);
	~Equation();

	void setCoefficients(vector<double> coefs);
	void setPowers(vector<double> pows);
	void printEquation();
	double derivativeAdjective(size_t adjectiveIndex, double point);
	double functionToDerevative(double coef, double power, double point);
	int getCoefficientsSize();
private:
	vector<double> coefficients; //the last element is a free member of the equation
	vector<double> powers; 
};


class EquationSystem {
public:
	EquationSystem(vector<vector<double>> coefs, vector<vector<double>> pows);
	~EquationSystem();

	void printSystem();
	void matrixJacobi(double point);
private:
	int size; //the number of equations in the system
	Equation* equations;
};
