#pragma once

#include "nonlinear_equations_system.h"

class Result {
public:
	bool is_solution;
	vector<double> matrix;
};

class Newton {
private:
	double Det(double**arr, int size, int p);

	void Print(vector<vector<double>> p);
	vector<double> EquationInPoint(vector<vector<double>> coefs, vector<vector<double>> pows, vector<double> points);
	double* Gaus(double **brr, int row, int col, double*det, bool& mark, double epsilon);
	double** InverseMatrix(vector<vector<double>>arr, int size, double*det, bool& mark, double epsilon);
	vector<double> Multiple(vector<vector<double>> matrix, vector<double> vect);
	vector<double> Difference(vector<double> a, vector<double> b);
	double Max(vector<double> m);

public:
	Newton();
	Result SolutionOfTheSystem(vector<vector<double>> coefs, vector<vector<double>> pows, vector<double> points, double epsilon);
};