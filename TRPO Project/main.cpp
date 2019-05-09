#include "Newton_method.h"
#include "Runge-Kutta_method.h"


int main()
{
	vector<vector<double>> coefficients = { {1,2,0}, {0,2,4} };
	vector<vector<double>> powers = { {3,4}, {3,2} };

	EquationSystem systemOfEquations(coefficients, powers);


	//Output of the system of equations
	cout << "The system of equations:" << endl; 
	systemOfEquations.printSystem();
	cout << endl;
	

	//Jacobi matrix at point 2
	double point = 2;
	cout << "Jacobi matrix:" << endl;
	systemOfEquations.matrixJacobi(point);



	system("pause");
	return 0;
}