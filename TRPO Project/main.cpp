#include "Newton_method.h"
#include "Runge-Kutta_method.h"


void Print(vector<vector<double>> m)
{
	for (size_t i = 0; i < m.size(); i++)
	{
		for (size_t j = 0; j < m[i].size(); j++)
		{
			cout << m[i][j] << " ";
		}
		cout << endl;
	}
}

int main()
{
	/*vector<vector<double>> coefficients = { { 1,2,0 },{ 0,2,4 } };
	vector<vector<double>> powers = { { 3,4 },{ 3,2 } };*/

	vector<vector<double>> coefficients = { { 1,1,1,-1 },{ 2,1,-4,0 },{3,-4,1,0} };
	vector<vector<double>> powers = { { 2,2,2 },{ 2,2,1 }, {2,1,2} };

	EquationSystem systemOfEquations(coefficients, powers);


	//Output of the system of equations
	cout << "The system of equations:" << endl;
	systemOfEquations.printSystem();
	cout << endl;


	//Jacobi matrix at point X0
	vector<double> point = { 0.5,0.5,0.5 };
	cout << "Jacobi matrix:" << endl;
	vector<vector<double>> matrixJacobi=systemOfEquations.matrixJacobi(point);
	Print(matrixJacobi);

	Newton *n = new Newton();
	n->SolutionOfTheSystem(coefficients, powers, point, 0.005);

	system("pause");
	return 0;
}