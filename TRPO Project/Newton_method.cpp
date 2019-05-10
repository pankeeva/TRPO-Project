#include "Newton_method.h"

void Newton::Print(vector<vector<double>> m)
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

Newton::Newton() {}

vector<double> Newton::EquationInPoint(vector<vector<double>> coefs, vector<vector<double>> pows, vector<double> points)
{
	vector<double> Fx0;
	for (int i = 0; i < coefs.size(); i++)
	{
		double result = 0;
		for (int j = 0; j < coefs[i].size() - 1; j++)
		{
			result += pow(points[i], pows[i][j]) * coefs[i][j];
		}
		Fx0.push_back(result + coefs[i][coefs[i].size() - 1]);
	}
	return Fx0;
}

Result Newton::SolutionOfTheSystem(vector<vector<double>> coefs, vector<vector<double>> pows, vector<double> points, double epsilon)
{
	Result result;
	EquationSystem e(coefs, pows);
	bool is_first = true;
	vector<double> multi;
	vector<double> r;
	while (true)
	{
		vector<vector<double>> jacobi = e.matrixJacobi(points);
		vector<double> Fx0 = this->EquationInPoint(coefs, pows, points);
		for (int i = 0; i < Fx0.size(); i++)
		{
			cout << "Fx0: " << Fx0[i] << " ";
		}
		cout << endl;

		bool mark = false;

		double*det = new double;
		if (is_first)
		{
			double**inverse = InverseMatrix(jacobi, jacobi.size(), det, mark, epsilon);
			if (*det == 0)
			{
				result.is_solution = false;
				return result;
			}
			vector < vector<double>> inv(jacobi.size(), vector<double>(jacobi.size()));
			for (int i = 0; i < jacobi.size(); i++)
			{
				for (int j = 0; j < jacobi.size(); j++)
				{
					inv[i][j] = inverse[i][j];
				}
			}
			Print(inv);

			cout << "R:" << endl;
			multi = Multiple(inv, Fx0);
			r = Difference(points, multi);
			for (int i = 0; i < r.size(); i++)
			{
				cout << r[i] << " ";
			}
			cout << endl;
		}
		else
		{
			double**br = new double *[jacobi.size()];
			for (int i = 0; i < jacobi.size(); i++)
			{
				br[i] = new double[jacobi.size() + 1];
				for (int j = 0; j < jacobi.size(); j++)
				{
					br[i][j] = jacobi[i][j];
				}
				br[i][jacobi.size()] = Fx0[i];
			}
			double* gaus = Gaus(br, jacobi.size(), jacobi.size() + 1, det, mark, epsilon);
			for (int i = 0; i < r.size(); i++)
			{
				r[i] += gaus[i];
				cout << "delta x: "<<gaus[i] << " ";
			}
			cout <<endl<< "R:" << endl;
			for (int i = 0; i < r.size(); i++)
			{
				cout << r[i] << " ";
			}
			cout << endl;
			vector<double> mm;
			for (int i = 0; i < r.size(); i++)
			{
				mm.push_back(gaus[i]);
			}
			multi = mm;
		}
		is_first = false;
		double max = Max(multi);
		if (max > epsilon)
		{
			//points = r;
			for (int i = 0; i < r.size(); i++)
			{
				points[i] = r[i];
			}
			continue;
		}
		result.matrix = r;
		result.is_solution = true;
		break;
	}

	return result;
}

vector<double> Newton::Multiple(vector<vector<double>> matrix, vector<double> vect)
{
	vector<double> res(vect.size());
	for (int i = 0; i < res.size(); i++)
	{
		res[i] = 0;
	}
	for (int j = 0; j < vect.size(); j++)
	{
		for (int i = 0; i < vect.size(); i++)
		{
			res[i] += matrix[i][j] * vect[j];
		}

	}
	return res;
}

vector<double> Newton::Difference(vector<double> a, vector<double> b)
{
	vector<double> res(a.size());
	for (int i = 0; i < res.size(); i++)
	{
		res[i] = 0;
	}
		for (int i = 0; i < a.size(); i++)
		{
			res[i] = a[i]-b[i];
		}
	return res;
}

double Newton::Max(vector<double> m)
{
	double max = m[0];
	for (int i = 1; i < m.size(); i++)
	{
		if (abs(m[i]) > max)
			max = abs(m[i]);
	}
	return max;
}

double Newton::Det(double**arr, int size, int p)
{
	double dob = 1;
	for (int i = 0; i < size; i++)
	{
		dob *= arr[i][i];
	}
	dob *= pow(-1, p);
	return dob;
}

double* Newton::Gaus(double **brr, int row, int col, double*det, bool& mark, double epsilon)
{
	double**arr = new double*[row];
	for (int i = 0; i < row; i++)
	{
		arr[i] = new double[col];
		for (int j = 0; j < col; j++)
		{
			arr[i][j] = brr[i][j];
		}
	}
	//Print(arr, row, col);
	int p = 0;
	for (int j = 0; j < col - 1; j++)
	{
		double max = arr[j][j];
		int indexrowmax = j;
		for (int i = j; i < row; i++)
		{
			if (abs(arr[i][j]) >= max)
			{
				max = abs(arr[i][j]);
				indexrowmax = i;
			}
		}
		if (max < epsilon)
			mark = true;
		//cout << "Determinant=0.Inverse matrix isn't exist";
		else
		{
			if (j != indexrowmax)
			{
				double tmp;
				for (int i = 0; i < col; i++)
				{
					tmp = arr[j][i];
					arr[j][i] = arr[indexrowmax][i];
					arr[indexrowmax][i] = tmp;
				}
				p++;
			}
			for (int i = j + 1; i < row; i++)
			{
				double m = -(arr[i][j] / arr[j][j]);
				double secondnumb = arr[j][j];
				//cout << m<<endl;
				for (int k = j; k < col; k++)
				{
					arr[i][k] = (arr[j][k] * m + arr[i][k])*secondnumb;
				}
			}

			//Print(arr, row, col);

		}

	}
	*det = Det(arr, row, p);
	double*results = new double[row];
	for (int i = 0; i < row; i++)
	{
		results[i] = 0;
	}
	for (int i = row - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = col - 2; j > i; j--)
		{
			sum += arr[i][j] * results[j];
		}
		results[i] = (arr[i][col - 1] - sum) / arr[i][i];
	}
	//cout << "..............." << endl;
	//for (int i = 0; i < row; i++)
	//{
	//	cout << results[i] << " ";

	//}
	cout << endl<< endl;
	return results;
}

double** Newton::InverseMatrix(vector<vector<double>>arr, int size, double*det, bool& mark, double epsilon)
{
	double**inverse = new double*[size];
	for (int i = 0; i < size; i++)
	{
		inverse[i] = new double[size];
	}
	double**matrix = new double*[size];
	for (int j = 0; j < size; j++)
	{
		matrix[j] = new double[size + 1];
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix[i][j] = arr[i][j];
		}
		matrix[i][size] = 0;
	}
	matrix[0][size] = 1;
	for (int i = 0; i < size; i++)
	{
		double*results = Gaus(matrix, size, size + 1, det, mark, epsilon);
		for (int j = 0; j < size; j++)
		{
			inverse[j][i] = results[j];
		}
		delete[]results;
		matrix[i][size] = 0;
		if (i + 1 != size)
			matrix[i + 1][size] = 1;
	}
	//Print(inverse, size, size);
	return inverse;
}