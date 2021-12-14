#include <iostream>
#include <algorithm>
#include <vector>
#include <iomanip>

int N = 5;
int L = 1;
const double delta = 1e-5;

using namespace std;

void SumSq(vector<double>& x)
{
	double sumsq = 0.0;
	for (int i = 0; i < N; i++)
	{
		sumsq += x[i] * x[i];
	}
	cout << "Vector sum of squares : " << sumsq << endl;
}

void Matrix_Zeroing(vector<vector<double>>& A1, vector<vector<double>>& A2, vector<vector<double>>& A3)
{
	for (int i = 0; i < N; i++)
	{
		vector<double> a;
		for (int j = 0; j < N; j++)
		{
			a.push_back(0.0);
		}
		A1.push_back(a);
		A2.push_back(a);
		A3.push_back(a);
	}
}

void Matrix_Generation_Off_Diag(vector<vector<double>>& A1, vector<vector<double>>& A2, vector<vector<double>>& A3)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = max(0, i - L); j <= min(i + L, N - 1); j++)
		{
			if (i != j)
			{
				double num;
				num = (double)rand() * 2.0 / RAND_MAX - 1.0;
				A1[i][j] = num;
				A2[i][j] = num;
				A3[i][j] = num;
			}
		}
	}
}

void Matrix_Generation_Diag(vector<vector<double>>& A, double q)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i != j)
			{
				A[i][i] += abs(A[i][j]);
			}
		}
		A[i][i] *= q;
	}
}

void Vector_Generation(vector<double>& x)
{
	double num;
	for (int i = 0; i < N; i++)
	{
		num = (double)rand() * 2.0 / RAND_MAX - 1.0;
		x.push_back(num);
	}
}

vector<double> Multiplication_MV(vector<vector<double>>& A, vector<double>& x)
{
	double sum;
	vector<double> b;
	for (int i = 0; i < N; i++)
	{
		sum = 0.0;
		for (int j = 0; j < N; j++)
		{
			sum += A[i][j] * x[j];
		}
		b.push_back(sum);
	}
	return b;
}

vector<double> Multiplication_MTV(vector<vector<double>>& A, vector<double>& x)
{
	double sum;
	vector<double> b;
	for (int i = 0; i < N; i++)
	{
		sum = 0.0;
		for (int j = 0; j < N; j++)
		{
			sum += A[j][i] * x[j];
		}
		b.push_back(sum);
	}
	return b;
}

vector<vector<double>> Multiplication_MTM(vector<vector<double>>& A)
{
	double sum;
	vector<vector<double>> Az;
	for (int i = 0; i < N; i++)
	{
		vector<double> b;
		for (int j = 0; j < N; j++)
		{
			sum = 0.0;
			for (int k = 0; k < N; k++)
			{
				sum += A[k][i] * A[k][j];
			}
			b.push_back(sum);
		}
		Az.push_back(b);
	}
	return Az;
}

void Matrix_Out(vector<vector<double>>& A, string name)
{
	cout << "Martix " << name << " : " << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << fixed << setprecision(3) << A[i][j] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

void Vector_Out(vector<double>& x, string name)
{
	cout << "Vector " << name << " : " << endl;
	for (int i = 0; i < N; i++)
	{
		cout << fixed << setprecision(3) << x[i] << "\t";
	}
	cout << endl << endl;
}

void Jacobi_Method(vector<vector<double>>& A, vector<double>& x, vector<double>& b, double q)
{
	int itr_count = 0;
	double max_diff;
	vector<double> x0, x1;
	for (int i = 0; i < N; i++)
	{
		x0.push_back(x[i]);
		x1.push_back(0.0);
	}
	do
	{
		itr_count++;
		max_diff = 0.0;
		for (int i = 0; i < N; i++)
		{
			x1[i] = b[i];
			for (int j = 0; j < N; j++)
			{
				if (i != j)
				{
					x1[i] -= A[i][j] * x0[j];
				}
			}
			x1[i] /= A[i][i];
			if (abs(x1[i] - x0[i]) > max_diff)
			{
				max_diff = abs(x1[i] - x0[i]);
			}
			x0[i] = x1[i];
		}
	} while (max_diff > delta);
	cout << "Jacobi method: " << endl;
	cout << "q: " << q << endl;
	cout << "Number of iterations: " << itr_count << endl;
	SumSq(x0);
	cout << "-----------------------------------------------" << endl;
}

void SOR_Method(vector<vector<double>>& A, vector<double>& x, vector<double>& b, double q, double w)
{
	int itr_count = 0;
	double max_diff;
	vector<double> x0, x1;
	for (int i = 0; i < N; i++)
	{
		x0.push_back(x[i]);
		x1.push_back(0.0);
	}
	do
	{
		itr_count++;
		max_diff = 0.0;
		for (int i = 0; i < N; i++)
		{
			x1[i] = b[i];
			for (int j = 0; j < i; j++)
			{
				x1[i] -= A[i][j] * x1[j];
			}
			for (int j = i + 1; j < N; j++)
			{
				x1[i] -= A[i][j] * x0[j];
			}
			x1[i] /= A[i][i];
			x1[i] = w * x1[i] + (1 - w) * x0[i];
			if (abs(x1[i] - x0[i]) > max_diff)
			{
				max_diff = abs(x1[i] - x0[i]);
			}
		}
		for (int i = 0; i < N; i++)
		{
			x0[i] = x1[i];
		}
	} while (max_diff > delta);
	cout << "Method SOR: " << endl;
	cout << "q: " << q << endl;
	cout << "w: " << w << endl;
	cout << "Number of iterations: " << itr_count << endl;
	SumSq(x0);
	cout << "-----------------------------------------------" << endl;
}

double Scalar_Multiplication(vector<double>& a, vector<double>& b)
{
	double sum = 0.0;
	for (int i = 0; i < N; i++)
	{
		sum += a[i] * b[i];
	}
	return sum;
}

void CGM_Method(vector<vector<double>>& A, vector<double>& x, vector<double>& b, double q)
{
	int itr_count = 0;
	double alpha, beta, max_diff;
	vector<double> r0, r1, p, x0, x1;
	for (int i = 0; i < N; i++)
	{
		x0.push_back(x[i]);
		x1.push_back(0.0);
		r0.push_back(b[i]);
		r1.push_back(0.0);
		for (int j = 0; j < N; j++)
		{
			r0[i] -= A[i][j] * x[j];
		}
		p.push_back(r0[i]);
	}
	do
	{
		itr_count++;
		max_diff = 0.0;
		vector<double> Ap = Multiplication_MV(A, p);
		alpha = Scalar_Multiplication(r0, r0) / Scalar_Multiplication(Ap, p);
		for (int i = 0; i < N; i++)
		{
			x1[i] = x0[i] + alpha * p[i];
			if (abs(x1[i] - x0[i]) > max_diff)
			{
				max_diff = abs(x1[i] - x0[i]);
			}
			r1[i] = r0[i] - alpha * Ap[i];
		}
		beta = Scalar_Multiplication(r1, r1) / Scalar_Multiplication(r0, r0);
		for (int i = 0; i < N; i++)
		{
			p[i] = r1[i] + beta * p[i];
			x0[i] = x1[i];
			r0[i] = r1[i];
		}
	} while (max_diff > delta);
	cout << "Method CGM: " << endl;
	cout << "q: " << q << endl;
	cout << "Number of iterations: " << itr_count << endl;
	SumSq(x0);
	cout << "-----------------------------------------------" << endl;
}

void Jacobi_Method_With_M_Steps(vector<vector<double>>& A, vector<double>& x0, vector<double>& b, double m)
{
	double max_diff;
	vector<double> x1;
	for (int i = 0; i < N; i++)
	{
		x1.push_back(0.0);
	}
	for (int k = 0; k < m; k++)
	{
		max_diff = 0.0;
		for (int i = 0; i < N; i++)
		{
			x1[i] = b[i];
			for (int j = 0; j < N; j++)
			{
				if (i != j)
				{
					x1[i] -= A[i][j] * x0[j];
				}
			}
			x1[i] /= A[i][i];
			if (abs(x1[i] - x0[i]) > max_diff)
			{
				max_diff = abs(x1[i] - x0[i]);
			}
			x0[i] = x1[i];
		}
		if (max_diff <= delta) break;
	}
}

void PCGM_Method(vector<vector<double>>& A, vector<double>& x, vector<double>& b, double q, double m)
{
	int itr_count = 0;
	double alpha, beta, max_diff;
	vector<double> r0, r1, rz0, rz1, p, x0, x1;
	for (int i = 0; i < N; i++)
	{
		x0.push_back(x[i]);
		x1.push_back(0.0);
		r0.push_back(b[i]);
		r1.push_back(0.0);
		rz0.push_back(0.0);
		rz1.push_back(0.0);
		for (int j = 0; j < N; j++)
		{
			r0[i] -= A[i][j] * x[j];
		}
	}
	Jacobi_Method_With_M_Steps(A, rz0, r0, m);
	for (int i = 0; i < N; i++)
	{
		p.push_back(rz0[i]);
	}
	do
	{
		itr_count++;
		max_diff = 0.0;
		vector<double> Ap = Multiplication_MV(A, p);
		alpha = -Scalar_Multiplication(rz0, r0) / Scalar_Multiplication(p, Ap);
		for (int i = 0; i < N; i++)
		{
			x1[i] = x0[i] - alpha * p[i];
			if (abs(x1[i] - x0[i]) > max_diff)
			{
				max_diff = abs(x1[i] - x0[i]);
			}
			r1[i] = r0[i] + alpha * Ap[i];
		}
		Jacobi_Method_With_M_Steps(A, rz1, r1, m);
		beta = Scalar_Multiplication(rz1, r1) / Scalar_Multiplication(rz0, r0);
		for (int i = 0; i < N; i++)
		{
			p[i] = rz1[i] + beta * p[i];
			x0[i] = x1[i];
			r0[i] = r1[i];
			rz0[i] = rz1[i];
			rz1[i] = 0.0;
		}
	} while (max_diff > delta);
	cout << "Method PCGM: " << endl;
	cout << "q: " << q << endl;
	cout << "m: " << m << endl;
	cout << "Number of iterations: " << itr_count << endl;
	SumSq(x0);
	cout << "-----------------------------------------------" << endl;
}

inline void generation(vector<vector<double>>& A1, vector<vector<double>>& A2, vector<vector<double>>& A3, vector<double>& x, vector<double>& y, vector<double>& b1, vector<double>& b2, vector<double>& b3) {
	cout << "Enter N and L" << endl;
	cin >> N >> L;
	cout << "generation begins" << endl;

	Matrix_Zeroing(A1, A2, A3);
	Matrix_Generation_Off_Diag(A1, A2, A3);
	Matrix_Generation_Diag(A1, 1.1);
	Matrix_Generation_Diag(A2, 2.0);
	Matrix_Generation_Diag(A3, 10.0);
	Vector_Generation(x);
	cout << "Exact answer: " << endl;
	SumSq(x);
	b1 = Multiplication_MV(A1, x);
	b2 = Multiplication_MV(A2, x);
	b3 = Multiplication_MV(A3, x);
	b1 = Multiplication_MTV(A1, b1);
	b2 = Multiplication_MTV(A2, b2);
	b3 = Multiplication_MTV(A3, b3);
	A1 = Multiplication_MTM(A1);
	A2 = Multiplication_MTM(A2);
	A3 = Multiplication_MTM(A3);
	Vector_Generation(y);
}

void Matrix_Zeroing_Smal (vector<vector<double>>& A1)
{
	for (int i = 0; i < N; i++)
	{
		vector<double> a;
		for (int j = 0; j < N; j++)
		{
			a.push_back(0.0);
		}
		A1.push_back(a);
	}
}
void Matrix_Generation_Off_Diag_Small(vector<vector<double>>& A1)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = max(0, i - L); j <= min(i + L, N - 1); j++)
		{
			if (i != j)
			{
				double num;
				num = (double)rand() * 2.0 / RAND_MAX - 1.0;
				A1[i][j] = num;
			}
		}
	}
}
void smallGeneration(vector<vector<double>>& A1, vector<double>& y, vector<double>& b1) {
	cout << "Enter N and L" << endl;
	cin >> N >> L;
	cout << "generation begins" << endl;

	Matrix_Zeroing_Smal(A1);
	Matrix_Generation_Off_Diag_Small(A1);
	Matrix_Generation_Diag(A1, 1.1);
	vector<double> x;
	Vector_Generation(x);
	cout << "Exact answer: " << endl;
	SumSq(x);
	b1 = Multiplication_MV(A1, x);
	b1 = Multiplication_MTV(A1, b1);
	A1 = Multiplication_MTM(A1);
	Vector_Generation(y);
}

void task_1() {
	vector<vector<double>> A1, A2, A3;
	Matrix_Zeroing(A1, A2, A3);
	Matrix_Generation_Off_Diag(A1, A2, A3);
	Matrix_Out(A1, "Before diag A1");
	Matrix_Generation_Diag(A1, 1.1);
	Matrix_Out(A1, "A1");
	Matrix_Generation_Diag(A2, 2.0);
	Matrix_Out(A2, "A2");
	Matrix_Generation_Diag(A3, 10.0);
	Matrix_Out(A3, "A3");
	vector<double> x, b1, b2, b3;
	Vector_Generation(x);
	Vector_Out(x, "x");
	cout << "Exact answer: " << endl;
	SumSq(x);
	cout << "-----------------------------------------------" << endl;
	b1 = Multiplication_MV(A1, x);
	Vector_Out(b1, "b1");
	b2 = Multiplication_MV(A2, x);
	Vector_Out(b2, "b2");
	b3 = Multiplication_MV(A3, x);
	Vector_Out(b3, "b3");
	b1 = Multiplication_MTV(A1, b1);
	b2 = Multiplication_MTV(A2, b2);
	b3 = Multiplication_MTV(A3, b3);
	A1 = Multiplication_MTM(A1);
	A2 = Multiplication_MTM(A2);
	A3 = Multiplication_MTM(A3);
	Matrix_Out(A1, "Result A1");
	Matrix_Out(A2, "Result A2");
	Matrix_Out(A3, "Result A3");
	Vector_Out(b1, "Result b1");
	Vector_Out(b2, "Result b2");
	Vector_Out(b3, "Result b3");
}

void task_2() {
	vector<vector<double>> A1, A2, A3;
	vector<double> x, y, b1, b2, b3;
	generation(A1, A2, A3, x, y, b1, b2, b3);
	Jacobi_Method(A1, y, b1, 1.1); 
	Jacobi_Method(A2, y, b2, 2.0);
	Jacobi_Method(A3, y, b3, 10.0);
}

void task_3() {
	vector<vector<double>> A1, A2, A3;
	vector<double> x, y, b1, b2, b3;
	generation(A1, A2, A3, x, y, b1, b2, b3);
	SOR_Method(A1, y, b1, 1.1, 0.25);
	SOR_Method(A1, y, b1, 1.1, 0.5);
	SOR_Method(A1, y, b1, 1.1, 1.0);
	SOR_Method(A1, y, b1, 1.1, 1.5);
	SOR_Method(A1, y, b1, 1.1, 1.75);
	SOR_Method(A2, y, b2, 2.0, 0.25);
	SOR_Method(A2, y, b2, 2.0, 0.5);
	SOR_Method(A2, y, b2, 2.0, 1.0);
	SOR_Method(A2, y, b2, 2.0, 1.5);
	SOR_Method(A2, y, b2, 2.0, 1.75);
	SOR_Method(A3, y, b3, 10.0, 0.25);
	SOR_Method(A3, y, b3, 10.0, 0.5);
	SOR_Method(A3, y, b3, 10.0, 1.0);
	SOR_Method(A3, y, b3, 10.0, 1.5);
	SOR_Method(A3, y, b3, 10.0, 1.75);
}

void task_4() {
	vector<vector<double>> A1, A2, A3;
	vector<double> x, y, b1, b2, b3;
	generation(A1, A2, A3, x, y, b1, b2, b3);
	CGM_Method(A1, y, b1, 1.1);
	CGM_Method(A2, y, b2, 2.0);
	CGM_Method(A3, y, b3, 10.0);
	cout << endl << "CGM Method done" << endl << "begin PCGM method" << endl;
	smallGeneration(A1, y, b1);
	PCGM_Method(A1, y, b1, 1.1, 1.0);
	PCGM_Method(A1, y, b1, 1.1, 2.0);
	PCGM_Method(A1, y, b1, 1.1, 3.0);
	PCGM_Method(A1, y, b1, 1.1, 4.0);
	PCGM_Method(A1, y, b1, 1.1, 5.0);
	PCGM_Method(A1, y, b1, 1.1, 6.0);
	PCGM_Method(A1, y, b1, 1.1, 7.0);
	PCGM_Method(A1, y, b1, 1.1, 8.0);
	PCGM_Method(A1, y, b1, 1.1, 9.0);
	PCGM_Method(A1, y, b1, 1.1, 10.0);
}

int main()
{
	setlocale(LC_ALL, "ru");
	srand(time(NULL));

	while (1) {
		auto taskNumber = -1;
		cout << "enter number for task(1,2,3,4), 0 for exit" << endl;
		cin >> taskNumber;
		switch (taskNumber)
		{
		case 1:
			task_1();
			break;
		case 2:
			task_2();
			break;
		case 3:
			task_3();
			break;
		case 4:
			task_4();


			break;
		case 0:
			return 0;
		}

	}

}
