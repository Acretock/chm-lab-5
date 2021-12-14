// NM_5_test_test.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <fstream>
const double eps = 1e-4;
const int N = 3, l = 1;

template <class T> const T& max(const T& a, const T& b);
template <class T> const T& min(const T& a, const T& b);
#pragma region Matrix
class Matrix
{
public:
	Matrix();
	Matrix(int rows, int cols, double init);
	Matrix(int N, double init);

	Matrix& operator=(const Matrix& m);
	Matrix(std::vector<std::vector<double>> A);
	~Matrix();
	Matrix transpose();
	double dotProduct(Matrix& a, Matrix& b);
	std::vector<std::vector<double>> m;
	int cols_;
	int rows_;
private:
};

Matrix::Matrix()
{
}

Matrix::Matrix(int rows, int cols, double init) {
	m.resize(rows);
	for (int i = 0; i < rows; i++) {
		m[i].resize(cols, init);
	}
	rows_ = m.size();
	cols_ = m[0].size();
}

Matrix::Matrix(int N, double init)
{
	m.resize(N);
	for (int i = 0; i < N; i++) {
		m[i].resize(N, init);
	}
	rows_ = m.size();
	cols_ = m[0].size();
}

Matrix& Matrix::operator=(const Matrix& m)
{
	this->m = m.m;
	this->cols_ = m.cols_;
	this->rows_ = m.rows_;
	return *this;
}

Matrix::Matrix(std::vector<std::vector<double>> A)
{
	m = A;
	rows_ = m.size();
	cols_ = m[0].size();
}
Matrix Matrix::transpose()
{
	Matrix ret(rows_, cols_, 0);
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			ret.m[j][i] = m[i][j];
		}
	}
	return ret;
}
double Matrix::dotProduct(Matrix& a, Matrix& b)
{
	double sum = 0;
	for (int i = 0; i < a.rows_; ++i) {
		sum += (a.m[i][0] * b.m[i][0]);
	}
	return sum;
}
Matrix::~Matrix()
{
}
Matrix operator+(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, const Matrix& b);
Matrix operator-(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, double b);
Matrix operator*(double b, const Matrix& a);
Matrix operator/(const Matrix& a, double b);

#pragma endregion
#pragma region Основые методы
//Генерирует случайное число в заданном диапазоне
double random_in_range(double a, double b);
//Генерирует матрицу по заданным условиям
void generateMatrix(std::vector<Matrix>& A, double* q);
//Генерирует вектор по заданным условиям
void generateVector(Matrix& x);
std::vector<Matrix> CGM(std::vector<Matrix>& A, std::vector<Matrix>& b, Matrix& prevX, double* iterations);
//Функция произведения ленточной матрицы на обычную
Matrix mul(Matrix& a, Matrix& b);
void matrixOut(Matrix& A);
#pragma endregion

std::vector <Matrix> SOR(std::vector <Matrix>& A, std::vector <Matrix>& b, Matrix& x, double* iterations, double w, std::ofstream& out1, std::ofstream& out2, std::ofstream& out3)
{
	//std::vector <Matrix> answer(3);
	std::vector <Matrix> answer(3);

	for (int i = 0; i < 3; i++) {

		generateVector(x);

		double delta = 10000;
		double iterCount = 0;
		Matrix result(N, 1);
		Matrix tempA = A[i];
		Matrix tempB = b[i];
		int iteration = 0;
		while (delta > eps) {
			delta = 0;
			for (int j = 0; j < N; j++) {
				double sum = (1 - w) * x.m[j][0] + (w / tempA.m[j][j]) * tempB.m[j][0];
				double temp = 0;
				for (int k = 0; k < j; k++) {
					temp += tempA.m[j][k] * result.m[k][0];
				}
				sum -= temp * (w / tempA.m[j][j]);
				temp = 0;
				for (int k = j + 1; k < N; k++) {
					temp += tempA.m[j][k] * x.m[k][0];
				}
				sum -= temp * (w / tempA.m[j][j]);
				result.m[j][0] = sum;



			}
			for (int k = 0; k < N; k++) {
				double m = result.m[k][0];
				double n = x.m[k][0];
				double error = abs(m - n);
				delta = max(error, delta);
			}
			iterations[i] += 1;
			iteration++;
			x = result;
		}
		if (i == 0) out1 << "Для w = " << w << "  количество итераций равно " << iteration << std::endl;
		if (i == 1) out2 << "Для w = " << w << "  количество итераций равно " << iteration << std::endl;
		if (i == 2) out3 << "Для w = " << w << "  количество итераций равно " << iteration << std::endl;
		answer[i] = result;
	}
	return answer;

	/*generateVector(x);
	for (int i = 2; i < 3; i++) {
		double delta = 10000;
		double iterCount = 0;
		Matrix result(N, 1,0);
		iterations[i] = 0;
		Matrix tempA = A[i];
		Matrix tempB = b[i];
		int iteration = 0;
		while (sqrt(delta) > eps) {
			delta = 0;
			for (int j = 0; j < N; j++) {
				double sum = (1 - w) * x.m[j][0] + (w / tempA.m[j][j]) * tempB.m[j][0];
				double temp = 0;
				for (int k = 0; k < j; k++) {
					temp += tempA.m[j][k] * result.m[k][0];
				}
				sum -= temp * (w / tempA.m[j][j]);
				temp = 0;
				for (int k = j +1; k < N; k++) {
					temp += tempA.m[j][k] * x.m[k][0];
				}
				sum -= temp * (w / tempA.m[j][j]);
				result.m[j][0] = sum;

				iterations[i] ++;
				iteration += 1;

			}

			Matrix error = result - x;
			for (int i = 0; i < N; i++) {
				delta += error.m[i][0] * error.m[i][0];
			}
			x = result;
		}
		std::cout << "Для w = " << w << " количество итераций равно " << iteration << std::endl;;
		answer[i] = result;
	}*/
	//return answer;

}

void task_3(std::vector <Matrix>& A, std::vector <Matrix>& b, Matrix& x) {
	std::ofstream out1("tsk_3_out_1.1.txt");
	std::ofstream out2("tsk_3_out_2.txt");
	std::ofstream out3("tsk_3_out_10.txt");
	std::vector<Matrix> SORResults(N);
	double iterations[3] = { 0,0,0 };
	generateVector(x);
	SORResults = SOR(A, b, x, iterations, 1, out1, out2, out3);
	/*for (double w = 0.1; w < 2; w += 0.1) {
		Matrix temp = Matrix(x);

		SORResults = SOR(A, b, x, iterations, 1);

	}*/
	//std::cout << "Коэффициент w = " << 1 << std::endl;
	//std::cout << "Решение системы 1 при q = 1.1: " << std::endl;// matrixOut(SORResults[0]);
	//std::cout << "Необходимое число итераций: " << iterations[0] << std::endl << std::endl << "Решение системы 2 при q = 2: " << std::endl;// matrixOut(SORResults[1]);
	//std::cout << "Необходимое число итераций: " << iterations[1] << std::endl << std::endl << "Решение системы 3 при q = 10: " << std::endl; //matrixOut(SORResults[2]);
	//std::cout << "Необходимое число итераций: " << iterations[2] << std::endl << std::endl;
	/*for (double w = 0.1; w < 2; w += 0.1) {
		Matrix temp = Matrix(x);

		SORResults = SOR(A, b, x, iterations, w, out1, out2, out3);

	}*/
	std::cout << "Коэффициент w = " << 1 << std::endl;
	std::cout << "Решение системы 1: " << std::endl; matrixOut(SORResults[0]);
	std::cout << "Необходимое число итераций: " << iterations[0] << std::endl << std::endl << "Решение системы 2: " << std::endl; matrixOut(SORResults[1]);
	std::cout << "Необходимое число итераций: " << iterations[1] << std::endl << std::endl << "Решение системы 3: " << std::endl; matrixOut(SORResults[2]);
	std::cout << "Необходимое число итераций: " << iterations[2] << std::endl << std::endl;
	out1.close(); out2.close(); out3.close();
}
#pragma endregion
#pragma region TASK4
Matrix PCJacobiMethod(Matrix A, Matrix b, int steps, Matrix& x, int& iter)
{

	// разность значений соседних итераций.
	// разность значений соседних итераций.
	double delta = 10000;
	Matrix result(N, 1, 0);
	Matrix tempA = A;
	Matrix tempB = b;
	srand(time(NULL));
	generateVector(x);
	for (int z = 0; z < steps; z++) {
		delta = 0;
		for (int j = 0; j < N; j++) {
			double sum = tempB.m[j][0];
			for (int g = 0; g < N; g++) {


				if (j != g) {
					sum -= tempA.m[j][g] * x.m[g][0];
				}
			}
			result.m[j][0] = sum / tempA.m[j][j];



		}
		iter = z;
		double sum = 0;
		Matrix error = result - x;
		for (int i = 0; i < N; i++) {
			sum += error.m[i][0] * error.m[i][0];
		}
		if (sqrt(sum) < eps)
			break;
		x = result;
	}
	return result;

}


std::vector <Matrix> CGM(std::vector <Matrix>& A, std::vector <Matrix>& b, Matrix& prevX, double* iterations)
{
	std::vector <Matrix> answer(3);


	for (int i = 0; i < 3; i++) {
		double delta = 10000;
		double iterCount = 0;
		Matrix tempA = A[i];
		Matrix tempB = b[i];
		generateVector(prevX);
		Matrix prevR(N, 1, 1), prevP(N, 1, 1);
		Matrix R(N, 1, 1), P(N, 1, 1), result(N, 1, 1);
		prevR = tempB - tempA * prevX;
		prevP = prevR;
		double product = prevP.dotProduct(prevR, prevR);




		while (delta > eps) {
			Matrix temp = tempA * prevP;
			double alpha = -prevP.dotProduct(prevR, prevR) / prevR.dotProduct(prevP, temp);
			result = prevX - alpha * prevP;
			R = prevR + alpha * temp;
			double product2 = prevP.dotProduct(R, R);
			if (sqrt(product2) >= eps) {
				double betha = product2 / prevP.dotProduct(prevR, prevR);
				P = R + betha * prevP;
				prevP = P;
				prevR = R;
				prevX = result;
				iterations[i]++;
			}
			else break;

		}
		answer[i] = result;
	}
	return answer;

}

std::vector <Matrix> PCGM(std::vector <Matrix>& A, std::vector <Matrix>& b, Matrix& prevX, double* iterations, std::ofstream& out)
{
	std::vector <Matrix> answer(3);

	Matrix m = prevX;

	for (int i = 0; i < 1; i++) {
		for (int iterCount = 50; iterCount > 5; iterCount--) {

			Matrix tempx(N, 1, 0);
			prevX = m;
			double delta = 10000;
			Matrix tempA = A[i];
			Matrix tempB = b[i];
			Matrix prevR(N, 1, 0), prevRW(N, 1, 0), prevP(N, 1, 0);
			Matrix R(N, 1, 0), RW(N, 1, 0), P(N, 1, 0), result(N, 1, 0);
			int iter = 0;
			Matrix temp = A[i] * prevX;
			prevR = tempB - temp;

			prevRW = PCJacobiMethod(A[i], prevR, iterCount, tempx, iter);

			prevP = prevR;

			double product = prevP.dotProduct(prevRW, prevR);

			while (delta > eps) {

				iterations[i] ++;
				double temp1 = prevP.dotProduct(prevRW, prevR);
				Matrix temp2 = A[i] * prevP;
				double temp3 = prevP.dotProduct(prevP, temp2);
				double alpha = -temp1 / temp3;
				temp = alpha * prevP;
				result = prevX - temp;
				R = prevR + alpha * temp2;
				product = prevP.dotProduct(R, R);
				delta = sqrt(product);
				generateVector(tempx);
				RW = PCJacobiMethod(A[i], R, iterCount, tempx, iter);

				double betha = prevP.dotProduct(RW, R) / (prevP.dotProduct(prevRW, prevR));

				P = RW + betha * prevP;

				prevX = result;
				prevP = P;
				prevRW = RW;
				prevR = R;

			}
			out << iterations[i] << std::endl;
			std::cout << "solved: iter = " << iterCount << "\tРезультат \t " << result.m[0][0] << "  " << result.m[1][0] << std::endl;
			iterations[i] = 0;
		}

		std::cout << "Количество итераций метода Якоби для задачи " << i + 1 << " : " << iter << std::endl;
		answer[i] = result;
	}

	return answer;
}



void task_4(std::vector <Matrix>& A, std::vector <Matrix>& b, Matrix& x) {
	std::ofstream out("task4.txt");
	std::vector<Matrix> CGMResults(N);
	std::vector<Matrix> PCGMResults(N);
	double iterations[3] = { 0,0,0 };

	CGMResults = CGM(A, b, x, iterations);
	std::cout << "     Method CGM    " << std::endl;
	std::cout << "Решение системы 1: " << std::endl; matrixOut(CGMResults[0]);
	std::cout << "Необходимое число итераций: " << iterations[0] << std::endl << std::endl << "Решение системы 2: " << std::endl; matrixOut(CGMResults[0]);
	std::cout << "необходимое число итераций: " << iterations[1] << std::endl << std::endl << "решение системы 3: " << std::endl; matrixOut(CGMResults[0]);
	std::cout << "Необходимое число итераций: " << iterations[2] << std::endl << std::endl;
	out.close();

	/*PCGMResults = PCGM(A, b, x, iterations, out);
	std::cout << "/*******************" << std::endl;
	std::cout << "*    Метод PCGM    *" << std::endl;
	std::cout << "********************//*" << std::endl;
	matrixOut(PCGMResults[0]);
	std::cout << "Решение системы 1: " << std::endl; matrixOut(PCGMResults[0]);
	std::cout << "Необходимое число итераций: " << iterations[0] << std::endl << std::endl << "Решение системы 2: " << std::endl; matrixOut(PCGMResults[0]);
	std::cout << "необходимое число итераций: " << iterations[1] << std::endl << std::endl << "решение системы 3: " << std::endl; matrixOut(PCGMResults[0]);
	std::cout << "Необходимое число итераций: " << iterations[2] << std::endl << std::endl;
	out.close();*/
}
#pragma endregion
int main()
{
	//srand(time(NULL));
	setlocale(LC_ALL, "Russian");
	/*******************
	* Начальные условия*
	*	  Задание 1    *
	********************/
	//Исходные матрицы
	std::vector<Matrix> A(3);
	std::vector<Matrix> At(3);
	std::vector<Matrix> Astar(3);
	std::vector<Matrix> B(3);
	std::vector<Matrix> Bstar(3);
	Matrix x = Matrix(N, 1, 0.0);

	double q[3] = { 1.1, 2, 10 };

	generateMatrix(A, q);

	std::cout << std::endl;
	generateVector(x);


	for (int i = 0; i < 3; i++) {


		B[i] = A[i] * x;


		At[i] = A[i].transpose();
		Astar[i] = At[i] * A[i];
		Bstar[i] = At[i] * B[i];


	}
	/*for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			std::cout << std::fixed << std::setprecision(2) <<  Astar[2].m[i][j] << "\t";
		}
		std::cout << std::endl;
	}*/
	//task_4(Astar, Bstar, x);
	while (1) {
		task_3(Astar, Bstar, x);
		int i;
		std::cin >> i;
		if (i == 0) {
			break;
		}
	}
}






Matrix mul(Matrix& a, Matrix& b)
{
	Matrix temp(a.rows_, b.cols_, 0);
	for (int i = 0; i < temp.rows_; ++i) {
		for (int j = 0; j < temp.cols_; ++j) {
			for (int k = 0; k < a.cols_; ++k) {
				temp.m[i][j] += (a.m[i][k] * b.m[k][j]);
			}
		}
	}
	return temp;
}

void matrixOut(Matrix& A) {
	for (auto a : A.m) {
		std::cout << std::fixed << std::setprecision(4) << a[0] << "\t";

		std::cout << std::endl;
	}
}

Matrix operator+(const Matrix& a, const Matrix& b)
{
	Matrix temp(a.rows_, b.cols_, 0);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] + b.m[i][j];
		}
	}
	return temp;
}

Matrix operator*(const Matrix& a, const Matrix& b)
{
	Matrix temp(a.rows_, b.cols_, 0);
	for (int i = 0; i < a.rows_; i++)
	{
		for (int j = 0; j < b.cols_; j++)
		{
			temp.m[i][j] = 0;
			for (int k = 0; k < b.rows_; k++)
			{
				temp.m[i][j] += a.m[i][k] * b.m[k][j];
			}
		}
	}
	return temp;
}

Matrix operator-(const Matrix& a, const Matrix& b)
{
	Matrix temp(a.rows_, b.cols_, 0);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] - b.m[i][j];
		}
	}
	return temp;
}

Matrix operator*(const Matrix& a, double b)
{
	Matrix temp(a.rows_, a.cols_, 0);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] * b;
		}
	}
	return temp;
}

Matrix operator*(double b, const Matrix& a)
{
	Matrix temp(a.rows_, a.cols_, 0);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] * b;
		}
	}
	return temp;
}

Matrix operator/(const Matrix& a, double b)
{
	Matrix temp(a.rows_, a.cols_);
	for (int i = 0; i < a.rows_; i++) {
		for (int j = 0; j < a.cols_; j++) {
			temp.m[i][j] = a.m[i][j] / b;
		}
	}
	return temp;
}

double random_in_range(double a, double b)
{
	double f = (double)rand() / RAND_MAX;
	return a + f * (b - a);
}

void generateMatrix(std::vector<Matrix>& A, double* q) {
	Matrix TEMP(N, 0);
	{


		for (signed i = 0; i < signed(N); ++i)
			for (signed j = 0; j < signed(N); j++) {
				if (j >= max(0, i - l) && j <= min(i + l, N)) {
					TEMP.m[i][j] = random_in_range(-1, 1);
				}
				else {
					TEMP.m[i][j] = 0;
				}
			}
		for (int i = 0; i < 3; i++) {
			for (signed j = 0; j < signed(N); j++) {
				double sum = 0;

				for (signed k = 0; k < signed(N); k++) {
					if (k != j)
						sum += abs(TEMP.m[k][j]);
				}
				TEMP.m[j][j] = q[i] * sum;
			}
			;
			A[i] = TEMP;
		}

	}




}

void generateVector(Matrix& x)
{
	for (int i = 0; i < N; i++) {
		x.m[i][0] = random_in_range(-1, 1);
	}
}


template <class T> const T& max(const T& a, const T& b) {
	return (a < b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
}
template <class T> const T& min(const T& a, const T& b) {
	return (a < b) ? a : b;     // or: return comp(a,b)?b:a; for version (2)
}