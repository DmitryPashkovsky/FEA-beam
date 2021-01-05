#include "matrix.h"

matrix::matrix(int _m, int _n) : m(_m), n(_n)
{

    M = new double* [m];
    for (int i = 0; i < m; i++)
        M[i] = new double[n];

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            M[i][j] = 0;
}

matrix::matrix(void)
{
    m = 0;
    n = 0;
    M = nullptr;
}


matrix::matrix(std::string filename)
{
    std::ifstream file;


    file.open(filename);
    if (!file.is_open())
    {
        std::cout << "File is not opened!\n";
        return;
    }
    file >> m;
    file >> n;
    M = new double* [m];

    for (int i = 0; i < m; i++)
        M[i] = new double[n];

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            file >> M[i][j];

    file.close();
}
matrix::matrix(const matrix& A)
{
    m = A.m;
    n = A.n;
    if ((M = new double* [A.m]) == nullptr)
        return;

    for (int i = 0; i < A.m; i++)
        if ((M[i] = new double[A.n]) == nullptr)
            return;

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            M[i][j] = A.M[i][j];
}

void matrix::rand_matrix(double k)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            if (i == j)
                M[i][j] = 1;
            else
                M[i][j] = 1;
}

void matrix::triang_matrix(void)
{

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            if (i == j)
                M[i][j] = 1;
            else if (j > i)
                M[i][j] = 1;
            else
                M[i][j] = 0.0;
}

void matrix::print(void)
{
    std::cout << "m = " << m << " n = " << n << "\n";

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%10.5lf ", M[i][j]);
        std::cout << "\n";
    }
    std::cout << "\n\n";
}
matrix& matrix::operator=(const matrix& A)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            M[i][j] = A.M[i][j];
    return *this;
}
matrix& matrix::operator*(double C)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            M[i][j] *= C;
    return *this;
}


matrix matrix::trans(void)
{
    matrix temp(n, m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            temp.M[j][i] = M[i][j];
    return temp;
}



void matrix::file_write(std::string filename)
{
    std::ofstream file;


    file.open(filename);
    if (!file.is_open())
    {
        std::cout << "File is not opened!\n";
        return;
    }
    file << m;
    file << "\n";
    file << n;
    file << "\n";

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            file << M[i][j] << " ";
        file << "\n";
    }
    file.close();

}
double matrix::norm2_vec(matrix& x)
{
    double s = 0.0;

    if (x.m >= 1 && x.n == 1)
    {
        for (int i = 0; i < x.m; i++)
            s += x.M[i][0] * x.M[i][0];

        return sqrt(s);
    }
    else if (x.m == 1 && x.n >= 1)
    {
        for (int i = 0; i < x.m; i++)
            s += x.M[0][i] * x.M[0][i];
        return sqrt(s);
    }
    else
        return -1;
}



matrix& matrix::operator+(matrix& A)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            M[i][j] += A.M[i][j];
    return *this;
}

matrix& matrix::operator-(matrix& A)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            M[i][j] -= A.M[i][j];
    return *this;
}


matrix matrix::operator*(const matrix& A)
{
    matrix C(m, A.n);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < A.n; j++)
        {
            C.M[i][j] = 0;
            for (int k = 0; k < n; k++)
                C.M[i][j] += M[i][k] * A.M[k][j];
        }

    return C;
}


matrix::~matrix(void)
{
    for (int i = 0; i < m; i++)
        delete[] M[i];
    delete[] M;
}


SLAE::SLAE(std::string f_A, std::string f_b) : A(f_A), b(f_b), x(A.m, 1), L(A.m, A.n), D(A.m, A.n), U(A.m, A.n)
{

}

SLAE::SLAE(matrix& A_, matrix& b_) : A(A_.m, A_.n), b(b_.m, 1), x(A.m, 1), L(A.m, A.n), D(A.m, A.n), U(A.m, A.n)
{
    A = A_;
    b = b_;
}


SLAE::SLAE(int n) : A(n, n), b(n, 1), L(n, n), D(n, n), U(n, n), x(n, 1)
{
    A.rand_matrix(1.221);
    b.rand_matrix(1);
}

void SLAE::LDU_decomposition(void)
{
    for (int i = 0; i < A.n; i++)
        for (int j = 0; j < A.n; j++)
        {
            if (i == j)
            {
                double s = 0.0;

                for (int k = 0; k <= i - 1; k++)
                    s += L.M[i][k] * D.M[k][k] * U.M[k][i];
                D.M[i][i] = A.M[i][i] - s;
            }

            if (i > j)
            {
                double s = 0.0;
                for (int k = 0; k <= j - 1; k++)
                    s += L.M[i][k] * D.M[k][k] * U.M[k][j];

                L.M[i][j] = (A.M[i][j] - s) / D.M[j][j];
            }
            if (j > i)
            {
                double s = 0.0;
                for (int k = 0; k <= i - 1; k++)
                    s += L.M[i][k] * D.M[k][k] * U.M[k][j];

                U.M[i][j] = (A.M[i][j] - s) / D.M[i][i];
            }
        }
}



void SLAE::reverse_substitution(void)
{
    matrix y(A.n, 1);

    y.M[0][0] = b.M[0][0];

    for (int i = 1; i < A.n; i++)
    {
        double s = 0.0f;

        for (int j = 0; j < A.n; j++)
            s += L.M[i][j] * y.M[j][0];

        y.M[i][0] = b.M[i][0] - s;
    }


    for (int i = 0; i < A.n; i++)
        y.M[i][0] /= D.M[i][i];

    x.M[A.n - 1][0] = y.M[A.n - 1][0];

    for (int i = A.n - 2; i >= 0; i--)
    {
        double s = 0.0;

        for (int j = i + 1; j < A.n; j++)
            s += U.M[i][j] * x.M[j][0];
        x.M[i][0] = (y.M[i][0] - s);
    }
}
