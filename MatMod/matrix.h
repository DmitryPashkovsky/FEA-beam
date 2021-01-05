

#ifndef __matr_h_
#define __matr_h_

#include <iostream>
#include <math.h>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <algorithm>


template <class Value>
int Sign(Value Val)
{
    if (Val == 0.0)
        return 0;
    if (Val > 0.0)
        return 1;
    else
        return -1;
}

class matrix
{
public:
    double **M; 
    int m, n; 

    matrix(int _m, int _n); 
    matrix(void);
    matrix(std::string filename);
    matrix(const matrix& A);
    void rand_matrix(double k);
    void triang_matrix(void);
    void print(void);
    matrix& operator=(const matrix& A);
    matrix& operator*(double C);
    matrix trans(void);
    void file_write(std::string filename);
    double norm2_vec(matrix& x);
    matrix& operator+(matrix& A);
    matrix& operator-(matrix& A);
    matrix operator*(const matrix& A);
    ~matrix(void);
};
class SLAE
{
public:
    matrix A; // matrix of system
    matrix b; // vector of free variables
    matrix L; // left uni triangle matrix
    matrix D; // diagonal matrix
    matrix U; // Right uni triagle matrix
    matrix x; // unknown variables vector

    SLAE(std::string f_A, std::string f_b);
    SLAE(matrix& A_, matrix& b_);
    SLAE(int n);
    void LDU_decomposition(void);
    void reverse_substitution(void);
}; 

#endif // !__matr_h_

