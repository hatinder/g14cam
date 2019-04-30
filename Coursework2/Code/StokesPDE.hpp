 //
// Created by hsingh9 on 24/04/2019.
//

#ifndef CODE_STOKESPDE_HPP
#define CODE_STOKESPDE_HPP
#include <Dense>
#include <Sparse>

using namespace Eigen;

class StokesPDE
{
public:
    SparseMatrix<double> createA (int N);
    VectorXd createB (double (*f) (double, double), double (*g) (double, double), int N, double a, double b);
    SparseMatrix<double> createBx(int N);
    SparseMatrix<double> createBy(int N);
    SparseMatrix<double> createZ (int N);
    SparseMatrix<double> createZN (int N);
    SparseMatrix<double>
    createC (SparseMatrix<double> A, SparseMatrix<double> Bx, SparseMatrix<double> By, SparseMatrix<double> Z,
             SparseMatrix<double> ZN);
    VectorXd createF(VectorXd Fu, VectorXd Fv, VectorXd Fp);

};


#endif //CODE_STOKESPDE_HPP
