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


};


#endif //CODE_STOKESPDE_HPP
