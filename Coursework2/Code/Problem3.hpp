//
// Created by hsingh9 on 24/04/2019.
//

#ifndef CODE_PROBLEM3_HPP
#define CODE_PROBLEM3_HPP


#include <Dense>
#include <Sparse>

using  namespace Eigen;

class Problem3
{
public:
    VectorXd findU (SparseMatrix<double> SpA, VectorXd F);
    VectorXd findBigU (SparseMatrix<double> SpA, VectorXd F);
    void A();
    void B();
};


#endif //CODE_PROBLEM3_HPP
