//
// Created by hsingh9 on 18/03/2019.
//

#ifndef COURSEWORK1_ROOTFINDING_HPP
#define COURSEWORK1_ROOTFINDING_HPP

#include <iostream>
#include <cmath>
#include <Dense>
#include <vector>

using namespace Eigen;
using namespace std;

class RootFinding
{
public:
    VectorXd calculateF(double x, double y, double z);
    MatrixXd calculateJ(double x, double y, double z);
    bool foundNewRoot(double x, double y, double z, ArrayXXd root);
    VectorXd calculateF2C (double x, double y,double x0, double y0,double utSqNorm);
    MatrixXd calculateJ2C (double x, double y,double x0, double y0);
    bool foundNewRoot(double x, double y, ArrayXXd root);
    void writeToFile (const string fNamePrefix, ArrayXXd roots, const int k, vector<string> colNames);



};


#endif //COURSEWORK1_ROOTFINDING_HPP
