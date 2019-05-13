//
// Created by hsingh9 on 17/04/2019.
//

#ifndef CODE_IUTILITY_HPP
#define CODE_IUTILITY_HPP

#include <vector>
#include <Dense>

using namespace std;
using namespace Eigen;

class IUtility
{
public:
    virtual vector<double> LinearSpaced (double min, double max, int size) = 0;

    virtual void
    writeToFile (string fNamePrefix, map<double, double> m, int k, vector<string> colNames) = 0;

    virtual void
    writeToFile (string fNamePrefix, map<double, vector<double>> m, int k, vector<string> colNames) = 0;

    virtual void
    writeToFile (string fNamePrefix, map<int, double> m, int k, vector<string> colNames) = 0;

    virtual void
    writeToFile (string fNamePrefix, VectorXd vxd, int k, vector<string> colNames) = 0;

    virtual void
    writeToFile (string fNamePrefix, ArrayXXd axd, int k, vector<string> colNames) = 0;

    virtual void
    writeToFile (string fNamePrefix, MatrixXd axd, int k) = 0;

};

#endif //CODE_IUTILITY_HPP
