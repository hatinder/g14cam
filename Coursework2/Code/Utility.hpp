//
// Created by hsingh9 on 17/04/2019.
//

#ifndef CODE_UTILITY_HPP
#define CODE_UTILITY_HPP


#include <map>
#include "IUtility.hpp"

class Utility :
        public IUtility
{
public:
    vector<double> LinearSpaced (double min, double max, int size) override;

    void writeToFile (string fNamePrefix, map<double, double> m, int k, vector<string> colNames) override;

    void writeToFile (string fNamePrefix, MatrixXd mxd, int k) override;

    void writeToFile (string fNamePrefix, ArrayXXd axd, int k, vector<string> colNames) override;

    void writeToFile (string fNamePrefix, VectorXd vxd, int k, vector<string> colNames) override;

    void writeToFile (string fNamePrefix, map<int, double> m, int k, vector<string> colNames) override;
};


#endif //CODE_UTILITY_HPP
