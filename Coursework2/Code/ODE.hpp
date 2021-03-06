//
// Created by hsingh9 on 22/04/2019.
//

#ifndef CODE_ODE_HPP
#define CODE_ODE_HPP

#include <vector>
#include <map>

using namespace std;

class ODE
{
public:

    map<double,vector<double>>
    applyRK2 (vector<double> (*F) (double, vector<double>), vector<double> initVal, double T, double N);

    map<double,vector<double>>
    applyIMP (vector<double> (*F) (double, vector<double>), vector<double> initVal, double T, double N);

};


#endif //CODE_ODE_HPP
