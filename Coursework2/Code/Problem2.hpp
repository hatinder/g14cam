//
// Created by hsingh9 on 20/04/2019.
//

#ifndef CODE_PROBLEM2_HPP
#define CODE_PROBLEM2_HPP

#include <iostream>
#include <map>
#include <vector>

using namespace std;

class Problem2
{
public:
    void ARungeKutta2 ();
    void AImplicitMidpoint();
    void BRungeKutta2 ();
    void BImplicitMidpoint();
    void CImplicitMidpoint ();
    map<double, double> computeHamilton (vector<map<double, double>> yVector, double T, double N);

};


#endif //CODE_PROBLEM2_HPP
