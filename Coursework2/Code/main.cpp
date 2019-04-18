#include <iostream>
#include <iomanip>
#include "GaussQuadrature.hpp"
#include "Utility.hpp"
#include "helper.hpp"

using namespace std;

void runProblem1a ();                                      //PROBLEM 1a Implementation
void runProblem1b ();                                      //PROBLEM 1b Implementation
void runProblem1c ();

int main ()
{
    std::cout << "Hello, World!" << std::endl;
//    runProblem1a();
//    runProblem1b();
    runProblem1c();
    return 0;
}

void runProblem1a ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 1 (a)" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"x", "y"};
    Utility utils;
    vector<double> nVec = utils.LinearSpaced(-1.0, 1.0, 201);
//    cout << nVec;
    GaussQuadrature gaussQuadrature;
    for (int i = 0; i < 5; ++i)
    {
        map<double, double> results = gaussQuadrature.LegendrePolynomial(i, nVec);
        utils.writeToFile("PROBLEM1A", results, i, colNames);
    }
//    cout<<gaussQuadrature.LegendrePolynomialDerivative(3,nVec);   //Derivatives Output, verified and disabled.
}

void runProblem1b ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 1 (b)" << endl;
    cout << "=====================" << endl;
    GaussQuadrature gaussQuadrature;
    cout << gaussQuadrature.NewtonMethod(5);

}

void runProblem1c ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 1 (c)" << endl;
    cout << "=====================" << endl;

}

