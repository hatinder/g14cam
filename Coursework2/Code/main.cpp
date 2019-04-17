#include <iostream>
#include <iomanip>
#include "GaussQuadrature.hpp"
#include "Utility.hpp"

using namespace std;

template<typename T>
ostream &operator<< (ostream &out, const vector<T> &v);    //cout operator overload for vector
template<typename T>
ostream &operator<< (ostream &out, const map<T, T> &m);    //cout operator overload for map
void runProblem1a ();                                      //PROBLEM 1a Implementation
void runProblem1b ();                                      //PROBLEM 1b Implementation

int main ()
{
    std::cout << "Hello, World!" << std::endl;
    runProblem1a();
    runProblem1b();
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

}

template<typename T>
ostream &operator<< (ostream &out, const vector<T> &v)
{
    for (int i = 0; i < v.size(); ++i)
    { out << setw(12) << v[i] << endl; }
    return out;
}

template<typename T>
ostream &operator<< (ostream &out, const map<T, T> &m)
{
    for (const pair<T, T> p:m)
    { out << setw(12) << p.first << " : " << setw(12) << p.second << endl; }
    return out;
}
