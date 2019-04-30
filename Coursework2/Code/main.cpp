#include <iostream>
#include <iomanip>
#include "GaussQuadrature.hpp"
#include "Utility.hpp"
#include "helper.hpp"
#include "Problem2.hpp"
#include "Problem3.hpp"

using namespace std;

void runProblem1a ();                                      //PROBLEM 1a Implementation
void runProblem1b ();                                      //PROBLEM 1b Implementation
void runProblem1c ();                                      //PROBLEM 1c Implementation
void runProblem1d ();                                      //PROBLEM 1d Implementation
void runProblem1e ();                                      //PROBLEM 1e Implementation
void runProblem1f ();                                      //PROBLEM 1f Implementation

int main ()
{
    std::cout << "Hello, World!" << std::endl;
//    runProblem1a();
//    runProblem1b();
//    runProblem1c();
//    runProblem1d();
//    runProblem1e();
//    runProblem1f();
//    Problem2 problem2;
//    problem2.ARungeKutta2();
//    problem2.AImplicitMidpoint();
//    problem2.BRungeKutta2();
//    problem2.BImplicitMidpoint();
//    problem2.CImplicitMidpoint();
    Problem3 problem3;
//    problem3.A();
    problem3.B();
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
    GaussQuadrature gaussQuadrature;
/*
    vector<vector<double >> pointsAndWeights=gaussQuadrature.findPointsAndWeights(3);
    cout<<"Points: "<<endl<<pointsAndWeights[0];
    cout<<"Weights: "<<endl<<pointsAndWeights[1];
    cout << "Approx Value: " << gaussQuadrature.getApproxValue1C(1.0, 3.0, pointsAndWeights[0], pointsAndWeights[1], 3, 5)
         << endl;
    cout<<"Exact  Value: "<<gaussQuadrature.getExactValueFor1C(1.0,3.0,5)<<endl;
*/
    double a = 1.0, b = 3.0;
    for (int i = 1; i <= 5; ++i)    //Setting n from 1 to 5
    {
        double exactValue, approxValue;
        for (int j = 1; j <= 2 * i; ++j)    //Setting degree from 1 to 2n
        {
            vector<vector<double >> pointsAndWeights = gaussQuadrature.findPointsAndWeights(i); //get roots, weights
            exactValue = gaussQuadrature.getExactValueFor1C(a, b, j);
            approxValue = gaussQuadrature.getApproxValue1C(a, b, pointsAndWeights[0], pointsAndWeights[1], i, j);
            cout << "Degree: " << setw(3) << j << " , n: " << setw(3) << i << " , Exact: " << setw(12)
                 << setprecision(10) << exactValue << " , Approx: " << setw(12) << approxValue << " , Error: "
                 << setw(12) << abs(exactValue - approxValue) << endl;
        }
    }
}

void runProblem1d ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 1 (d)" << endl;
    cout << "=====================" << endl;
    GaussQuadrature gaussQuadrature;
    vector<string> colNames = {"n", "error"};
    Utility utility;
    double a = 0.0, b = 3 * M_PI / 4.0, approxValue, exactValue;
    exactValue = gaussQuadrature.getExactValueFor1D(a, b);
    map<int,double> results;
    for (int i = 1; i <= 10; ++i)
    {
        vector<vector<double >> pointsAndWeights = gaussQuadrature.findPointsAndWeights(i); //get roots, weights
        approxValue = gaussQuadrature.getApproxValue1D(a, b, pointsAndWeights[0], pointsAndWeights[1], i);
        cout << "n: " << setw(3) << i << " , Exact: " << setw(12) << setprecision(10) << exactValue << " , Approx: "
             << setw(12) << approxValue << " , Error: " << setw(12) << abs(exactValue - approxValue) << endl;
        results[i]=abs(exactValue - approxValue);
    }
    utility.writeToFile("PROBLEM1D",results,0,colNames);
}

void runProblem1e()
{
    cout << "=====================" << endl;
    cout << "Running Problem 1 (e)" << endl;
    cout << "=====================" << endl;
    GaussQuadrature gaussQuadrature;
    double a=0.0,b=1.0,approxValue;
    auto f=[](double x, double y) {return sin(x*x+y*y);};
    for (int i = 1; i <= 5; ++i)
    {
        vector<vector<double >> pointsAndWeights = gaussQuadrature.findPointsAndWeights(i);
        approxValue=gaussQuadrature.getApproxValue1E(a,b,f,i,pointsAndWeights[0],pointsAndWeights[1]);
        cout<<"n: "<<i <<" , approxValue: "<<approxValue<<endl;

    }
}

void runProblem1f()
{
    cout << "=====================" << endl;
    cout << "Running Problem 1 (f)" << endl;
    cout << "=====================" << endl;
    GaussQuadrature gaussQuadrature;
    double a=0.0,b=1.0,approxValue;
    auto f=[](double x, double y) {return sin(x*x+y*y);};   //input function using lambda feature of C++ 11
    for (int i = 1; i <= 5; ++i)
    {
        vector<vector<double >> pointsAndWeights = gaussQuadrature.findPointsAndWeights(i);
        approxValue=gaussQuadrature.getApproxValue1F(a,b,f,i,pointsAndWeights[0],pointsAndWeights[1]);
        cout<<"n: "<<i <<" , approxValue: "<<approxValue<<endl;

    }
}