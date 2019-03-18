#include <iostream>
#include <Dense>
#include <cmath>
#include <vector>
#include "CubicSpline.hpp"
#include "CubicSplinePeriodic.hpp"
#include "RootFinding.hpp"

using namespace std;
using namespace Eigen;

void runProblem1c();
void runProblem1f();
void runProblem2b();

int main ()
{
//    runProblem1c();
    //runProblem1f();
    //VectorXd test=VectorXd::LinSpaced(9,0,1);
    //cout<<test<<endl;
    runProblem2b();
    return 0;
}

void runProblem1c()
{

    cout << "Problem 1 a to c" << std::endl;
    vector<string> uValuesNames={"uSpacePoints", "6PIx", "COS(6PIx)", "q3n"};
    vector<string> aeColNames={"n", "ApproxError", "log10AE", "log10n" };
    ArrayXi n(5);
    n <<   8, 16, 32, 64, 128 ;
    ArrayXXd infinityNorm(5 ,4);
    for (int j=0 ; j < n.size(); j++)
    {
        const double a=0,b=1;
        VectorXd nodalPoints=VectorXd::LinSpaced(n[j]+1,a,b);
        VectorXd fValue(n[j]+1);
        for (int i = 0; i < n[j]+1; ++i)
        {
            fValue[i]=cos(6*M_PI*nodalPoints[i]);
        }
//        cout<<"Function Value"<<endl;
//        cout<<fValue<<endl;
//        cout<<"Function Value"<<endl;
        CubicSpline cs;
        VectorXd coefficents=cs.findCoefficients(fValue);
//        cout<<coefficents<<endl;
        ArrayXXd uValues=cs.findSplineValues(coefficents,nodalPoints);
        VectorXd ApproxError=uValues.col(2)-uValues.col(3);
        infinityNorm(j,0)=n[j];
        infinityNorm(j,1)=ApproxError.lpNorm<Infinity>();
        infinityNorm(j,2)=log10(ApproxError.lpNorm<Infinity>());
        infinityNorm(j,3)=log10(n[j]);
        //cout<<"Approx Error: "<< ApproxError.lpNorm<Infinity>() <<endl;
        cs.writeToFile("1C.", uValues, n[j], uValuesNames);
    }
    CubicSpline t;
    t.writeToFile("1C.APPROXERROR.", infinityNorm, 0, aeColNames);

}

void runProblem1f()
{

    cout << "Problem 1 d to e" << std::endl;
    vector<string> uValuesNames={"uSpacePoints", "6PIx", "COS(6PIx)", "q3n"};
    vector<string> aeColNames={"n", "ApproxError", "log10AE", "log10n" };
    ArrayXi n(5);
    n <<  8, 16, 32, 64, 128 ;
    ArrayXXd infinityNorm(5,4);
    for (int j=0 ; j < n.size(); j++)
    {
        const double a=0,b=1;
        const double h=(b-a)/(double)(n[j]);
        //cout<<"a: "<<a<<" ,b: "<<b<<" ,b-h: "<<b-h<<endl;
        VectorXd nodalPoints=VectorXd::LinSpaced(n[j],a,b-h);
        ////cout<<"NODALPOINTS: " <<endl<< nodalPoints<<endl;
        VectorXd fValue(n[j]);
        for (int i = 0; i < n[j]; ++i)
        {
            fValue[i]=cos(6*M_PI*nodalPoints[i]);
        }
        //cout<<"FVALUE: "<<endl<<fValue<<endl;
        CubicSplinePeriodic csp;
        VectorXd coefficents=csp.findCoefficients(fValue);
        //cout<<"COEFFICIENTS: "<<endl<<coefficents<<endl;
        VectorXd splineNodalPoints=VectorXd::LinSpaced(n[j]+1,a,b);
        //cout<<"SplineNodalPoints:"<<endl<<splineNodalPoints<<endl;
        ArrayXXd uValues=csp.findSplineValues(coefficents,splineNodalPoints);
        VectorXd ApproxError=uValues.col(2)-uValues.col(3);
        infinityNorm(j,0)=n[j];
        infinityNorm(j,1)=ApproxError.lpNorm<Infinity>();
        infinityNorm(j,2)=log10(ApproxError.lpNorm<Infinity>());
        infinityNorm(j,3)=log10(n[j]);
        ////cout<<"Approx Error: "<< ApproxError.lpNorm<Infinity>() <<endl;
        csp.writeToFile("1E.", uValues, n[j], uValuesNames);
    }
    CubicSpline t;
    t.writeToFile("1E.APPROXERROR.", infinityNorm, 0, aeColNames);
}

void runProblem2b()
{
    RootFinding rootFinding;
    ArrayXXd initialGuess=ArrayXXd::Zero(8,3);
    ArrayXXd roots=ArrayXXd::Zero(4,3);
    initialGuess<<0.1,0.1,0.1,
            -0.1,0.1,0.1,
            0.1,-0.1,0.1,
            0.1,0.1,-0.1,
            -0.1,-0.1,0.1,
            0.1,-0.1,-0.1,
            -0.1,0.1,-0.1,
            -0.1,-0.1,-0.1;
//    cout<<rootFinding.calculateF(0.1,0.1,0.1)<<endl;
//    cout<<rootFinding.calculateJ(0.1,0.1,0.1)<<endl;
    int nRoots=0;
    for (int i = 0; i < initialGuess.rows(); ++i)
    {
        VectorXd xOld=initialGuess.row(i);
        VectorXd xNew=VectorXd::Zero(3);
        double TOL=pow(10,-12);
        double tolerance=1;
        int j=0;
        for (; j < 100 && tolerance > TOL; ++j)
        {
//            cout<<"xOld: "<< xOld<<endl;
            VectorXd F=rootFinding.calculateF(xOld(0),xOld(1),xOld(2));
            MatrixXd J=rootFinding.calculateJ(xOld(0),xOld(1),xOld(2));
            FullPivLU<MatrixXd> JInv(J);
            VectorXd y=-JInv.solve(F);
            xNew=xOld+y;
//            cout<<"xNew: "<< xNew<<endl;
            tolerance=(xNew-xOld).lpNorm<2>();
            xOld=xNew;
        }
        if(rootFinding.foundNewRoot(xNew(0),xNew(1),xNew(2),roots))
        {
            cout<<"Initial Guess: ("<<initialGuess(i,0)<<", "<<initialGuess(i,1)<<", "<<initialGuess(i,2)<<")"<<endl;
            cout<<"Iteration: "<<j<<", New Root Found: ("<<xNew(0)<<", "<<xNew(1)<<", "<<xNew(2)<<")"<<endl;
            roots.row(nRoots)=xNew;
            nRoots++;
            if(nRoots==4)
                break;
        }
        else
        {
//            cout<<"Initial Guess: ("<<initialGuess(i,0)<<", "<<initialGuess(i,1)<<", "<<initialGuess(i,2)<<")"<<endl;
//            cout<<"Iteration: "<<j<<", Not New Root : ("<<xNew(0)<<", "<<xNew(1)<<", "<<xNew(2)<<")"<<endl;
        }
    }
    cout<<"All Roots: "<<endl<<roots<<endl;



}