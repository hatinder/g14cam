//
// Created by hsingh9 on 24/04/2019.
//
#include <Dense>
#include <Core>
#include <Sparse>
#include <iostream>
#include "Problem3.hpp"
#include "StokesPDE.hpp"
#include "Utility.hpp"
#include <SparseCholesky>
#include <ctime>

using namespace std;
using namespace Eigen;

auto gA = [] (double x, double y)
{
    double res=cos(4*M_PI*x)*cos(4*M_PI*y);
//        cout<<"u("<<x<<","<<y<<"): "<<res<<endl;
    return res;
};

VectorXd Problem3::findU (int N)
{
    initParallel();
    VectorXd uVec=VectorXd::Zero((N-1)*(N-1));
    double a=0,b=1;
    double h=(b-a)/N;
    StokesPDE sPDE;
    SparseMatrix<double> SpA = sPDE.createA(N);
//    cout << SpA << endl;
    auto f = [] (double x, double y)
    { return 32 * M_PI *M_PI* cos(4 * M_PI * x) * cos(4 * M_PI * y); };

    VectorXd F = sPDE.createB(f, gA, N, a, b);
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(SpA);
    if (solver.info() != Success)
    {
        cout << "Decomposition Failed. " << endl;
        return uVec;
    }
    uVec = solver.solve(F);
    if (solver.info() != Success)
    {
        cout << "Solving Failed. " << endl;
        return uVec;
    }
    return uVec;
//    cout << Z << endl;

}

void Problem3::A ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 3 (A)" << endl;
    cout << "=====================" << endl;
    vector<string> colNames={"N","ERROR"};
    vector<int> N={16,32,64,128,256,512};
    vector<double> error(N.size(),0.0);
    Utility utility;
    clock_t begin,end;
    for (int k = 0; k < N.size(); ++k)
    {
        cout<<"Running for N: "<<N[k]<<endl;
        begin = clock();
        VectorXd nodePoint=VectorXd::LinSpaced(N[k]+1,0,1);
        VectorXd uHat=findU(N[k]);
        MatrixXd Z = MatrixXd::Zero(N[k]+1,N[k]+1);
        for (int i = 0; i < N[k]+1; ++i)       //column
        {
            for (int j = 0; j < N[k]+1; ++j)       //row
            {
//            cout<<"i*h: "<<nodePoint(i)<<" , j*h: "<<nodePoint(j)<<endl;
//            cout<<"i*(N-1)+j: "<<(i-1)*(N-1)+(j-1)<<endl;
                if(i==0 || i ==N[k] || j==0 || j==N[k])
                    Z(j,i)=gA(nodePoint(i),nodePoint(j));
                else
                    Z(j,i)=uHat[(i-1)*(N[k]-1)+(j-1)];
            }
//        cout<<"Z"<<endl<<Z<<endl;
        }
        end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for N: " << N[k] << " is :" << elapsed_secs << endl;
        utility.writeToFile("VECTOR",nodePoint,N[k]);
        utility.writeToFile("MATRIX",Z,N[k]);
        VectorXd uOrig=VectorXd::Zero((N[k]-1)*(N[k]-1));
        for (int l = 0; l < (N[k]-1); ++l)
        {
            for (int i = 0; i < (N[k]-1); ++i)
            {
                uOrig(l*(N[k]-1)+i)=gA(nodePoint(i+1),nodePoint(l+1));
            }
        }
        utility.writeToFile("VECTORUHAT",uHat,N[k]);
        utility.writeToFile("VECTORU",uOrig,N[k]);
        error[k]=(uOrig-uHat).norm();
    }
    map<int,double> results;
    for (int m = 0; m < error.size(); ++m)
    {
        results[N[m]]=error[m];
    }
    utility.writeToFile("PRB3AERROR",results,0,colNames);
}
