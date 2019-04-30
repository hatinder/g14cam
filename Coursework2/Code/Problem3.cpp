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
    double res = cos(4 * M_PI * x) * cos(4 * M_PI * y);
//        cout<<"u("<<x<<","<<y<<"): "<<res<<endl;
    return res;
};

auto gB = [] (double x, double y)
{
    double res = (y == 1) ? 1.0 : 0.0;
    return res;
};

auto f = [] (double x, double y)
{ return 32 * M_PI * M_PI * cos(4 * M_PI * x) * cos(4 * M_PI * y); };

VectorXd Problem3::findU (SparseMatrix<double> SpA, VectorXd F)
{
    VectorXd uVec = VectorXd::Zero(SpA.rows());
//    cout << SpA << endl;

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(SpA);
    if (solver.info() != Success)
    {
        cout << "Decomposition Failed. " << endl;
        return uVec;
    }
    setNbThreads(16);
    cout<<"Threads used: "<<nbThreads()<<endl;
//    ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver(SpA);
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
    vector<string> colNames = {"N", "ERROR"};
    vector<int> N = {16, 32, 64, 128}; //,256,512};
    vector<double> error(N.size(), 0.0);
    double a = 0, b = 1;
    Utility utility;
    StokesPDE sPDE;
    clock_t begin, end;
    for (unsigned k = 0; k < N.size(); ++k)
    {
        cout << "Running for N: " << N[k] << endl;
        begin = clock();
        VectorXd nodePoint = VectorXd::LinSpaced(N[k] + 1, 0, 1);
        SparseMatrix<double> SpA = sPDE.createA(N[k]);
        VectorXd F = sPDE.createB(f, gA, N[k], a, b);
        VectorXd uHat = findU(SpA, F);
        MatrixXd Z = MatrixXd::Zero(N[k] + 1, N[k] + 1);
        for (int i = 0; i < N[k] + 1; ++i)       //column
        {
            for (int j = 0; j < N[k] + 1; ++j)       //row
            {
//            cout<<"i*h: "<<nodePoint(i)<<" , j*h: "<<nodePoint(j)<<endl;
//            cout<<"i*(N-1)+j: "<<(i-1)*(N-1)+(j-1)<<endl;
                if (i == 0 || i == N[k] || j == 0 || j == N[k])
                {
                    Z(j, i) = gA(nodePoint(i), nodePoint(j));
                }
                else
                {
                    Z(j, i) = uHat[(i - 1) * (N[k] - 1) + (j - 1)];
                }
            }
//        cout<<"Z"<<endl<<Z<<endl;
        }
        end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken for N: " << N[k] << " is :" << elapsed_secs << endl;
        utility.writeToFile("VECTOR", nodePoint, N[k]);
        utility.writeToFile("MATRIX", Z, N[k]);
        VectorXd uOrig = VectorXd::Zero((N[k] - 1) * (N[k] - 1));
        for (int l = 0; l < (N[k] - 1); ++l)
        {
            for (int i = 0; i < (N[k] - 1); ++i)
            {
                uOrig(l * (N[k] - 1) + i) = gA(nodePoint(i + 1), nodePoint(l + 1));
            }
        }
        utility.writeToFile("VECTORUHAT", uHat, N[k]);
        utility.writeToFile("VECTORU", uOrig, N[k]);
        error[k] = (uOrig - uHat).norm();
    }
    map<int, double> results;
    for (unsigned m = 0; m < error.size(); ++m)
    {
        results[N[m]] = error[m];
    }
    utility.writeToFile("PRB3AERROR", results, 0, colNames);
}

void Problem3::B ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 3 (B)" << endl;
    cout << "=====================" << endl;
    int N = 4;
    double a = 0, b = 1;
//    double h = (b - a) / N;
    StokesPDE stokesPde;
    SparseMatrix<double> A = stokesPde.createA(N);                  // creates A (N-1) x (N-1)
    SparseMatrix<double> Z = stokesPde.createZ(N);                  // creates 0 (N-1) x (N-1)
    SparseMatrix<double> ZN = stokesPde.createZN(N);                // creates 0 (N) x (N)
    SparseMatrix<double> Bx = stokesPde.createBx(N);                // creates Bx (N-1) x (N)
    SparseMatrix<double> By = stokesPde.createBy(N);                // creates By (N-1) x (N)
    SparseMatrix<double> C = stokesPde.createC(A, Bx, By, Z, ZN);   // creates By 2x(N-1) + (N) X 2x(N-1) + (N)
//    cout<<A<<endl;
//    cout<<Z<<endl;
//    cout<<C<<endl;
//    cout<<stokesPde.createBy(N)<<endl;
    VectorXd Fu = stokesPde.createB(f, gB, N, a, b); //Creating F_u gB=function based on boundary y=1
    VectorXd Fv = VectorXd::Zero((N - 1) * (N - 1));
    VectorXd Fp = VectorXd::Zero(N * N);
    VectorXd F = stokesPde.createF(Fu, Fv, Fp);
    clock_t begin,end;
    begin = clock();
    VectorXd U = findBigU(C, F);
    end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time Taken to Solve for N: " << N << " is :" << elapsed_secs << endl;
//    cout<<U<<endl;
    Utility utility;
    VectorXd nodePoint = VectorXd::LinSpaced(N + 1, 0, 1);
    MatrixXd matU = MatrixXd::Zero(N + 1, N + 1);
    for (int i = 0; i < N + 1; ++i)       //column
    {
        for (int j = 0; j < N + 1; ++j)       //row
        {
            if (i == 0 || i == N || j == 0 || j == N)
            {
                matU(j, i) = gB(nodePoint(i), nodePoint(j));
            }
            else
            {
                matU(j, i) = U[(i - 1) * (N - 1) + (j - 1)];
            }
        }
    }
    utility.writeToFile("PRB3BMATRIXU",matU,0);
    MatrixXd matV = MatrixXd::Zero(N + 1, N + 1);
    for (int i = 0; i < N + 1; ++i)       //column
    {
        for (int j = 0; j < N + 1; ++j)       //row
        {
            if (i == 0 || i == N || j == 0 || j == N)
            {
                matV(j, i) = 0;
            }
            else
            {
                matV(j, i) = U[(N - 1) * (N - 1) + (i - 1) * (N - 1) + (j - 1)];
            }
        }
    }
    utility.writeToFile("PRB3BMATRIXV",matV,0);
}

VectorXd Problem3::findBigU (SparseMatrix<double> SpA, VectorXd F)
{
    VectorXd uVec = VectorXd::Zero(SpA.rows());
//    cout << SpA << endl;
    SparseQR <SparseMatrix<double>, COLAMDOrdering<int> > solver;
    solver.compute(SpA);
    if (solver.info() != Success)
    {
        cout << "Decomposition Failed. " << endl;
        return uVec;
    }
//    ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver(SpA);
    uVec = solver.solve(F);
    if (solver.info() != Success)
    {
        cout << "Solving Failed. " << endl;
        return uVec;
    }
    return uVec;
//    cout << Z << endl;
}
