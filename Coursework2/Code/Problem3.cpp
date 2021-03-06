//
// Created by hsingh9 on 24/04/2019.
//
#include <Dense>
#include <Sparse>
#include <iostream>
#include "Problem3.hpp"
#include "StokesPDE.hpp"
#include "Utility.hpp"
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

auto fA = [] (double x, double y)
{ return 32 * M_PI * M_PI * cos(4 * M_PI * x) * cos(4 * M_PI * y); };

auto fB = [] (double x, double y)
{
    double res = 0.0;
    return res;
};

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
//    setNbThreads(16);
//    cout << "Threads used: " << nbThreads() << endl;
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
    vector<int> N = {16, 32, 64, 128 , 256, 512,1024};
    vector<double> error(N.size(), 0.0);
    double a = 0, b = 1;
    Utility utility;
    StokesPDE sPDE;
    clock_t begin, end;
    for (unsigned k = 0; k < N.size(); ++k)
    {
        cout << "Running for N: " << N[k] << endl;
        VectorXd nodePoint = VectorXd::LinSpaced(N[k] + 1, 0, 1);
        SparseMatrix<double> SpA = sPDE.createANew(N[k]);
        VectorXd F = sPDE.createB(fA, gA, N[k], a, b);
        begin = clock();
        VectorXd uHat = findU(SpA, F);
        end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Time Taken to find U for N: " << N[k] << " is :" << elapsed_secs << endl;
        MatrixXd Z = MatrixXd::Zero(N[k] + 1, N[k] + 1);
        for (int i = 0; i < N[k] + 1; ++i)       //column
        {
            for (int j = 0; j < N[k] + 1; ++j)       //row
            {
                if (i == 0 || i == N[k] || j == 0 || j == N[k])
                {
                    Z(j, i) = gA(nodePoint(i), nodePoint(j));
                }
                else
                {
                    Z(j, i) = uHat[(i - 1) * (N[k] - 1) + (j - 1)];
                }
            }
        }
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
    int N = 64  ;
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
    VectorXd Fu = stokesPde.createBU(gB, N, a, b); //Creating F_u gB=function based on boundary y=1
//    cout<<"Fu: "<<Fu<<endl;
    VectorXd Fv = VectorXd::Zero((N - 1) * (N - 1));
    VectorXd Fp = VectorXd::Zero(N * N);
    VectorXd F = stokesPde.createF(Fu, Fv, Fp);
    clock_t begin, end;
    begin = clock();
    VectorXd U = findBigU(C, F);
    end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time Taken to Solve for N: " << N << " is :" << elapsed_secs << endl;
//    cout<<U<<endl;
    Utility utility;
//    utility.writeToFile("PRB3BMATRIXC",C,0);
//    utility.writeToFile("PRB3BMATRIXF",F,0);
//    utility.writeToFile("PRB3BVECTORU",U,0);
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
    begin=clock();
    utility.writeToFile("PRB3BMATRIXU", matU, 0);
    end=clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time Taken to Write Matrix U for N: " << N << " is :" << elapsed_secs << endl;
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
    begin=clock();
    utility.writeToFile("PRB3BMATRIXV", matV, 0);
    end=clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time Taken to Write Matrix V for N: " << N << " is :" << elapsed_secs << endl;
}

VectorXd Problem3::findBigU (SparseMatrix<double> SpA, VectorXd F)
{
//    SpA.insertBack(SpA.rows()-1,SpA.cols()-1)=1;
    VectorXd uVec = VectorXd::Zero(SpA.rows());
//    cout << SpA << endl;
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(SpA);
    if (solver.info() != Success)
    {
        cout << "Decomposition Failed. " << endl;
        return uVec;
    }
//    SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
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
