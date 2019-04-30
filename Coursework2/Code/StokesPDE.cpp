//
// Created by hsingh9 on 24/04/2019.
//

#include "StokesPDE.hpp"
#include <iostream>

using namespace std;

SparseMatrix<double> StokesPDE::createA (int N)
{
    int n = N - 1;
    SparseMatrix<double> A(n * n, n * n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; ++j)
        {
            A.insert(i + j * n, i + j * n) = 4;
            if ((i + 1) < n)
            {
                A.insert(i + j * n, i + 1 + j * n) = -1;
            }
            if ((i - 1) >= 0)
            {
                A.insert(i + j * n, i - 1 + j * n) = -1;
            }
            if ((j + 1) < n)
            {
                A.insert(i + j * n, i + (j + 1) * n) = -1;
            }
            if ((j - 1) >= 0)
            {
                A.insert(i + j * n, i + (j - 1) * n) = -1;
            }
        }
//        cout<<"i: "<<i<<"j: "<<n-1<<endl;
//        cout<<A<<endl;
    }
    return A;
}

VectorXd StokesPDE::createB (double (*f) (double, double), double (*g) (double, double), int N, double a, double b)
{
    int n = N - 1;
    double h = (b - a) / (double) N;
    VectorXd v(n * n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double u = 0.0;
//            if(j==(n-1))
//                u=g(0,1);
//            else
//                u=g(0,0);
            if (i == 0)        //Left Boundary + Top + Bottom Outer
            {
                if (j == 0)
                {
                    u = g(a + (i + 1) * h, a + j * h) + g(a + i * h, a + (j + 1) * h);
                }
                else if (j == (n - 1))
                {
                    u = g(a + i * h, a + (j + 1) * h) + g(a + (i + 1) * h, a + (j + 2) * h);
                }
                else
                {
                    u = g(a + i * h, a + (j + 1) * h);
                }
            }
            else if (i == (n - 1))        //Right Boundary Top + Bottom Outer
            {
                if (j == 0)
                {
                    u = g(a + (i + 1) * h, a + j * h) + g(a + (i + 2) * h, a + (j + 1) * h);
                }
                else if (j == (n - 1))
                {
                    u = g(a + (i + 2) * h, a + (j + 1) * h) + g(a + (i + 1) * h, a + (j + 2) * h);
                }
                else
                {
                    u = g(a + (i + 2) * h, a + (j + 1) * h);
                }
            }
            else if (j == 0)        //Bottom Internal
            {
                if (i > 0 or i < (n - 1))
                {
                    u = g(a + (i + 1) * h, a + j * h);
                }
            }
            else if (j == (n - 1))    //Top Internal
            {
                if (i > 0 or i < (n - 1))
                {
                    u = g(a + (i + 1) * h, a + (j + 2) * h);
                }
            }
            v(i * n + j) = h * h * f(a + (i + 1) * h, a + (j + 1) * h) + u;
//            cout<<"("<<i<<","<<j<<")"<< " , i*n+j: "<<i*n+j<<" , f: "<<h*h*f(a+(i+1)*h,a+(j+1)*h)<<" ,g:"<<u<<endl;
        }
    }
    return v;
}

SparseMatrix<double> StokesPDE::createBx (int N)
{
    int n = N - 1;
    SparseMatrix<double> Bx(n * n, N * N);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < N; ++j)
        {
            if (j < (N - 1))
            {
                Bx.insert(i + j * n, i + j * N) = -1;
                Bx.insert(i + j * n, i + 1 + j * N) = -1;
                Bx.insert(i + j * n, i + N + j * N) = 1;
                Bx.insert(i + j * n, i + N + 1 + j * N) = 1;
            }
        }
//        cout<<"i: "<<i<<"j: "<<n-1<<endl;
//        cout<<A<<endl;
    }
    return Bx;

}

SparseMatrix<double> StokesPDE::createBy (int N)
{
    int n = N - 1;
    SparseMatrix<double> By(n * n, N * N);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < N; ++j)
        {
            if (j < (N - 1))
            {
                By.insert(i + j * n, i + j * N) = -1;
                By.insert(i + j * n, i + 1 + j * N) = 1;
                By.insert(i + j * n, i + N + j * N) = -1;
                By.insert(i + j * n, i + N + 1 + j * N) = 1;
            }
        }
//        cout<<"i: "<<i<<"j: "<<n-1<<endl;
//        cout<<A<<endl;
    }
    return By;
}

SparseMatrix<double> StokesPDE::createZ (int N)
{
    int n = N - 1;
    SparseMatrix<double> A(n * n, n * n);
    return A;
}

SparseMatrix<double>
StokesPDE::createC (SparseMatrix<double> A, SparseMatrix<double> Bx, SparseMatrix<double> By, SparseMatrix<double> Z,
                    SparseMatrix<double> ZN)
{
    SparseMatrix<double> BxT = Bx.transpose();
    SparseMatrix<double> ByT = By.transpose();
    SparseMatrix<double> C(A.rows() + A.rows() + ByT.rows(), A.cols() + Z.cols() + Bx.cols());
    C.reserve(
            A.nonZeros() + Z.nonZeros() + Bx.nonZeros() + Z.nonZeros() + A.nonZeros() + By.nonZeros() + BxT.nonZeros() +
            ByT.nonZeros() + ZN.nonZeros());
    // Set First Column  i.e. A , Z, BxT
    for (Index c = 0; c < A.cols(); ++c)
    {
        C.startVec(c);
        for (SparseMatrix<double>::InnerIterator itA(A, c); itA; ++itA)
        {
            C.insertBack(itA.row(), c) = itA.value();
        }
        for (SparseMatrix<double>::InnerIterator itZ(Z, c); itZ; ++itZ)
        {
            C.insertBack(itZ.row() + A.rows(), c) = itZ.value();
        }
        for (SparseMatrix<double>::InnerIterator itBxT(BxT, c); itBxT; ++itBxT)
        {
            C.insertBack(itBxT.row() + A.rows() + Z.rows(), c) = itBxT.value();
        }
    }

    // Set Second Column  i.e. 0 , A, ByT

    for (Index c = 0; c < Z.cols(); ++c)
    {
        C.startVec(A.cols() + c);
        for (SparseMatrix<double>::InnerIterator itZ(Z, c); itZ; ++itZ)
        {
            C.insertBack(itZ.row(), A.cols() + c) = itZ.value();
        }
        for (SparseMatrix<double>::InnerIterator itA(A, c); itA; ++itA)
        {
            C.insertBack(itA.row() + Z.rows(), A.cols() + c) = itA.value();
        }
        for (SparseMatrix<double>::InnerIterator itByT(ByT, c); itByT; ++itByT)
        {
            C.insertBack(itByT.row() + Z.rows() + A.rows(), A.cols() + c) = itByT.value();
        }
    }
//
//    // Set Third Column  i.e. Bx , By, 0
//
    for (Index c = 0; c < Bx.cols(); ++c)
    {
        C.startVec(A.cols() + Z.cols() + c);
        for (SparseMatrix<double>::InnerIterator itBx(Bx, c); itBx; ++itBx)
        {
            C.insertBack(itBx.row(), A.cols() + Z.cols() + c) = itBx.value();
        }
        for (SparseMatrix<double>::InnerIterator itBy(By, c); itBy; ++itBy)
        {
            C.insertBack(itBy.row() + Z.rows(), A.cols() + Z.cols() + c) = itBy.value();
        }
        for (SparseMatrix<double>::InnerIterator itZN(ZN, c); itZN; ++itZN)
        {
            C.insertBack(itZN.row() + Z.rows() + A.rows(), A.cols() + Z.cols() + c) = itZN.value();
        }
    }

    C.finalize();
    return C;
}

SparseMatrix<double> StokesPDE::createZN (int N)
{
    SparseMatrix<double> ZN(N * N, N * N);
    return ZN;
}

VectorXd StokesPDE::createF (VectorXd Fu, VectorXd Fv, VectorXd Fp)
{
    VectorXd F(Fu.size() + Fv.size() + Fp.size());
    F << Fu, Fv, Fp;
    return F;
}

VectorXd StokesPDE::createBU (double (*g) (double, double), int N, double a, double b)
{
    int n = N - 1;
    double h = (b - a) / (double) N;
    VectorXd v(n * n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double u = 0.0;
//            if(j==(n-1))
//                u=g(0,1);
//            else
//                u=g(0,0);
            if (i == 0)        //Left Boundary + Top + Bottom Outer
            {
                if (j == 0)
                {
                    u = g(a + (i + 1) * h, a + j * h) + g(a + i * h, a + (j + 1) * h);
                }
                else if (j == (n - 1))
                {
                    u = g(a + i * h, a + (j + 1) * h) + g(a + (i + 1) * h, a + (j + 2) * h);
                }
                else
                {
                    u = g(a + i * h, a + (j + 1) * h);
                }
            }
            else if (i == (n - 1))        //Right Boundary Top + Bottom Outer
            {
                if (j == 0)
                {
                    u = g(a + (i + 1) * h, a + j * h) + g(a + (i + 2) * h, a + (j + 1) * h);
                }
                else if (j == (n - 1))
                {
                    u = g(a + (i + 2) * h, a + (j + 1) * h) + g(a + (i + 1) * h, a + (j + 2) * h);
                }
                else
                {
                    u = g(a + (i + 2) * h, a + (j + 1) * h);
                }
            }
            else if (j == 0)        //Bottom Internal
            {
                if (i > 0 or i < (n - 1))
                {
                    u = g(a + (i + 1) * h, a + j * h);
                }
            }
            else if (j == (n - 1))    //Top Internal
            {
                if (i > 0 or i < (n - 1))
                {
                    u = g(a + (i + 1) * h, a + (j + 2) * h);
                }
            }
            v(i * n + j) = u;
//            cout<<"("<<i<<","<<j<<")"<< " , i*n+j: "<<i*n+j<<" , f: "<<h*h*f(a+(i+1)*h,a+(j+1)*h)<<" ,g:"<<u<<endl;
        }
    }
    return v;
}
