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
    int n=N-1;
    double h=(b-a)/(double)N;
    VectorXd v(n*n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double u=0.0;
//            if(j==(n-1))
//                u=g(0,1);
//            else
//                u=g(0,0);
            if(i==0)        //Left Boundary + Top + Bottom Outer
            {
                 if(j==0)
                    u=g(a+(i+1)*h,a+j*h)+g(a+i*h,a+(j+1)*h);
                 else if(j==(n-1))
                    u=g(a+i*h,a+(j+1)*h)+g(a+(i+1)*h,a+(j+2)*h);
                else
                    u=g(a+i*h,a+(j+1)*h);
            }
            else if(i==(n-1))        //Right Boundary Top + Bottom Outer
            {
                if(j==0)
                    u=g(a+(i+1)*h,a+j*h)+g(a+(i+2)*h,a+(j+1)*h);
                else if(j==(n-1))
                    u=g(a+(i+2)*h,a+(j+1)*h)+g(a+(i+1)*h,a+(j+2)*h);
                else
                    u=g(a+(i+2)*h,a+(j+1)*h);
            }
            else if(j==0)        //Bottom Internal
            {
                if(i>0 or i<(n-1))
                    u=g(a+(i+1)*h,a+j*h);
            }
            else if(j==(n-1))    //Top Internal
            {
                if(i>0 or i<(n-1))
                    u=g(a+(i+1)*h,a+(j+2)*h);
            }
            v(i*n+j)=h*h*f(a+(i+1)*h,a+(j+1)*h) + u;
//            cout<<"("<<i<<","<<j<<")"<< " , i*n+j: "<<i*n+j<<" , f: "<<h*h*f(a+(i+1)*h,a+(j+1)*h)<<" ,g:"<<u<<endl;
        }
    }
    return  v;
}
