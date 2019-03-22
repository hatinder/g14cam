//
// Created by hsingh9 on 18/03/2019.
//
#include <fstream>
#include <iomanip>
#include "RootFinding.hpp"

VectorXd RootFinding::calculateF (double x, double y, double z)
{
    VectorXd F=VectorXd::Zero(3);
    F[0]=pow(x,2)+pow(y,2)+pow(z,2)-1;
    F[1]=4*(pow(x,2)+pow(y,2))-pow(z,2);
    F[2]=3*x-y;
    return F;
}

MatrixXd RootFinding::calculateJ (double x, double y, double z)
{
    MatrixXd J=MatrixXd::Zero(3,3);
    J(0,0)=2*x;
    J(0,1)=2*y;
    J(0,2)=2*z;
    J(1,0)=8*x;
    J(1,1)=8*y;
    J(1,2)=-2*z;
    J(2,0)=3;
    J(2,1)=-1;
    J(2,2)=0;
    return J;
}

bool RootFinding::foundNewRoot (double x, double y, double z, ArrayXXd root)
{
    for (int i = 0; i < root.rows(); ++i)
    {
        if(root(i,0)==x && root(i,1)==y && root(i,2)==z)
            return false;
    }
    return true;
}

bool RootFinding::foundNewRoot (double x, double y, ArrayXXd root)
{
    for (int i = 0; i < root.rows(); ++i)
    {
        if(root(i,0)==x && root(i,1)==y)
            return false;
    }
    return true;
}

VectorXd RootFinding::calculateF2C (double x, double y, double x0, double y0,double l2Norm)
{
    VectorXd F = VectorXd::Zero(2);
    F[0] = pow(x, 2) + pow(y, 2)/4 - 1;
    //F[1] = (-4 * pow(x, 2))/(sqrt(1-pow(x,2))) - 0.05;
    //F[1] = y-0.1031;
    //F[1] = 2*x0*(x-x0)+ (0.5)*y0*(y-y0)-l2Norm*0.05;   //TODO : Pending Review
    F[1] = y0*(x-x0) -x0*(y-y0)-l2Norm*0.05;   //TODO : Pending Review
    return F;
}

MatrixXd RootFinding::calculateJ2C (double x, double y, double x0, double y0)
{
    MatrixXd J=MatrixXd::Zero(2,2);
    J(0,0)=2*x;
    J(0,1)=y/2;
//    J(1,0)=(8*x)/(pow(pow(x,2)-1,2));
    //J(1,0)=((-4)*pow(x,3)*sqrt(-pow(x,2)+1)+8*x*pow(-pow(x,2)+1,1.5))/(pow((pow(x,2)-1),2));
    J(1,0)=y0;
    J(1,1)=-x0;
    return J;
}

VectorXd RootFinding::calculateF2D (VectorXd u,double h,double lambda)
{

    VectorXd F=VectorXd::Zero(u.size()-2);
    for (int i = 1; i < u.size()-1; ++i)    //i =1 => u(1)
    {
        F[i-1]=u[i-1]-2*u[i]+u[i+1]+pow(h,2)*lambda*exp(u[i]);
    }
    return F;
}

MatrixXd RootFinding::calculateJ2D (VectorXd u, double h, double lambda)
{
    MatrixXd J=MatrixXd::Zero(u.size()-2,u.size()-2);
    for (int i = 1; i < u.size()-2; ++i)
    {
        J(i-1,i)=1;
        J(i-1,i-1)=-2+pow(h,2)*lambda*exp(u[i]);
        J(i,i-1)=1;
    }
    J(u.size()-3,u.size()-3)=-2+pow(h,2)*lambda*exp(u[u.size()-3]);
    return J;
}


void RootFinding::writeToFile (const string fNamePrefix, ArrayXXd roots, const int k, vector<string> colNames)
{

    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    for(auto v:colNames)
        oFileStream<<setw(12)<<v;
    oFileStream<<endl;
    oFileStream << roots <<endl;
    oFileStream.close();
}

MatrixXd RootFinding::createFuLambda2E(VectorXd u, double h, double lambda)
{
    MatrixXd J=MatrixXd::Zero(u.size()-2,u.size()-2+1);
    for (int i = 1; i < u.size()-2; ++i)
    {
        J(i-1,i)=1;
        J(i-1,i-1)=-2*u[i]+pow(h,2)*lambda*exp(u[i]);
//        J(i-1,i-1)=-2;
        J(i,i-1)=1;
        J(i-1,u.size()-2)=pow(h,2)*lambda*exp(u[i]);
    }
    J(u.size()-3,u.size()-3)=-2*u[u.size()-2]+pow(h,2)*lambda*exp(u[u.size()-2]);
    J(u.size()-3,u.size()-2)=pow(h,2)*lambda*exp(u[u.size()-2]);
    return J;

}

VectorXd RootFinding::calculateF2E (VectorXd u, double h, double lambda, VectorXd u0, VectorXd tk, double lambdaTilda,
                                    double lambda0,double ds)
{
    VectorXd F=VectorXd::Zero(u.size()-1);
    for (int i = 1; i < u.size()-1; ++i)    //i =1 => u(1)
    {
        F[i-1]=u[i-1]-2*u[i]+u[i+1]+pow(h,2)*lambda*exp(u[i]);
    }
    F[u.size()-2]=(u-u0).segment(1,u.size()-2).dot(tk) + lambdaTilda*(lambda-lambda0)-ds;
    return F;

}

MatrixXd RootFinding::calculateJ2E (VectorXd u, double h, double lambda, VectorXd tk, double lambdaTilda)
{
    MatrixXd J=MatrixXd::Zero(u.size()-1,u.size()-1);
    for (int i = 1; i < u.size()-2; ++i)
    {
        J(i-1,i)=1;
        J(i-1,i-1)=-2+pow(h,2)*lambda*exp(u[i]);
        J(i,i-1)=1;
        J(i-1,u.size()-2)=pow(h,2)*exp(u[i]);
        J(u.size()-2,i-1)=tk[i-1];
    }
    J(u.size()-3,u.size()-3)=-2+pow(h,2)*lambda*exp(u[u.size()-3]);
    J(u.size()-3,u.size()-2)=pow(h,2)*exp(u[u.size()-2]);
    J(u.size()-2,u.size()-3)=tk[tk.size()-1];
    J(u.size()-2,u.size()-2)=lambdaTilda;

    return J;
}

VectorXd RootFinding::findInitialU (double lambda, int n)
{
    double h = (double)1/n;
    ArrayXXd roots = ArrayXXd::Zero(n + 1, 2);
    roots.col(0) = ArrayXd::LinSpaced(n + 1, 0, 1);
    VectorXd uInitial = VectorXd::Zero(n + 1);
    double TOL = pow(10, -12);
    double tolerance = 1;
    int i = 0;
    for (; i < 100 && tolerance > TOL; ++i)
    {
        //        cout<<"uInitial: "<<endl<<uInitial<<endl;
        VectorXd uNew = VectorXd::Zero(n + 1);
        VectorXd y = VectorXd::Zero(n - 1);
        VectorXd F_U = calculateF2D(uInitial, h, lambda);
        //        cout<<"F_U: "<<endl<<F_U<<endl;
        MatrixXd J_U = calculateJ2D(uInitial, h, lambda);
        //        cout<<"J_U: "<<endl<<J_U<<endl;
        PartialPivLU<MatrixXd> JInv(J_U);
        y = -JInv.solve(F_U);
        uNew.segment(1, n - 1) = uInitial.segment(1, n - 1) + y.segment(0, n - 1);
        tolerance = (uNew - uInitial).lpNorm<2>();
        //        cout<<"tolerance: "<< tolerance<<endl;
        uInitial = uNew;
    }
    cout << "Iteration: " << i << " for lambda: " << lambda << endl;
    return uInitial;
}
