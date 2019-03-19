//
// Created by hsingh9 on 18/03/2019.
//

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

VectorXd RootFinding::calculateF2C (double x, double y, double x0, double y0)
{
    VectorXd F = VectorXd::Zero(2);
    F[0] = pow(x, 2) + pow(y, 2)/4-1;
    //F[1] = (-4 * pow(x, 2))/(sqrt(1-pow(x,2))) - 0.05;
    //F[1] = y-0.1031;
    F[1] = 2*x0*(x-x0)+ (0.5)*y0*(y-y0)-0.05;   //TODO : Pending Review
    return F;
}

MatrixXd RootFinding::calculateJ2C (double x, double y, double x0, double y0)
{
    MatrixXd J=MatrixXd::Zero(2,2);
    J(0,0)=2*x;
    J(0,1)=y/2;
//    J(1,0)=(8*x)/(pow(pow(x,2)-1,2));
    //J(1,0)=((-4)*pow(x,3)*sqrt(-pow(x,2)+1)+8*x*pow(-pow(x,2)+1,1.5))/(pow((pow(x,2)-1),2));
    J(1,0)=2*x0;
    J(1,1)=(0.5)*y0;
    return J;
}
