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
