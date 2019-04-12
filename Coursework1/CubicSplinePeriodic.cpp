//
// Created by hsingh9 on 15/03/2019.
//

#include "CubicSplinePeriodic.hpp"

VectorXd CubicSplinePeriodic::findCoefficients (VectorXd fValue)
{
    const int mSize=fValue.size();
    VectorXd cValue=VectorXd::Zero(mSize);
    VectorXd localFValue=VectorXd::Zero(mSize);
    for (int i = 0; i < mSize ; ++i)
    {
        localFValue[i]=fValue[i];
    }
    MatrixXd M=MatrixXd::Zero(mSize,mSize);
    M(0,0)=4;
    M(0,1)=1;
    M(0,mSize-1)=1;
    M(mSize-1,0)=1;
    M(mSize-1,mSize-2)=1;
    M(mSize-1,mSize-1)=4;

    for (int j = 1; j < mSize-1; ++j)
    {
        M(j,j-1)=1;
        M(j,j)=4;
        M(j,j+1)=1;
    }
//    cout<<"MATRIX: "<<endl<<M<<endl;
//    cout<<"LOCAL F VALUE:"<<endl<<localFValue<<endl;
    FullPivLU<MatrixXd> fullPivLU(M);
    cValue=fullPivLU.solve(localFValue);
//    cout<<cValue<<endl;
    return cValue;
}

//ArrayXXd CubicSplinePeriodic::findSplineValues (VectorXd coeff, VectorXd nPoints)
//{
//    VectorXd extnPoints(nPoints.size()+2);
//    extnPoints[0]=-1;
//    extnPoints[extnPoints.size()-1]=-2;
//    for (int j = 0; j < nPoints.size(); ++j)
//    {
//        extnPoints[j+1]=nPoints[j];
//    }
//    cout<<"extnPoints: "<<endl<<extnPoints<<endl;
//    VectorXd extnCoeff(coeff.size()+2);
//    extnCoeff[0]=coeff[coeff.size()-1];
//    extnCoeff[extnCoeff.size()-1]=coeff[0];
//    for (int k = 0; k < coeff.size(); ++k)
//    {
//        extnCoeff[k+1]=coeff[k];
//    }
//    cout<<"extnCoeff: "<<endl<<extnCoeff<<endl;
//    const int n=nPoints.size();
//    int uSpacePointsSize=20*n+1;
//    ArrayXXd uniformValues(uSpacePointsSize,4);
//    uniformValues.col(0)=ArrayXd::LinSpaced(uSpacePointsSize,0,1);
//    uniformValues.col(1)=6*M_PI*uniformValues.col(0);
//    uniformValues.col(2)=uniformValues.col(1).cos();
//    uniformValues.col(3)=ArrayXd::Zero(uSpacePointsSize);
//    for (int i = 0; i < uSpacePointsSize; ++i)
//    {
////        cout<<uniformValues(i,0)<<endl;
//        VectorXd vXdBi=computeBi(extnPoints,uniformValues(i,0));
////        cout<<"vXdBi size: "<<vXdBi.size()<<endl;
////        cout<<"coeff size: "<<coeff.size()<<endl;
//        double q3;
////        if(i==0)
////        {
////            q3 = vXdBi.dot(minusOneCoeff);
////        }
////        else
////        {
////            q3 = vXdBi.dot(coeff);
////        }
////        cout<<"q3: "<<q3<<endl;
//        q3 = vXdBi.dot(extnCoeff);
//        uniformValues(i,3)=q3;
//    }
//
//    return uniformValues;
//}
//

ArrayXXd CubicSplinePeriodic::findSplineValues (VectorXd coeff, VectorXd nPoints)
{
    const int n=nPoints.size();
    int uSpacePointsSize=20*(n-1)+1;
    ArrayXXd uniformValues(uSpacePointsSize,4);
    uniformValues.col(0)=ArrayXd::LinSpaced(uSpacePointsSize,0,1);
    uniformValues.col(1)=6*M_PI*uniformValues.col(0);
    uniformValues.col(2)=uniformValues.col(1).cos();
    uniformValues.col(3)=ArrayXd::Zero(uSpacePointsSize);
    VectorXd extnCoeff=VectorXd::Zero(coeff.size()+1);
    extnCoeff.head(coeff.size())=coeff;
    extnCoeff[extnCoeff.size()-1]=coeff[0];
    //cout<<"extnCoeff: "<<endl<<extnCoeff<<endl;
    for (int i = 0; i < uSpacePointsSize; ++i)
    {
        double q3;
//        if(i < 1)
//        {
//            VectorXd vXdBi=computeBi(nPoints,uniformValues(i,0));
//            double sum=0;
//            sum=vXdBi[1]*extnCoeff[n-1];
//            q3=sum+vXdBi.dot(extnCoeff);
//        }
//        else if (i > 160)
//        {
//            VectorXd vXdBi=computeBi(nPoints,uniformValues(i,0));
//            double sum=0;
//            sum=vXdBi[n-2]*coeff[n-1];
//            q3=sum+vXdBi.dot(extnCoeff);;
//        }
//        else
//        {
            VectorXd vXdBi=computeBi(nPoints,uniformValues(i,0));
//        cout<<"vXdBi size: "<<vXdBi.size()<<endl;
//        cout<<"coeff size: "<<extnCoeff.size()<<endl;

            q3 = vXdBi.dot(extnCoeff);
//        }
//        cout<<uniformValues(i,0)<<endl;

//        cout<<"vXdBi size: "<<vXdBi.size()<<endl;
//        cout<<"coeff size: "<<coeff.size()<<endl;

//        if(i==0)
//        {
//            q3 = vXdBi.dot(minusOneCoeff);
//        }
//        else
//        {
//            q3 = vXdBi.dot(coeff);
//        }
//        cout<<"q3: "<<q3<<endl;

        uniformValues(i,3)=q3;
    }

    return uniformValues;
}

VectorXd CubicSplinePeriodic::computeBi (VectorXd nPoints, double x)
{
    int n=nPoints.size();
    double a=nPoints[0],b=nPoints[n-1];
    double h=(nPoints[n-1]-nPoints[0])/(double)(n-1);
    //cout<<"n = "<<n<<", h = " <<h<<", a ="<<a<<", b ="<< b <<", x:"<<x<<endl<<"nPoints: "<<nPoints<<endl;
    VectorXd Bi=VectorXd::Zero(n);
    for (int i = 0; i < n; ++i)
    {
        double tempB;
        double xi;
        xi = nPoints[i];
        tempB=(x-xi)/h;
        Bi[i]=getBx(tempB);
        //cout<<"B("<<i<<", "<<x<<"): "<< tempB<<" => Bx("<<i<<", "<<x<<"): "<< Bi[i]<<endl;
    }
    //cout<<"Bi :"<<endl<<Bi<<endl;
    return Bi;
}


