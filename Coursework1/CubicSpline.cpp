//
// Created by hsingh9 on 08/03/2019.
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "CubicSpline.hpp"
using namespace std;
using namespace Eigen;


VectorXd CubicSpline::findCoefficients (VectorXd fValue)
{
    const int mSize=fValue.size() + 2;
    VectorXd cValue=VectorXd::Zero(mSize);
    VectorXd localFValue=VectorXd::Zero(mSize);
    localFValue(0)=0;
    localFValue(mSize-1)=0;
    for (int i = 1; i < mSize - 1; ++i)
    {
        localFValue[i]=fValue[i-1];
    }
    //cout<<localFValue<<endl;
    MatrixXd M=MatrixXd::Zero(mSize,mSize);
    M(0,0)=1;
    M(0,1)=-2;
    M(0,2)=1;
    M(mSize-1,mSize-3)=1;
    M(mSize-1,mSize-2)=-2;
    M(mSize-1,mSize-1)=1;

    for (int j = 1; j < mSize-1; ++j)
    {
        M(j,j-1)=1;
        M(j,j)=4;
        M(j,j+1)=1;
    }
    //cout<<M<<endl;
    //cout<<fValue<<endl;
    PartialPivLU<MatrixXd> partialPivLU(M);
    cValue=partialPivLU.solve(localFValue);
    //cout<<cValue<<endl;
    return cValue;
}

ArrayXXd CubicSpline::findSplineValues (VectorXd coeff, VectorXd nPoints)
{
    const int n=nPoints.size() - 1;
    double a=nPoints[0],b=nPoints[n];
    int uSpacePointsSize=20*n+1;
    ArrayXXd uniformValues(uSpacePointsSize,4);
    uniformValues.col(0)=ArrayXd::LinSpaced(uSpacePointsSize,a,b);
    uniformValues.col(1)=6*M_PI*uniformValues.col(0);
    uniformValues.col(2)=uniformValues.col(1).cos();
    uniformValues.col(3)=ArrayXd::Zero(uSpacePointsSize);
    for (int i = 0; i < uSpacePointsSize; ++i)
    {
//        cout<<uniformValues(i,0)<<endl;
        VectorXd vXdBi=computeBi(nPoints,uniformValues(i,0));
//        cout<<"vXdBi size: "<<vXdBi.size()<<endl;
//        cout<<"coeff size: "<<coeff.size()<<endl;
        double q3=vXdBi.dot(coeff);
//        cout<<"q3: "<<q3<<endl;
        uniformValues(i,3)=q3;
    }

    return uniformValues;
}

double CubicSpline::getBx (double xVal)
{
    double Bx;
    if(xVal >= -2 && xVal < -1)
    {
        Bx=pow((xVal + 2), 3);
    }
    else if(xVal >= -1 && xVal < 0)
    {
        Bx=(1.0 + 3 * (xVal + 1) + 3 * (pow((xVal + 1), 2)) - 3 * (pow((xVal + 1), 3)));
    }
    else if(xVal >= 0 && xVal < 1)
    {
        Bx=(1.0 + 3 * (1 - xVal) + 3 * (pow((1 - xVal), 2)) - 3 * (pow((1 - xVal), 3)));
    }
    else if(xVal >= 1 && xVal < 2)
    {
        Bx=pow((2 - xVal), 3);
    }
    else
    {
        Bx=0;
    }
    return Bx;
}

VectorXd CubicSpline::computeBi (VectorXd nPoints, double x)
{
    int n=nPoints.size();
    double a=nPoints[0],b=nPoints[n-1];
    double h=(b-a)/(double)(n-1);
//    cout<<"n = "<<n<<", h = "<<h<<", a ="<<a<<", b ="<< b <<endl;
    VectorXd Bi=VectorXd::Zero(n+2);
    for (int i = 0; i < n+2; ++i)
    {
        double tempB;
        if(i==0)
        {
            tempB=(x-(a-h))/h;
        }
        else if(i == n+1)
        {
            tempB=(x-(b+h))/h;
        }
        else
        {
            tempB=(x-nPoints[i-1])/h;
        }
        Bi[i]=getBx(tempB);
//        cout<<"B("<<i<<", "<<x<<"): "<< tempB<<" => Bx("<<i<<", "<<x<<"): "<< Bi[i]<<endl;
    }
    return Bi;
}

void CubicSpline::writeToFile (const string fNamePrefix, ArrayXXd uniEvalPoints, const int k, vector<string> colNames)
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
    oFileStream << uniEvalPoints <<endl;
    oFileStream.close();
}

ArrayXXd CubicSpline::readSpline (const string &filename, int size)
{
    ArrayXXd splineData=ArrayXXd::Zero(size,2);
    ifstream ifs;
    ifs.open(filename.c_str(), ios_base::in);
    if (!ifs)
    {
        cout << "Cant Open input file: " << filename << endl;
        exit(1);
    }
    string s;
    int i=0;
    while(getline(ifs,s))
    {
        istringstream iss(s);
        double a,b;
        iss >>a>>b;
        splineData.row(i)<<a,b;
        i++;
    }
    return splineData;
}

