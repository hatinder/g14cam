#include <iostream>
#include <Dense>
#include <cmath>
#include <vector>
#include "CubicSpline.hpp"
#include "CubicSplinePeriodic.hpp"
#include "RootFinding.hpp"

using namespace std;
using namespace Eigen;

void runProblem1c();
void runProblem1E ();
void runProblem1F();
void runProblem2b();
void runProblem2c();
void runProblem2D();
void runProblem2E();

int main ()
{
//    runProblem1c();
//    runProblem1E();
//    runProblem1F();
//    VectorXd test=VectorXd::LinSpaced(8,0.5,4);
//    cout<<test<<endl;
//    runProblem2b();
//    runProblem2c();
//    runProblem2D();
    runProblem2E();
    return 0;
}

void runProblem1c()
{

    cout << "Problem 1 a to c" << std::endl;
    vector<string> uValuesNames={"uSpacePoints", "6PIx", "COS(6PIx)", "q3n"};
    vector<string> aeColNames={"n", "ApproxError", "log10AE", "log10n" };
    ArrayXi n(5);
    n <<   8, 16, 32, 64, 128 ;
    ArrayXXd infinityNorm(5 ,4);
    for (int j=0 ; j < n.size(); j++)
    {
        const double a=0,b=1;
        VectorXd nodalPoints=VectorXd::LinSpaced(n[j]+1,a,b);
        VectorXd fValue(n[j]+1);
        for (int i = 0; i < n[j]+1; ++i)
        {
            fValue[i]=cos(6*M_PI*nodalPoints[i]);
        }
//        cout<<"Function Value"<<endl;
//        cout<<fValue<<endl;
//        cout<<"Function Value"<<endl;
        CubicSpline cs;
        VectorXd coefficents=cs.findCoefficients(fValue);
//        cout<<coefficents<<endl;
        ArrayXXd uValues=cs.findSplineValues(coefficents,nodalPoints);
        VectorXd ApproxError=uValues.col(2)-uValues.col(3);
        infinityNorm(j,0)=n[j];
        infinityNorm(j,1)=ApproxError.lpNorm<Infinity>();
        infinityNorm(j,2)=log10(ApproxError.lpNorm<Infinity>());
        infinityNorm(j,3)=log10(n[j]);
        //cout<<"Approx Error: "<< ApproxError.lpNorm<Infinity>() <<endl;
        cs.writeToFile("1C.", uValues, n[j], uValuesNames);
    }
    CubicSpline t;
    t.writeToFile("1C.APPROXERROR.", infinityNorm, 0, aeColNames);

}

void runProblem1E ()
{

    cout << "Problem 1 d to e" << std::endl;
    vector<string> uValuesNames={"uSpacePoints", "6PIx", "COS(6PIx)", "q3n"};
    vector<string> aeColNames={"n", "ApproxError", "log10AE", "log10n" };
    ArrayXi n(5);
    n <<  8, 16, 32, 64, 128 ;
    ArrayXXd infinityNorm(5,4);
    for (int j=0 ; j < n.size(); j++)
    {
        const double a=0,b=1;
        const double h=(b-a)/(double)(n[j]);
        //cout<<"a: "<<a<<" ,b: "<<b<<" ,b-h: "<<b-h<<endl;
        VectorXd nodalPoints=VectorXd::LinSpaced(n[j],a,b-h);
        ////cout<<"NODALPOINTS: " <<endl<< nodalPoints<<endl;
        VectorXd fValue(n[j]);
        for (int i = 0; i < n[j]; ++i)
        {
            fValue[i]=cos(6*M_PI*nodalPoints[i]);
        }
        //cout<<"FVALUE: "<<endl<<fValue<<endl;
        CubicSplinePeriodic csp;
        VectorXd coefficents=csp.findCoefficients(fValue);
        //cout<<"COEFFICIENTS: "<<endl<<coefficents<<endl;
        VectorXd splineNodalPoints=VectorXd::LinSpaced(n[j]+1,a,b);
        //cout<<"SplineNodalPoints:"<<endl<<splineNodalPoints<<endl;
        ArrayXXd uValues=csp.findSplineValues(coefficents,splineNodalPoints);
        VectorXd ApproxError=uValues.col(2)-uValues.col(3);
        infinityNorm(j,0)=n[j];
        infinityNorm(j,1)=ApproxError.lpNorm<Infinity>();
        infinityNorm(j,2)=log10(ApproxError.lpNorm<Infinity>());
        infinityNorm(j,3)=log10(n[j]);
        ////cout<<"Approx Error: "<< ApproxError.lpNorm<Infinity>() <<endl;
        csp.writeToFile("1E.", uValues, n[j], uValuesNames);
    }
    CubicSpline t;
    t.writeToFile("1E.APPROXERROR.", infinityNorm, 0, aeColNames);
}

void runProblem1F()
{
    //Right Wing
    CubicSpline cs1F ;
    ArrayXXd sn1 =cs1F.readSpline("spline_natural1.dat",8);
    ArrayXXd sn2 =cs1F.readSpline("spline_natural2.dat",8);
    ArrayXXd sp1 =cs1F.readSpline("spline_periodic.dat",6);
    vector<string> xyValuesNames={"x", "y"};
    ArrayXXd xyValue=ArrayXXd::Zero(20*(sn1.rows()-1)+1,2);
    VectorXd nodalPoints=VectorXd::LinSpaced(sn1.rows(),0,1);
    VectorXd xValue=sn1.col(0);
    VectorXd yValue=sn1.col(1);
    CubicSpline cs;
    VectorXd coefficents=cs.findCoefficients(xValue);
    ArrayXXd xValues=cs.findSplineValues(coefficents,nodalPoints);
    cout<<"xValues.rows(): "<<xValues.rows()<<endl;
    cout<<"xyValue.rows(): "<<xyValue.rows()<<endl;
    for (int i = 0; i < xValues.rows(); ++i)
    {
        xyValue(i,0)=xValues(i,3);
    }
    coefficents=cs.findCoefficients(yValue);
    ArrayXXd yValues=cs.findSplineValues(coefficents,nodalPoints);
    for (int i = 0; i < yValues.rows(); ++i)
    {
        xyValue(i,1)=yValues(i,3);
    }
    ostringstream iterate_label1;
    iterate_label1.width(3);
    iterate_label1.fill('0');
    iterate_label1 << 1;
    string file_name = "1FSN1." + iterate_label1.str();
    cs.writeToFile(file_name, xyValue, 1, xyValuesNames);

    //Left Wing
    xValue=sn2.col(0);
    yValue=sn2.col(1);
    coefficents=cs.findCoefficients(xValue);
    xValues=cs.findSplineValues(coefficents,nodalPoints);
    for (int i = 0; i < xValues.rows(); ++i)
    {
        xyValue(i,0)=xValues(i,3);
    }
    coefficents=cs.findCoefficients(yValue);
    yValues=cs.findSplineValues(coefficents,nodalPoints);
    for (int i = 0; i < yValues.rows(); ++i)
    {
        xyValue(i,1)=yValues(i,3);
    }
    ostringstream iterate_label2;
    iterate_label2.width(3);
    iterate_label2.fill('0');
    iterate_label2 << 2;
    file_name = "1FSN1." + iterate_label2.str();
    cs.writeToFile(file_name, xyValue, 2, xyValuesNames);

    //Centre
    ArrayXXd xyPValue=ArrayXXd::Zero(20*sp1.rows()+1,2);
    VectorXd pNodalPoints=VectorXd::LinSpaced(sp1.rows()+1,0,1);
    VectorXd xpValue=sp1.col(0);
    VectorXd ypValue=sp1.col(1);
    CubicSplinePeriodic csp;
    coefficents=csp.findCoefficients(xpValue);
    xValues=csp.findSplineValues(coefficents,pNodalPoints);
    cout<<"xValues.rows(): "<<xValues.rows()<<endl;
    cout<<"xyPValue.rows(): "<<xyPValue.rows()<<endl;
    for (int i = 0; i < xValues.rows(); ++i)
    {
        xyPValue(i,0)=xValues(i,3);
    }
    coefficents=csp.findCoefficients(ypValue);
    yValues=csp.findSplineValues(coefficents,pNodalPoints);
    for (int i = 0; i < yValues.rows(); ++i)
    {
        xyPValue(i,1)=yValues(i,3);
    }
    ostringstream iterate_label3;
    iterate_label3.width(3);
    iterate_label3.fill('0');
    iterate_label3 << 3;
    file_name = "1FSN1." + iterate_label3.str();
    cs.writeToFile(file_name, xyPValue, 3, xyValuesNames);

}

void runProblem2b()
{
    RootFinding rootFinding;
    ArrayXXd initialGuess=ArrayXXd::Zero(8,3);
    ArrayXXd roots=ArrayXXd::Zero(4,3);
    initialGuess<<0.1,0.1,0.1,
            -0.1,0.1,0.1,
            0.1,-0.1,0.1,
            0.1,0.1,-0.1,
            -0.1,-0.1,0.1,
            0.1,-0.1,-0.1,
            -0.1,0.1,-0.1,
            -0.1,-0.1,-0.1;
//    cout<<rootFinding.calculateF(0.1,0.1,0.1)<<endl;
//    cout<<rootFinding.calculateJ(0.1,0.1,0.1)<<endl;
    int nRoots=0;
    for (int i = 0; i < initialGuess.rows(); ++i)
    {
        VectorXd xOld=initialGuess.row(i);
        VectorXd xNew=VectorXd::Zero(3);
        double TOL=pow(10,-12);
        double tolerance=1;
        int j=0;
        for (; j < 100 && tolerance > TOL; ++j)
        {
//            cout<<"xOld: "<< xOld<<endl;
            VectorXd F=rootFinding.calculateF(xOld(0),xOld(1),xOld(2));
            MatrixXd J=rootFinding.calculateJ(xOld(0),xOld(1),xOld(2));
            FullPivLU<MatrixXd> JInv(J);
            VectorXd y=-JInv.solve(F);
            xNew=xOld+y;
//            cout<<"xNew: "<< xNew<<endl;
            tolerance=(xNew-xOld).lpNorm<2>();
            xOld=xNew;
        }
        if(rootFinding.foundNewRoot(xNew(0),xNew(1),xNew(2),roots))
        {
            cout<<"Initial Guess: ("<<initialGuess(i,0)<<", "<<initialGuess(i,1)<<", "<<initialGuess(i,2)<<")"<<endl;
            cout<<"Iteration: "<<j<<", New Root Found: ("<<xNew(0)<<", "<<xNew(1)<<", "<<xNew(2)<<")"<<endl;
            roots.row(nRoots)=xNew;
            nRoots++;
            if(nRoots==4)
                break;
        }
        else
        {
//            cout<<"Initial Guess: ("<<initialGuess(i,0)<<", "<<initialGuess(i,1)<<", "<<initialGuess(i,2)<<")"<<endl;
//            cout<<"Iteration: "<<j<<", Not New Root : ("<<xNew(0)<<", "<<xNew(1)<<", "<<xNew(2)<<")"<<endl;
        }
    }
    cout<<"All Roots: "<<endl<<roots<<endl;
}

void runProblem2c()
{
    RootFinding rootFinding;
    vector<string> rootColNames={"x", "y"};
    ArrayXXd roots=ArrayXXd::Zero(200,2);
//    cout<<rootFinding.calculateF(0.1,0.1,0.1)<<endl;
//    cout<<rootFinding.calculateJ(0.1,0.1,0.1)<<endl;
    int nRoots=0;
    double ds=0.05;
    VectorXd xInitial(2);
    VectorXd unitTangent=VectorXd::Zero(2);
    xInitial(0)=0;xInitial(1)=2;
    cout<<"Intial Guess: "<<endl<<xInitial<<endl;
    cout<<"xInitial.normalized(): "<<endl<<xInitial.normalized()<<endl;
    VectorXd xNew=VectorXd::Zero(2);
    VectorXd xNew2=VectorXd::Zero(2);
    while (nRoots < roots.rows())
    {
//        unitTangent=ds*xOld.normalized();
//        cout<<"unitTangent: "<<endl<<unitTangent<<endl;
        //double TOL=0.005;
        unitTangent(0)=xInitial(1)/xInitial.lpNorm<2>();
        unitTangent(1)=-xInitial(0)/xInitial.lpNorm<2>();
        xNew=xInitial+ds*unitTangent;
//        xInitial=xInitial.normalized();
        double TOL=pow(10,-2);
        double tolerance=1;
        int j=0;
        for (; j < 100 && tolerance > ds; ++j)
        {
            cout<<"xNew: "<<endl<<xNew<<endl;
            cout<<"xInitial: "<<endl<<xInitial<<endl;
            VectorXd F=rootFinding.calculateF2C(xNew(0),xNew(1),xInitial(0),xInitial(1),xInitial.lpNorm<2>());
            cout<<"F(x): "<<endl<<F<<endl;
            MatrixXd J=rootFinding.calculateJ2C(xNew(0),xNew(1),xInitial(0),xInitial(1));
            cout<<"Matrix J: "<<endl<<J<<endl;
            FullPivLU<MatrixXd> JInv(J);
            VectorXd y=-JInv.solve(F);
            cout<<"y: "<<endl<<y<<endl;
            xNew2=xNew+y;
            cout<<"xNew2: "<<endl<<xNew2<<endl;
            tolerance=(xNew2-xNew).lpNorm<2>();
            cout<<"tolerance: "<< tolerance<<endl;
            //tolerance=(xNew2-xNew).squaredNorm();
            xNew=xNew2;

        }
        if(rootFinding.foundNewRoot(xNew2(0),xNew2(1),roots))
        {
//            cout<<"xOld: ("<<xOld(0)<<", "<<xOld(1)<<")"<<endl<<"xNew: ("<<xOld(0)<<", "<<xOld(1)<<")"<<endl;
            cout<<"Iteration: "<<j<<", New Root Found: ("<<xNew2(0)<<", "<<xNew2(1)<<")"<<endl;
            roots.row(nRoots)=xNew2;
            nRoots++;
        }
        else
        {
//            cout<<"Initial Guess: ("<<initialGuess(i,0)<<", "<<initialGuess(i,1)<<", "<<initialGuess(i,2)<<")"<<endl;
//            cout<<"Iteration: "<<j<<", Not New Root : ("<<xNew(0)<<", "<<xNew(1)<<", "<<xNew(2)<<")"<<endl;
        }
        xInitial=xNew2;
    }
    cout<<"All Roots: "<<endl<<roots<<endl;
    rootFinding.writeToFile("2C",roots,0,rootColNames);
}

void runProblem2D()
{
    cout << "Problem 2D" << std::endl;
    vector<string> rootColNames={"h", "u"};
    int n=64;
    //double lambda=0.5;
    VectorXd lambda=VectorXd::LinSpaced(8,0.5,4);
    double h = (double)1/n;
    RootFinding rootFinding;
    for (int j = 0; j < lambda.size(); ++j)
    {
        ArrayXXd roots=ArrayXXd::Zero(n+1,2);
        roots.col(0)=ArrayXd::LinSpaced(n+1,0,1);
        VectorXd uInitial=VectorXd::Zero(n+1);
        double TOL=pow(10,-12);
        double tolerance=1;
        int i=0;
        for (;i < 100 && tolerance > TOL ; ++i)
        {
    //        cout<<"uInitial: "<<endl<<uInitial<<endl;
            VectorXd uNew=VectorXd::Zero(n+1);
            VectorXd y=VectorXd::Zero(n-1);
            VectorXd F_U=rootFinding.calculateF2D(uInitial,h,lambda(j));
    //        cout<<"F_U: "<<endl<<F_U<<endl;
            MatrixXd J_U=rootFinding.calculateJ2D(uInitial,h,lambda(j));
    //        cout<<"J_U: "<<endl<<J_U<<endl;
            PartialPivLU<MatrixXd> JInv(J_U);
            y=-JInv.solve(F_U);
            uNew.segment(1,n-1)=uInitial.segment(1,n-1)+y.segment(0,n-1);
            tolerance=(uNew-uInitial).lpNorm<2>();
    //        cout<<"tolerance: "<< tolerance<<endl;
            uInitial=uNew;
        }
        cout<<"Iteration: "<<i<< " for lambda: "<<lambda(j)<<endl;
        roots.col(1)=uInitial;
    //    cout<<"All Roots: "<<endl<<roots<<endl;
        rootFinding.writeToFile("2D",roots,lambda(j)*10,rootColNames);
    }

}

void runProblem2E()
{
    cout << "Problem 2E" << std::endl;
    vector<string> rootColNames={"h", "u"};
    int n=64;
    int nRoots=0;
    ArrayXXd roots=ArrayXXd::Zero(n+1,201);
    roots.col(0)=ArrayXd::LinSpaced(n+1,0,1);
    ArrayXXd otherInfo=ArrayXXd::Zero(4,201);
    double ds=0.125,TOL=0.001;
    RootFinding rootFinding;
    double lambdaInitial=0.5;
    double h = (double)1/n;
    double tolerance=1,lbdTilda,lbdNew,lbdyNew;
    VectorXd uInitial= rootFinding.findInitialU(lambdaInitial, n);
    VectorXd tk = VectorXd::Zero(n-1);
    VectorXd uNew = VectorXd::Zero(n + 1);
    VectorXd uyNew = VectorXd::Zero(n + 1);
    VectorXd lastRow;
    while (nRoots < roots.cols()-1)
    {
//        cout << "uInitial: " << endl << uInitial << endl;
        MatrixXd FuFlambda = rootFinding.createFuLambda2E(uInitial, h, lambdaInitial);
//        cout << "FuFlambda: " << endl << FuFlambda << endl;
        JacobiSVD<MatrixXd> svd(FuFlambda, ComputeFullU | ComputeFullV);
        MatrixXd Vstar = svd.matrixV();
        lastRow = Vstar.col(Vstar.cols() - 1);
//        cout << "lastRow:" << endl << lastRow << endl;
//        lastRow=lastRow/lastRow.lpNorm<2>(); Already a unit vector!! //TODO: Remove it
//        cout << "lastRow Unit Vector:" << endl << lastRow << endl;
        tk = lastRow.segment(0, lastRow.size() - 1);
        lbdTilda = lastRow.tail(1)(0);
//      cout<<"svd.matrixV(): "<<endl<<svd.matrixV()<<endl;
//      cout<<"svd.singularValues(): "<<endl<<svd.singularValues()<<endl;
//      cout<<"svd.matrixU(): "<<endl<<svd.matrixU()<<endl;
//      cout<<"Vstar: "<<endl<<Vstar<<endl;
//        cout<<"tk:"<<endl<<tk<<endl;
//        cout<<"lbdTilda:"<<endl<<lbdTilda<<endl;
        uNew.segment(1, n - 1) = uInitial.segment(1, n - 1) + ds * (tk);
        lbdNew = lambdaInitial + ds * lbdTilda;
//        cout<<"uNew: "<<endl<<uNew<<endl;
//        cout<<"lbdNew: "<<endl<<lbdNew<<endl;
        int j=0;
        tolerance=1;
        for (; j < 100 && tolerance > ds; ++j)
        {
            VectorXd F = rootFinding.calculateF2E(uNew, h, lbdNew, uInitial, tk, lbdTilda, lambdaInitial,ds);
//            cout << "F: " << endl << F << endl;
            MatrixXd J = rootFinding.calculateJ2E(uNew, h, lbdNew, tk, lbdTilda);
//            cout << "J: " << endl << J << endl;
            FullPivLU<MatrixXd> JInv(J);
            VectorXd y = -JInv.solve(F);
//            cout << "y: " << endl << y << endl;
//            cout<< "uInitial.segment(1, n - 2): "<<endl<<uInitial.segment(1, n - 1)<<endl;
//            cout<< " y.segment(0, n - 2): "<<endl<<y.segment(0, n - 1)<<endl;
            uyNew.segment(1, n - 1) = uInitial.segment(1, n - 1) + y.segment(0, n - 1);
//            cout<< "uyNew: "<<endl<<uyNew<<endl;
            lbdyNew = y.tail(1)(0);
            VectorXd tolCalcNew(n + 2), tolCalcOld(n + 2);
            tolCalcNew.segment(0,uyNew.size()) = uyNew;
            tolCalcNew[uyNew.size()] = lbdyNew;
//            cout<<"tolCalcNew: "<<endl<<tolCalcNew<<endl;
            tolCalcOld.segment(0,uNew.size()) = uNew;
            tolCalcOld[uyNew.size()] = lbdNew;
//            cout<<"tolCalcOld: "<<endl<<tolCalcOld<<endl;
            tolerance = (tolCalcNew - tolCalcOld).lpNorm<2>();
//            cout << "tolerance: " << tolerance << endl;
            uNew=uyNew;
            lbdNew=lbdyNew;
        }
        for (int i = 0; i < uyNew.size(); ++i)
        {
            roots(i,nRoots+1)=uyNew[i];
        }
        otherInfo(0,nRoots)=nRoots;
        otherInfo(1,nRoots)=lbdyNew;
        otherInfo(1,nRoots)=tolerance;
        otherInfo(1,nRoots)=j;
        nRoots++;
        uInitial=uNew;
        lambdaInitial=lbdNew;
        cout<<"nRoots: "<<nRoots<<endl;
    }
//    cout<<"Roots: "<<endl<<roots<<endl;
    rootFinding.writeToFile("2EROOTS",roots,0,rootColNames);
    rootFinding.writeToFile("2EOTHERINFO",otherInfo,0,rootColNames);
}