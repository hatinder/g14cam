//
// Created by hsingh9 on 20/04/2019.
//

#include "Problem2.hpp"
#include "helper.hpp"
#include "ODE.hpp"
#include "Utility.hpp"
#include <vector>
#include <cmath>
//global function F for Problem 2a
auto gFA = [] ( double t, vector<double> x)
{
    vector<double> results(2, 0.0);
    results[0]=x[1];
    results[1]=1.0+x[1]*sin(t)-x[0]*(x[0]+1.0);
    return results;
};

//global function F for Problem 2b
auto gFB = [] (double theta1, double theta2,double p1, double p2)
{
    vector<double> results(4, 0.0);
    results[0] = (p1 - p2 * cos(theta1 - theta2)) / (1.0 + pow(sin(theta1 - theta2),2));
    results[1] = (2.0 * p2 - p1 * cos(theta1 - theta2)) / (1.0 + pow(sin(theta1 - theta2),2));
    double a = (2.0 * p1 * p2 * sin(theta1 - theta2) * (1.0 + pow(sin(theta1 - theta2),2)));
    double b = (sin(2.0 * (theta1 - theta2)) * (p1 * p1 - 2.0 * p1 * p2 * cos(theta1 - theta2) + 2 * p2 * p2));
    double c = (2.0 * pow((1.0 + pow(sin(theta1 - theta2),2)),2));
    double g = 9.8;
    results[2] = -((a - b) / c + 2.0 * g * sin(theta1));
    results[3] = -((b - a) / c + g * sin(theta2));
    return results;
};

//global function F for Problem 2b
auto gFBN = [] (double t, vector<double> x)
{
    double theta1=x[0],theta2=x[1],p1=x[2],p2=x[3];
    vector<double> results(4, 0.0);
    results[0] = (p1 - p2 * cos(theta1 - theta2)) / (1.0 + pow(sin(theta1 - theta2),2));
    results[1] = (2.0 * p2 - p1 * cos(theta1 - theta2)) / (1.0 + pow(sin(theta1 - theta2),2));
    double a = (2.0 * p1 * p2 * sin(theta1 - theta2) * (1.0 + pow(sin(theta1 - theta2),2)));
    double b = (sin(2.0 * (theta1 - theta2)) * (p1 * p1 - 2.0 * p1 * p2 * cos(theta1 - theta2) + 2 * p2 * p2));
    double c = (2.0 * pow((1.0 + pow(sin(theta1 - theta2),2)),2));
    double g = 9.8;
    results[2] = -((a - b) / c + 2.0 * g * sin(theta1));
    results[3] = -((b - a) / c + g * sin(theta2));
    return results;
};

/*
// * Old one
void Problem2:: ARungeKutta2 ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (ARungeKutta2)" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"t", "y"};
    double u0 = 1, v0 = 0, T = 4 * M_PI;
    auto f1 = [] (double v, double u, double t)
    { return (1 + v * sin(t) - u * (u + 1)); };
    auto f2 = [] (double v)
    { return v; };
    vector<int> N = {128, 256, 512, 1024, 2048, 4096};
    ODE ode;
    Utility utility;
    map<double, double> keyNValueYPair;
    for (unsigned i = 0; i < N.size(); ++i)
    {
        vector<map<double, double>> rk2 = ode.applyRungeKutta2(f1, f2, u0, v0, T, N[i]);
        utility.writeToFile("PROBLEM2ARK2V", rk2[0], i, colNames);
        utility.writeToFile("PROBLEM2ARK2U", rk2[1], i, colNames);
        keyNValueYPair[N[i]] = (cos(T) - rk2[1].rbegin()->second);
    }
    utility.writeToFile("PROBLEM2ARK2ERROR", keyNValueYPair, 0, colNames);
}
*/

// * New one
void Problem2:: ARungeKutta2 ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (A) Runge Kutta 2" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"t", "x1", "x2"};
    vector<double> initVal={1.0,0.0};
    double T = 4 * M_PI;
    vector<int> N = {128, 256, 512, 1024, 2048, 4096};
    ODE ode;
    Utility utility;
    map<double, double> keyNValueYPair;
    for (unsigned i = 0; i < N.size(); ++i)
    {
        map<double,vector<double>> rk2 = ode.applyRK2(gFA, initVal, T, N[i]);
        utility.writeToFile("PROBLEM2ARK2",rk2, i, colNames);
        keyNValueYPair[N[i]] = (cos(T) - rk2.rbegin()->second[0]);
    }
    utility.writeToFile("PROBLEM2ARK2ERROR", keyNValueYPair, 0, colNames);
}

/*
 * old one
void Problem2::BRungeKutta2 ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (B) Runge Kutta 2" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"a", "b"};
    vector<double> iv1={M_PI/6.0,M_PI/6.0,0.0,0.0};         //Initial Values
    vector<double> iv2={(3.0*M_PI)/4.0,(3.0*M_PI)/4.0,0.0,0.0}; //Initial Value s
    double T=100, N=5000;
    ODE ode;
    vector<map<double,double>> keyValuePair=ode.applyRungeKutta2(gFB, iv1, T, N);
    Utility utility;
    utility.writeToFile("PROBLEM2BRK2THETA1",keyValuePair[0],0,colNames);
    utility.writeToFile("PROBLEM2BRK2THETA2",keyValuePair[1],0,colNames);
    map<double,double> THPair1=computeHamilton(keyValuePair,T,N);
    utility.writeToFile("PRB2BRK2HMT",THPair1,1,colNames);
//    utility.writeToFile("PROBLEM2BRK2P",keyValuePair[1],0,colNames);
    vector<map<double,double>> keyValuePair2=ode.applyRungeKutta2(gFB, iv2, T, N);
    utility.writeToFile("PROBLEM2B2RK2THETA1",keyValuePair2[0],0,colNames);
    utility.writeToFile("PROBLEM2B2RK2THETA2",keyValuePair2[1],0,colNames);
    map<double,double> THPair2=computeHamilton(keyValuePair2,T,N);
    utility.writeToFile("PRB2BRK2HMT",THPair2,2,colNames);

}
*/


void Problem2::BRungeKutta2 ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (B) Runge Kutta 2" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"THETA1", "THETA2","P1","P2"};
    vector<double> iv1={M_PI/6.0,M_PI/6.0,0.0,0.0};         //Initial Values
    vector<double> iv2={(3.0*M_PI)/4.0,(3.0*M_PI)/4.0,0.0,0.0}; //Initial Value s
    double T=100, N=5000;
    ODE ode;
    map<double,vector<double>> keyValuePair=ode.applyRK2(gFBN, iv1, T, N);
    Utility utility;
    utility.writeToFile("PROBLEM2BRK2",keyValuePair, 1, colNames);
    map<double,double> THPair1=computeHamilton(keyValuePair);
    utility.writeToFile("PRB2BRK2HMT",THPair1,1,colNames);
    map<double,vector<double>> keyValuePair2=ode.applyRK2(gFBN, iv2, T, N);
    utility.writeToFile("PROBLEM2BRK2",keyValuePair2, 2, colNames);
    map<double,double> THPair2=computeHamilton(keyValuePair2);
    utility.writeToFile("PRB2BRK2HMT",THPair2,2,colNames);

}


void Problem2::CImplicitMidpoint ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (C) Implicit Midpoint" << endl;
    cout << "=====================" << endl;

    vector<string> colNames = {"a", "b"};
    map<double,double> theta2TcPair;
    int m=100;
    for (int i = 0; i < m; ++i)
    {
        vector<double> iv1={(3.0*M_PI)/6.0,(M_PI/3.0)+(double)i*M_PI/200.0,0.0,0.0};         //Initial Values
        double T=100, N=5000;
        ODE ode;
        map<double,vector<double>> keyValuePair=ode.applyIMP(gFBN, iv1, T, N);
        double pi=M_PI;
        double time=findByValue(keyValuePair,pi);
//        cout<<"Theta2: "<<iv1[1]<<" , time: "<<time<<endl;
        theta2TcPair[iv1[1]]=time;   //TODO: Pending Review overloaded function only accepts same data type
    }
    Utility utility;
    utility.writeToFile("PRB2CTT",theta2TcPair,0,colNames);
}

/*
 * Old One
void Problem2::AImplicitMidpoint ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (A) Implicit Midpoint" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"t", "y"};
    double u0 = 1, v0 = 0, T = 4 * M_PI;
    auto f1 = [] (double v, double u, double t, double vi, double ui)
    { return (1 + (v + vi) / (2.0) * sin(t) - (u + ui) / (2.0) * ((u + ui) / (2.0) + 1)); };
    auto f2 = [] (double v, double vi)
    { return (v + vi) / 2.0; };
    vector<int> N = {128, 256, 512, 1024, 2048, 4096};
    ODE ode;
    Utility utility;
    map<double, double> keyNValueYPair;
    for (unsigned i = 0; i < N.size(); ++i)
    {
        vector<map<double, double>> rk2 = ode.applyImplicitMidpoint(f1, f2, u0, v0, T, N[i]);
        utility.writeToFile("PROBLEM2AIMV", rk2[0], i, colNames);
        utility.writeToFile("PROBLEM2AIMU", rk2[1], i, colNames);
        keyNValueYPair[N[i]] = (cos(T) - rk2[1].rbegin()->second);
    }
    utility.writeToFile("PROBLEM2AIMERROR", keyNValueYPair, 0, colNames);

}
*/

void Problem2::AImplicitMidpoint ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (A) Implicit Midpoint" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"t", "x1", "x2"};
    vector<double> initVal={1.0,0.0};
    double T = 4 * M_PI;
    vector<int> N = {128, 256, 512, 1024, 2048, 4096};
    ODE ode;
    Utility utility;
    map<double, double> keyNValueYPair;
    for (unsigned i = 0; i < N.size(); ++i)
    {
        map<double,vector<double>> rk2 = ode.applyIMP(gFA,initVal,T,N[i]);
        utility.writeToFile("PROBLEM2AIMP",rk2, i, colNames);
        keyNValueYPair[N[i]] = (cos(T) - rk2.rbegin()->second[0]);
    }
    utility.writeToFile("PROBLEM2AIMERROR", keyNValueYPair, 0, colNames);

}

void Problem2::BImplicitMidpoint ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (B) Implicit Midpoint" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"a", "b"};
    vector<double> iv1={M_PI/6.0,M_PI/6.0,0.0,0.0};         //Initial Values
    vector<double> iv2={(3.0*M_PI)/4.0,(3.0*M_PI)/4.0,0.0,0.0}; //Initial Values
    double T=100, N=5000;
    ODE ode;
    map<double,vector<double>> keyValuePair=ode.applyIMP(gFBN, iv1, T, N);
    Utility utility;
    utility.writeToFile("PROBLEM2BIMP",keyValuePair, 1, colNames);
    map<double,double> THPair1=computeHamilton(keyValuePair);
    utility.writeToFile("PRB2BIMHMT",THPair1,1,colNames);
//    utility.writeToFile("PROBLEM2BRK2P",keyValuePair[1],0,colNames);
    map<double,vector<double>> keyValuePair2=ode.applyIMP(gFBN, iv2, T, N);
    utility.writeToFile("PROBLEM2BIMP",keyValuePair2, 2, colNames);
    map<double,double> THPair2=computeHamilton(keyValuePair2);
    utility.writeToFile("PRB2BIMHMT",THPair2,2,colNames);


}

/*
 * Old one
map<double, double> Problem2::computeHamilton (vector<map<double, double>> yVector, double T, double N)
{
    map<double, double> THPair;
    double dt=T/N, t0=0;
    for (unsigned i = 0; i < yVector[0].size(); ++i)
    {
        double g=9.8;
        double p1=yVector[0][t0];
        double p2=yVector[1][t0];
        double t1=yVector[2][t0];
        double t2=yVector[3][t0];
        double H=(p1*p1+2*p2*p2-2*p1*p2*cos(t1-t2))/(2*(1+sin(t1-t2)*sin(t1-t2)))-2*g*cos(t1)-g*cos(t2);
        THPair[t0]=H;
        t0+=dt;
    }
    return THPair;
}
*/

map<double, double> Problem2::computeHamilton (map<double,vector<double>> yVector)
{
    map<double, double> THPair;
    double g=9.8;
    for(auto x:yVector)
    {
        double t1=x.second[0];
        double t2=x.second[1];
        double p1=x.second[2];
        double p2=x.second[3];
        double H=(p1*p1+2*p2*p2-2*p1*p2*cos(t1-t2))/(2*(1+sin(t1-t2)*sin(t1-t2)))-2*g*cos(t1)-g*cos(t2);
        THPair[x.first]=H;
    }
    return THPair;
}
