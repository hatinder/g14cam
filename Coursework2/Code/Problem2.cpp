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
// Please note first argument t is not used but kept it for design consistency
auto gFB = [] (double t, vector<double> x)
{
    double theta1=x[0],theta2=x[1],p1=x[2],p2=x[3];
    vector<double> results(4, 0.0);
    results[0] = (p1 - p2 * cos(theta1 - theta2)) / (1.0 + pow(sin(theta1 - theta2),2)); //theta1'
    results[1] = (2.0 * p2 - p1 * cos(theta1 - theta2)) / (1.0 + pow(sin(theta1 - theta2),2)); //theta2'
    double a = (2.0 * p1 * p2 * sin(theta1 - theta2) * (1.0 + pow(sin(theta1 - theta2),2)));
    double b = (sin(2.0 * (theta1 - theta2)) * (p1 * p1 - 2.0 * p1 * p2 * cos(theta1 - theta2) + 2 * p2 * p2));
    double c = (2.0 * pow((1.0 + pow(sin(theta1 - theta2),2)),2));
    double g = 9.8;
    results[2] = -((a - b) / c + 2.0 * g * sin(theta1));    //p1'
    results[3] = -((b - a) / c + g * sin(theta2));  //p2'
    return results;
};

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


void Problem2::BRungeKutta2 ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (B) Runge Kutta 2" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"t","THETA1", "THETA2","P1","P2"};
    vector<double> iv1={M_PI/6.0,M_PI/6.0,0.0,0.0};         //Initial Values
    vector<double> iv2={(3.0*M_PI)/4.0,(3.0*M_PI)/4.0,0.0,0.0}; //Initial Value s
    double T=100, N=5000;
    ODE ode;
    map<double,vector<double>> keyValuePair=ode.applyRK2(gFB, iv1, T, N);
    Utility utility;
    utility.writeToFile("PROBLEM2BRK2",keyValuePair, 1, colNames);
    map<double,double> THPair1=computeHamilton(keyValuePair);
    utility.writeToFile("PRB2BRK2HMT",THPair1,1,colNames);
    map<double,vector<double>> keyValuePair2=ode.applyRK2(gFB, iv2, T, N);
    utility.writeToFile("PROBLEM2BRK2",keyValuePair2, 2, colNames);
    map<double,double> THPair2=computeHamilton(keyValuePair2);
    utility.writeToFile("PRB2BRK2HMT",THPair2,2,colNames);
}


void Problem2::CImplicitMidpoint ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (C) Implicit Midpoint" << endl;
    cout << "=====================" << endl;

    vector<string> colNames = {"t","THETA1", "THETA2","P1","P2"};
    vector<string> colNames2 = {"THETA2","time"};
    map<double,double> theta2TcPair;
    Utility utility;
    int m=100;
    for (int i = 0; i < m; ++i)
    {
        vector<double> iv1={(3.0*M_PI)/4.0,(M_PI/3.0)+(double)i*M_PI/200.0,0.0,0.0};         //Initial Values
        double T=100, N=5000;
        ODE ode;
        map<double,vector<double>> keyValuePair=ode.applyIMP(gFB, iv1, T, N);
//        utility.writeToFile("PROBLEM2CIMP",keyValuePair, i, colNames);
        double pi=M_PI;
        double time=findByValue(keyValuePair,pi);
//        cout<<"Theta2: "<<iv1[1]<<" , time: "<<time<<endl;
        theta2TcPair[iv1[1]]=time;
    }
    utility.writeToFile("PRB2CTT",theta2TcPair,0,colNames2);
}


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
        map<double,vector<double>> imp = ode.applyIMP(gFA,initVal,T,N[i]);
        utility.writeToFile("PROBLEM2AIMP",imp, i, colNames);
        keyNValueYPair[N[i]] = (cos(T) - imp.rbegin()->second[0]);  //rbegin=reverse sort,second=value[0]=u
    }
    utility.writeToFile("PROBLEM2AIMERROR", keyNValueYPair, 0, colNames);
}

void Problem2::BImplicitMidpoint ()
{
    cout << "=====================" << endl;
    cout << "Running Problem 2 (B) Implicit Midpoint" << endl;
    cout << "=====================" << endl;
    vector<string> colNames = {"t","THETA1", "THETA2","P1","P2"};
    vector<double> iv1={M_PI/6.0,M_PI/6.0,0.0,0.0};         //Initial Values
    vector<double> iv2={(3.0*M_PI)/4.0,(3.0*M_PI)/4.0,0.0,0.0}; //Initial Values
    double T=100, N=5000;
    ODE ode;
    map<double,vector<double>> keyValuePair=ode.applyIMP(gFB, iv1, T, N);
    Utility utility;
    utility.writeToFile("PROBLEM2BIMP",keyValuePair, 1, colNames);
    map<double,double> THPair1=computeHamilton(keyValuePair);
    utility.writeToFile("PRB2BIMHMT",THPair1,1,colNames);
//    utility.writeToFile("PROBLEM2BRK2P",keyValuePair[1],0,colNames);
    map<double,vector<double>> keyValuePair2=ode.applyIMP(gFB, iv2, T, N);
    utility.writeToFile("PROBLEM2BIMP",keyValuePair2, 2, colNames);
    map<double,double> THPair2=computeHamilton(keyValuePair2);
    utility.writeToFile("PRB2BIMHMT",THPair2,2,colNames);
}


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
        double H=(p1*p1+2.0*p2*p2-2.0*p1*p2*cos(t1-t2))/(2.0*(1.0+sin(t1-t2)*sin(t1-t2)))-2.0*g*cos(t1)-g*cos(t2);
        THPair[x.first]=H;
    }
    return THPair;
}
