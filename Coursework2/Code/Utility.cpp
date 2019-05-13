//
// Created by hsingh9 on 17/04/2019.
//

#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include "Utility.hpp"

vector<double> Utility::LinearSpaced(const double min, const double max, const int size)
{
    vector<double> v;
    for (int i = 0; i < size; ++i)
    {
        v.push_back(min + i * ((max - min) / (size - 1)));
    }
    return v;
}

void
Utility::writeToFile(const string fNamePrefix, map<double, double> m, const int k, vector<string> colNames)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    for (const auto &v:colNames)
        oFileStream << setw(12) << v;
    oFileStream << endl;
    for (auto v:m)
        oFileStream << setw(12) << v.first << setw(12) << v.second << endl;
    oFileStream.close();
}

void Utility::writeToFile(string fNamePrefix, map<int, double> m, int k, vector<string> colNames)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    for (const auto &v:colNames)
        oFileStream << setw(12) << v;
    oFileStream << endl;
    for (auto v:m)
        oFileStream << setw(12) << v.first << setw(12) << v.second << endl;
    oFileStream.close();
}

void Utility::writeToFile(string fNamePrefix, VectorXd vxd, int k, vector<string> colNames)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    for (const auto &v:colNames)
        oFileStream << setw(12) << v;
    oFileStream << endl;
    for (int i = 0; i < vxd.rows(); ++i)
    {
        oFileStream << vxd(i);
        oFileStream << endl;
    }
    oFileStream.close();
}

void Utility::writeToFile(string fNamePrefix, MatrixXd mxd, int k)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    oFileStream << mxd;
    oFileStream.close();
}

void Utility::writeToFile(string fNamePrefix, ArrayXXd axd, int k, vector<string> colNames)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    for (const auto &v:colNames)
        oFileStream << setw(12) << v;
    oFileStream << endl;
    for (int i = 0; i < axd.rows(); ++i)
    {
        for (int j = 0; j < axd.cols(); ++j)
        {
            oFileStream << axd(i, j);
        }
        oFileStream << endl;
    }
    oFileStream.close();
}

void Utility::writeToFile(string fNamePrefix, map<double, vector<double>> m, int k, vector<string> colNames)
{
    ostringstream iterate_label;
    iterate_label.width(3);
    iterate_label.fill('0');
    iterate_label << k;
    string file_name = fNamePrefix + iterate_label.str() + ".txt";
    ofstream oFileStream;
    oFileStream.open(file_name.c_str());
    assert(oFileStream.is_open());
    for (const auto &v:colNames)
        oFileStream << setw(12) << v;
    oFileStream << endl;
    for (auto v:m)
    {
        oFileStream << setw(12) << v.first;
        for(auto w:m[v.first])
        {
            oFileStream << setw(12) << w;
        }
        oFileStream<<endl;
    }
    oFileStream.close();
}

