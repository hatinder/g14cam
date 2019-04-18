//
// Created by hsingh9 on 17/04/2019.
//

#ifndef CODE_HELPER_HPP
#define CODE_HELPER_HPP

#include <vector>
#include <cmath>

using namespace std;

//L2 NORM GENERIC function for vector
template<class T> T l2_norm (const vector<T>& x)
{
    double result=0.0;
    for (int i = 0; i < x.size(); ++i)
    {
        result+=x[i]*x[i];
    }
    return sqrt(result);
}

//cout operator for formatted vector output
template<typename T>
ostream &operator<< (ostream &out, const vector<T> &v)
{
    for (int i = 0; i < v.size(); ++i)
    { out << setw(12) << v[i] << endl; }
    return out;
}

//cout operator for formatted map output
template<typename T>
ostream &operator<< (ostream &out, const map<T, T> &m)
{
    for (const pair<T, T> p:m)
    { out << setw(12) << p.first << " : " << setw(12) << p.second << endl; }
    return out;
}

#endif //CODE_HELPER_HPP