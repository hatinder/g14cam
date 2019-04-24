//
// Created by hsingh9 on 17/04/2019.
//

#ifndef CODE_HELPER_HPP
#define CODE_HELPER_HPP

#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <cassert>

using namespace std;

//L2 NORM GENERIC function for vector
template<class T>
T l2_norm (const vector<T> &x)
{
    double result = 0.0;
    for (int i = 0; i < x.size(); ++i)
    {
        result += x[i] * x[i];
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

template<typename T>
vector<T> operator* (T &a, const vector<T> &v)
{
    vector<T> results(v.size());
    for (int i = 0; i < v.size(); ++i)
    {
        results[i] = a * v[i];
    }
    return results;
}

template<typename T>
vector<T> operator/ (const vector<T> &lhs, T &rhs)
{
    vector<T> results(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
    {
        results[i] = lhs[i] / rhs;
    }
    return results;
}

template<typename T>
vector<T> operator+ (const vector<T> &lhs, const vector<T> &rhs)
{
    assert(lhs.size() == rhs.size());
    vector<T> results(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
    {
        results[i] = lhs[i] + rhs[i];
    }
    return results;
}

template<typename T>
vector<T> operator- (const vector<T> &lhs, const vector<T> &rhs)
{
    assert(lhs.size() == rhs.size());
    vector<T> results(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
    {
        results[i] = lhs[i] - rhs[i];
    }
    return results;
}

template<typename T>
T findByValue(const map<T,T> &m,T &value)
{
    T result;
    for (auto it=m.begin(); it!=m.end(); ++it)
    {
        if(abs(it->second)>=value)
        {
            result = it->first;
            break;
        }
    }
    return result;
}
#endif //CODE_HELPER_HPP
