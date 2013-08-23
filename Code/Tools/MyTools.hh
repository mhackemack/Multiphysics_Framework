//
//
//
//  A toolkit of functions for global use
//
//

#ifndef mytools_hh_
#define mytools_hh_

// Containers
#include <string>
#include <vector>
#include <map>
//#include <tr1/unordered_map>
//#include <tr1/unordered_set>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iterator>
#include <queue>
#include <stack>
#include <list>
#include <deque>

// Numerical
#include <float.h>
#include <cfloat>
#include <math.h>
#include <cmath>
#include <functional>
#include <numeric>
#include <algorithm>
#include <limits>
#include <random>

// IOS Libraries
#include <ios>
#include <iomanip>

// Stream Libraries
#include <sstream>
#include <istream>
#include <fstream>
#include <ostream>
#include <iostream>

// C/C++ Libraries
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <time.h>

// RTTI
#include <typeinfo>

// Memory tools
#include <memory>

// Exceptions
#include <stdexcept>

// Options
#include <getopt.h>
#include <stdarg.h>



#include "MyHelpers.hh"

const double Avogad = 6.022E+23;

typedef unsigned int uint32;
typedef unsigned long uint64;
typedef unsigned short uint16;

namespace MyTools
{
    inline void MakeLower(std::string&);
    inline void MakeUpper(std::string&);
    inline std::string lower(std::string);
    inline std::string upper(std::string);
    
    inline bool checkpos(std::string,std::string);
    inline bool checkposic(std::string,std::string);
    
    inline int StringToInt(std::string);
    inline double StringToDouble(std::string);
    inline int s2i(std::string);
    inline double s2d(std::string);
    template <typename T>
    inline void StringToNumber(std::string,T&);
    template <typename T>
    inline void s2n(std::string,T&);
    template <typename T>
    inline T StringToNumber(std::string);
    template <typename T>
    inline T s2n(std::string);
    template <typename T>
    inline std::vector<T> s2n(std::vector<std::string>);
    inline std::string IntToString(int);
    inline std::string DoubleToString(double);
    template <typename T>
    inline std::string NumberToString(T);
    template <typename T>
    inline std::string n2s(T);
    template <typename T>
    inline std::vector<std::string> NumberToString(std::vector<T>);
    template <typename T>
    inline std::vector<std::string> n2s(std::vector<T>);
    
    inline std::string Capitalize(std::string);
    inline std::string Uncapitalize(std::string);
    inline void RemoveWhitespace(std::string&);
    inline bool strcmpic(std::string,std::string);
    
    inline std::vector<std::string> delimit(std::string,char);
    inline std::vector<std::string> delimit(std::string,std::string);
    inline std::vector<std::string> delimit(std::vector<std::string>,std::string);
    inline std::string undelimit(std::vector<std::string>,std::string);
    template <typename T>
    inline std::string undelimit(std::vector<T>,std::string);
    inline void undelimit(std::vector<std::string>,char,std::string&);
    inline std::string undelimit(std::vector<std::string>, std::string&,
                                 std::vector<std::string>::iterator,
                                 std::vector<std::string>::iterator);
    template <typename T>
    inline void RemoveElement(int,std::vector<T>&);
    template <typename T>
    inline void RemoveElement(std::vector<T>&,T t);
    
    template <typename T>
    inline void Print(std::ostream&,std::vector<T>,std::string);
    
    template <typename T>
    inline void Print(std::ostream&,std::vector<T>,int);
    
    template <typename T>
    inline void Print(std::ostream&,std::vector<T>,int,int);
    
    template <typename T>
    inline std::vector<T> subvector(const std::vector<T>&,int,int);
    
    inline bool checkforID(std::string,std::string);
    inline bool checkforID(std::string,std::vector<std::string>);
    
    inline std::string removepos(std::string,std::string);
    
    inline void Purge(std::string&,std::string);
    inline void Purge(std::string&,std::string,std::string);
    
    template <typename T>
    inline void RemoveDuplicates(std::vector<T>& vec);
    
    inline bool count(std::string,std::vector<std::string>);
    
    template <typename T>
    inline bool IsDuplicate(T,std::vector<T>);
    
    template <typename T>
    inline void Scale(T&,T,T,T,T);
    
    template <typename T>
    inline void Scale(std::vector<std::vector<T> >&,T,T,T,T);
    
    template <typename T>
    inline T mean(const std::vector<T>&);
    
    template <typename T>
    inline T StdDev(const std::vector<T>&);

    template <typename T>
    inline T RelativeError(const std::vector<T>&);
    
    template <typename T, typename U>
    inline void ExitWithoutParams(T,U,U);
    
    template <typename T, typename U>
    inline void ExitWithoutParams(T,U,U,std::string);
    
    inline bool OpenFile(std::ifstream&,std::string,bool);
    inline void OpenFile(std::ifstream&,std::string);

    inline bool OpenFile(std::ofstream&,std::string,bool);
    inline void OpenFile(std::ofstream&,std::string);
    
    inline void RemoveComment(std::string&,std::string);
    
    inline std::vector<std::string> GetLineAndDelimit(std::ifstream&,std::string);
    inline std::vector<std::string> GetLineAndDelimit(std::ifstream&,std::string,
                                                      std::string);
    inline bool GetLineAndDelimit(std::ifstream&,std::string,
                                  std::vector<std::string>&);
    inline bool GetLineAndDelimit(std::ifstream&,std::string,std::string,
                                  std::vector<std::string>&);
    
    template <typename T>
    inline void ThrowError(T);
    template <typename T,typename U>
    inline void ThrowError(T,U);
    template <typename T,typename U,typename V>
    inline void ThrowError(T,U,V);
    
    inline std::string xmlopen(std::string);
    inline std::string xmlopen(std::string,std::string);
    
    inline std::string xmlclose(std::string);
    inline std::string xmlclose(std::string,std::string);
    
    inline std::string tablevel(int);
    inline int DetermineClose(int,std::vector<int>);
    
    template <typename T>
    inline void AddValue(std::vector<T>&,T);
    
    template <typename T>
    inline T GetMin(const std::vector<T>&);
    
    template <typename T>
    inline T GetMax(const std::vector<T>&);
    
    template <typename T>
    inline T Sum(const std::vector<T>&);
        
    template <typename T,typename U>
    inline void GenerateXMLblock(std::ostream&,const std::vector<T>&,
                                 std::vector<int>,
                                 const std::vector<int>&,
                                 const std::vector<U>&);
    inline void AddLine(std::ifstream&, std::string&, std::string);
    
    template <typename T>
    inline void lineXML(std::ostream&, std::string, T, std::string);
    
    template <typename T>
    inline void lineXML(std::ostream&, std::string, T, std::string, int);
    
    inline bool GetBoolFromString(std::string);
    inline std::vector<bool> GetBoolFromString(const std::vector<std::string>&);
    
    template <typename T>
    inline void RequireSize(const std::vector<T>&,int,std::string);
    
    template <typename T>
    inline bool CheckSize(const std::vector<T>&,int,std::string);
    
    inline std::string GetParam(const std::vector<std::string>&, std::string);
    
    template <typename T>
    inline std::pair<bool,uint32> CheckForParam(const std::vector<T>&, T);
    inline std::pair<bool,uint32> CheckForParam(const std::vector<std::string>&,
                                                std::string);
    
    template <typename T, typename U>
    inline std::vector<T> ConvertVector(const std::vector<U>&);
    
    template <typename T>
    inline std::vector<T> ConvertStringVector(const std::vector<std::string>&);
    
    inline void SetDoublePrecision(double&,int,double);
    
    template <typename T>
    inline typename std::vector<T>::iterator
    GetIteratorPosition(std::vector<T>,T);
    
    template <typename T>
    inline uint32 GetIteratorIndex(std::vector<T>,T);
    
    inline bool AddDirectory(std::string);

    inline std::string MakeDirectoryString(std::string&);
    //inline std::string MakeDirectoryString(const char*&);
    inline std::string MakeDirectoryString(const char[]);
    
    template <typename T>
    inline T interpolate(T,T,T,T,T);
    template <typename T>
    inline T loginterpolate(T,T,T,T,T);
    
    inline double ExplicitEuler(double,double,double);
    inline double ImplicitEuler(double,double,double);
    inline double RungeKutta2nd(double,double,double);
    inline double RungeKutta4th(double,double,double);
    inline double RungeKutta4thUnits(double,double,double);
    
    //template <typename T> inline T factorial(uint32);
};

// ============================================================================
//                                                      Make String Lower
inline void MyTools::MakeLower(std::string& str)
{
    for(std::string::size_type i = 0; i < str.size(); ++i) {
        str[i] = tolower(str[i]);
    }
}
// ============================================================================
//                                                      Make String Upper
inline void MyTools::MakeUpper(std::string& str)
{
    for(std::string::size_type i = 0; i < str.size(); ++i) {
        str[i] = toupper(str[i]);
    }
}
// ============================================================================
//                                                      Make String Lower
inline std::string MyTools::lower(std::string str)
{
    for(std::string::size_type i = 0; i < str.size(); ++i) {
        str[i] = tolower(str[i]);
    }
    return str;
}
// ============================================================================
//                                                      Make String Upper
inline std::string MyTools::upper(std::string str)
{
    for(std::string::size_type i = 0; i < str.size(); ++i) {
        str[i] = toupper(str[i]);
    }
    return str;
}
// ============================================================================
//                                                      Check for Position
inline bool MyTools::checkpos(std::string str, std::string param)
{
    if(str.find(param) != std::string::npos) {
        return true;
    } else {
        return false;
    }
}
// ============================================================================
//                                              Check for Position - Ignore case
inline bool MyTools::checkposic(std::string str, std::string param)
{
    MakeLower(str);
    MakeLower(param);
    if(str.find(param) != std::string::npos) {
        return true;
    } else {
        return false;
    }
}
// ============================================================================
//                                                      String to Int
inline int MyTools::StringToInt(std::string str)
{
    CheckAlpha checkalpha;
    std::find_if(str.begin(),str.end(),checkalpha);
    
    std::istringstream converter(str);
    int result;
    converter >> result;
    return result;
}
// ============================================================================
//                                                      String to Double
inline double MyTools::StringToDouble(std::string str)
{
    CheckAlpha checkalpha;
    std::find_if(str.begin(),str.end(),checkalpha);

    std::istringstream converter(str);
    double result;
    converter >> result;
    return result;
}
// ============================================================================
inline int MyTools::s2i(std::string str) { return StringToInt(str); }
// ============================================================================
inline double MyTools::s2d(std::string str) { return StringToDouble(str); }
// ============================================================================
//                                                      String to Number
template <typename T>
inline void MyTools::StringToNumber(std::string str, T& num)
{
    CheckAlpha checkalpha;
    std::find_if(str.begin(),str.end(),checkalpha);
    if(str.find("\\") == 0) { str.erase(0,1); }
    
    std::istringstream converter(str);
    converter >> num;
}
// ============================================================================
template <typename T>
inline void MyTools::s2n(std::string s, T& n) { StringToNumber(s,n); }
// ============================================================================
//                                                      String to Number
template <typename T>
inline T MyTools::StringToNumber(std::string str)
{
    T value;
    CheckAlpha checkalpha;
    std::find_if(str.begin(),str.end(),checkalpha);
    if(str.find("\\") == 0) { str.erase(0,1); }

    std::istringstream converter(str);
    converter >> value;
    return value;
}
// ============================================================================
//                                                      String to Number
template <typename T>
inline T MyTools::s2n(std::string str)
{
    return StringToNumber<T>(str);
}
// ============================================================================
//                                                      String to Number
template <typename T>
inline std::vector<T> MyTools::s2n(std::vector<std::string> v)
{
    std::vector<T> nv(v.size(),0);
    for(uint32 i = 0; i < v.size(); ++i) {
        nv.at(i) = s2n<T>(v.at(i));
    }
    return nv;
}
// ============================================================================
//                                                      Int to String
inline std::string MyTools::IntToString(int num)
{
    std::stringstream out;
    out << num;
    return out.str();
}
// ============================================================================
//                                                      Double to String
inline std::string MyTools::DoubleToString(double num)
{
    std::stringstream out;
    out << num;
    return out.str();
}
// ============================================================================
//                                                      Number to String
template <typename T>
inline std::string MyTools::NumberToString(T num)
{
    std::stringstream out;
    out << num;
    return out.str();
}
// ============================================================================
template <typename T>
inline std::string MyTools::n2s(T num)
{
    return NumberToString<T>(num);
}
// ============================================================================
template <typename T>
inline std::vector<std::string> MyTools::NumberToString(std::vector<T> num)
{
    std::vector<std::string> v;
    for(uint32 i = 0; i < num.size(); ++i) {
        v.push_back(NumberToString(num[i]));
    }
    return v;
}
// ============================================================================
//                                                      Capitalize String
inline std::string MyTools::Capitalize(std::string str)
{
    if(str.size() > 0) { str[0] = toupper(str[0]); }
    return str;
}
// ============================================================================
//                                                      Uncapitalize String
inline std::string MyTools::Uncapitalize(std::string str)
{
    if(str.size() > 0) { str[0] = tolower(str[0]); }
    return str;
}
// ============================================================================
//                                                      Remove String whitespace
inline void MyTools::RemoveWhitespace(std::string& str)
{
    for(uint32 i = 0; i < str.length(); ++i) {
        if(str[i] == ' ') { str.erase(i,1); }
    }
    
    if(str.length() == 0) {
        std::cout << "Warning tools::RemoveWhitespace() "
        << "was passed an empty string" << std::endl;
    }
}
// ============================================================================
//                                                  String Compare - Ignore Case
inline bool MyTools::strcmpic(std::string str1, std::string str2)
{
    MakeLower(str1);
    MakeLower(str2);
    RemoveWhitespace(str1);
    RemoveWhitespace(str2);
    if(str1 == str2) { return true; }
    else { return false; }
}
// ============================================================================
//                                                      Delimit Line
inline std::vector<std::string> MyTools::delimit(std::string line, char delimiter)
{
    std::vector<std::string> delimitedLine;
    std::string token;
    std::istringstream iss(line);
    while( getline(iss,token,delimiter) ) {
        if( !token.empty() ) { delimitedLine.push_back(token); }
    }
    return delimitedLine;
}
// ============================================================================
//                                                      Delimit Line
inline std::vector<std::string> MyTools::delimit(std::string line, std::string delimiter)
{
    for(std::string::size_type j = 0; j < delimiter.length(); ++j) {
        for(std::string::size_type i = 0; i < line.length(); ++i) {
            if(line[i] == delimiter[j]) { line[i] = ' '; }
        }
    }
    return delimit(line,' ');
    
}
// ============================================================================
//                                                      Delimit Line
inline std::vector<std::string> MyTools::delimit(std::vector<std::string> line,
                                                 std::string delimiter)
{
    std::vector<std::string> newDelimited;
    for(uint32 i = 0; i < line.size(); ++i) {
        std::vector<std::string> newline = delimit(line[i],delimiter);
        for(uint32 j = 0; j < newline.size(); ++j) {
            newDelimited.push_back(newline[j]);
        }
    }
    return newDelimited;
}
// ============================================================================
//                                                      Undelimit Line
template <typename T>
inline std::string MyTools::undelimit(std::vector<T> delimV, std::string delimiter)
{
    std::vector<std::string> delimStr = NumberToString(delimV);
    std::string str = "";
    for(uint32 i = 0; i < delimStr.size(); ++i) {
        str += delimStr[i];
        if(i < delimStr.size()-1) { str += delimiter; }
    }
    return str;
}
// ============================================================================
//                                                      Undelimit Line
inline std::string MyTools::undelimit(std::vector<std::string> delimStr, std::string delimiter)
{
    std::string str = "";
    for(uint32 i = 0; i < delimStr.size(); ++i) {
        str += delimStr[i];
        if(i < delimStr.size()-1) { str += delimiter; }
    }
    return str;
}
// ============================================================================
//                                                      Undelimit Line
inline void MyTools::undelimit(std::vector<std::string> dl ,char delim, std::string& s)
{
    s="";
    if(dl.size() == 0) { return; }
    for(uint32 i = 0; i < dl.size(); ++i) {
        s += dl[i];
        if(i < dl.size()-1) { s += delim; }
    }
}
// ============================================================================
//                                                      Undelimit Line
inline std::string MyTools::undelimit(std::vector<std::string> dl, std::string& s,
                                      std::vector<std::string>::iterator start,
                                      std::vector<std::string>::iterator stop)
{
    std::string str = "";
    for(uint32 i = std::distance(dl.begin(),start);
        i < std::distance(dl.begin(),stop)+1; ++i) {
        str += dl.at(i);
        if(i < dl.size()-1) { str += s; }
    }
    return str;
}
// ============================================================================
//                                                      Remove Vector Element
template <typename T>
void MyTools::RemoveElement(int pos, std::vector<T>& vec)
{
    vec.erase(vec.begin()+pos);
}

// ============================================================================
//                                                      Remove Vector Element
template <typename T>
void MyTools::RemoveElement(std::vector<T>& vec, T t)
{
    typename std::vector<T>::iterator ite;
    for(ite = vec.begin(); ite != vec.end(); ++ite) {
        if( *ite == t ) {
            vec.erase(ite);
            break;
        }
    }
}
// ============================================================================
//                                                      Print Vector
template <typename T>
inline void MyTools::Print(std::ostream& out, std::vector<T> vec, std::string c)
{
    for(uint32 i = 0; i < vec.size(); ++i) {
        out << vec.at(i);
        if(i+1 != vec.size()) { out << c; }
    }
    
}
// ============================================================================
//                                                      Print Vector
template <typename T>
inline void MyTools::Print(std::ostream& out, std::vector<T> vec, int width)
{
    out.width(width);
    for(uint32 i = 0; i < vec.size(); ++i) {
        out << std::setw(width) << vec.at(i);
    }
}
// ============================================================================
//                                                      Print Vector
template <typename T>
inline void MyTools::Print(std::ostream& out, std::vector<T> vec, int width, int width2)
{
    out.width(width);
    out.precision(width2);
    for(uint32 i = 0; i < vec.size(); ++i) {
        out << std::setw(width) << std::setprecision(width2) << vec.at(i);
    }
}
// ============================================================================
//                                                      Get portion of vector
template <typename T>
inline std::vector<T> MyTools::subvector(const std::vector<T>& v,int first, int last)
{
    std::vector<T> nv;
    for(uint32 i = first; i <= last; ++i) {
        if(i >= v.size()) { break; }
        nv.push_back(v.at(i));
    }
    return nv;
}
// ============================================================================
//                                                      Check for ID
inline bool MyTools::checkforID(std::string str, std::string param)
{
    for(std::string::size_type i = 0; i < param.size(); ++i) {
        if(str.find(param[i]) == std::string::npos) { return false; }
    }
    return true;
}
// ============================================================================
//                                                      Check for ID
inline bool MyTools::checkforID(std::string str, std::vector<std::string> param)
{
    for(uint32 j = 0; j < param.size(); ++j) {
        for(std::string::size_type i = 0; i < param[j].size(); ++i) {
            if(str.find(param[i]) != std::string::npos) { return true; }
        }
    }
    return false;
}
// ============================================================================
//                                                      Remove Position
inline std::string MyTools::removepos(std::string orig, std::string param)
{
    return (orig.find(param) != std::string::npos) ? orig.erase(orig.find(param),param.size()) : orig;
}
// ============================================================================
//                                                      Purge Line
inline void MyTools::Purge(std::string& str,std::string delimiter)
{
    bool fixed = false;
    while (!fixed) {
        fixed = true;
        if(checkposic(str,delimiter)) {
            str.replace(str.find(delimiter),delimiter.size()," ");
            fixed = false;
        }
    }
}
// ============================================================================
//                                                      Purge Line
inline void MyTools::Purge(std::string& str,std::string delimiter, std::string replace)
{
    bool fixed = false;
    while (!fixed) {
        fixed = true;
        if(checkposic(str,delimiter)) {
            str.replace(str.find(delimiter),delimiter.size(),replace);
            fixed = false;
        }
    }
}
// ============================================================================
//                                                      Remove Duplicates
template<typename T>
inline void MyTools::RemoveDuplicates(std::vector<T>& vec)
{
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}
// ============================================================================
//                                                      Count
inline bool MyTools::count(std::string str,std::vector<std::string> allstrs)
{
    for(uint32 i = 0; i < allstrs.size(); ++i) {
        if((lower(allstrs[i])).find(lower(str)) != std::string::npos) {
            return true;
        }
    }
    return false;
}
// ============================================================================
//                                                      Is Duplicate
template <typename T>
inline bool MyTools::IsDuplicate(T val, std::vector<T> v)
{
    for(uint32 i = 0; i < v.size(); ++i) {
        if(v.at(i) == val) { return true; }
    }
    return false;
}
// ============================================================================
//                                                      Scale number
template <typename T>
inline void MyTools::Scale(T& val, T min, T max, T a, T b)
{
    // scale val from [min,max] to [a,b]
    val = (b-a)*(val-min)/(max-min)+a;
}
// ============================================================================
//                                                      Scale range
template <typename T>
inline void MyTools::Scale(std::vector<std::vector<T> >& v, T min, T range, T scaleMin, T scaleMax)
{
    T scaleRange = scaleMax - scaleMin;
    
    double avg = 0;
    assert(range != 0);
    for(uint32 j = 0; j < v.size(); ++j) {
        for(uint32 i = 0; i < v[j].size(); ++i) {
            v[j][i] = ((scaleRange*(v[j][i]-min))/range) + scaleMin;
            double inverse = fabs(scaleMax-v[j][i])/(v[j][i]);
            if(j > 4 && fabs(inverse-avg) > 10000.) { inverse = avg; }
            else { avg = (avg+inverse)/2.; }
            v[j][i] *= pow(inverse,2);
            if(fabs(v[j][i] - avg) > 1e6) { v[j][i] *= 1.e6; }
        }
    }
}
// ============================================================================
//                                                      Mean
template <typename T>
inline T MyTools::mean(const std::vector<T>& v) {
    T accum = std::accumulate(v.begin(),v.end(),accum);
    return accum/static_cast<T>(v.size());
}
// ============================================================================
//                                                      Standard Deviation
template <typename T>
inline T MyTools::StdDev(const std::vector<T>& v)
{
    T mean = mean(v);
    T var = 0;
    for(uint32 i = 0; i < v.size(); ++i) {
        var += pow(v.at(i)-mean,2);
    }
    var /= (static_cast<T>(v.size())-static_cast<T>(1));
    return sqrt(var);
}
// ============================================================================
//                                                      Standard Deviation
template <typename T>
inline T MyTools::RelativeError(const std::vector<T>& v) {
    T mean = mean(v);
    T relErr = 0;
    for(uint32 i = 0; i < v.size(); ++i) {
        relErr += pow(v.at(i)-mean,2);
    }
    relErr /= (static_cast<T>(v.size())-1);
    relErr = sqrt(relErr);
    return relErr / mean / sqrt(static_cast<T>(v.size()));
}
// ============================================================================
//                                                      Exit w/o params
template <typename T, typename U>
inline void MyTools::ExitWithoutParams(T arg, U args, U nParams) {
    if(args < nParams) {
        std::cout << arg << " needs " << nParams << " to continue"
        << ". Exitting...\n" << std::endl;
        exit(-1);
    }
}
// ============================================================================
//                                                      Exit w/o params
template <typename T, typename U>
inline void MyTools::ExitWithoutParams(T arg, U args, U nParams, std::string message) {
    if(args < nParams) {
        std::cout << arg << " needs " << nParams << " to continue"
        << ". Please specify " << message << " in the command line.\n" << std::endl;
        exit(-1);
    }
}
// ============================================================================
//                                                      Open File
inline bool MyTools::OpenFile(std::ifstream& inf, std::string fname, bool doExit) {
    inf.open(fname.c_str());
    if(!inf) {
        if(doExit) { Exception("Tools","OpenFile - unable to open file",fname,FATAL); }
        else {
            Exception("Tools","OpenFile - unable to open file",fname,WARNING);
            return false; }
    } else {
        return true;
    }
}
// ============================================================================
//                                                      Open File
inline void MyTools::OpenFile(std::ifstream& inf, std::string fname) {
    inf.open(fname.c_str());
    if(!inf) {
        Exception("Tools","OpenFile - unable to open file",fname,FATAL);
    }
}
// ============================================================================
//                                                      Open File
inline bool MyTools::OpenFile(std::ofstream& inf, std::string fname, bool doExit) {
    inf.open(fname.c_str());
    if(!inf) {
        if(doExit) { Exception("Tools","OpenFile - unable to open file",fname,FATAL); }
        else {
            Exception("Tools","OpenFile - unable to open file",fname,WARNING);
            return false;
        }
    } else {
        return true;
    }
}
// ============================================================================
//                                                      Open File
inline void MyTools::OpenFile(std::ofstream& inf, std::string fname) {
    inf.open(fname.c_str());
    if(!inf) {
        Exception("Tools","OpenFile - unable to open file",fname,FATAL);
    }
}
// ============================================================================
//                                                      Remove Comment
inline void MyTools::RemoveComment(std::string& line,std::string commentTag)
{
    if(checkposic(line,commentTag)) {
        line.erase(lower(line).find(commentTag));
    }
}
// ============================================================================
//                                                      Get Line and Delimit
inline std::vector<std::string> MyTools::GetLineAndDelimit(std::ifstream& inf, std::string str) {
    while(!inf.eof()) {
        std::string line;
        getline(inf,line,'\n');
        std::vector<std::string> delimited = delimit(line,str);
        if(!delimited.empty()) { return delimited; }
    }
    std::vector<std::string> d;
    std::cout << "WARNING! Function GetLineAndDelimit reached end of file and ";
    std::cout << "is returning empty vector" << std::endl;
    return d;
}
// ============================================================================
//                                                      Get Line and Delimit
inline std::vector<std::string> MyTools::GetLineAndDelimit(std::ifstream& inf, std::string str,
                                                           std::string commentTag) {
    while(!inf.eof()) {
        std::string line;
        getline(inf,line,'\n');
        RemoveComment(line,commentTag);
        std::vector<std::string> delimited = delimit(line,str);
        if(!delimited.empty()) { return delimited; }
    }
    std::vector<std::string> d;
    std::cout << "WARNING! Function GetLineAndDelimit reached end of file and ";
    std::cout << "is returning empty vector" << std::endl;
    return d;
}
// ============================================================================
//                                                      Get Line and Delimit
inline bool MyTools::GetLineAndDelimit(std::ifstream& inf, std::string str,
                                       std::string commentTag, std::vector<std::string>& v) {
    v.clear();
    while(!inf.eof()) {
        std::string line;
        getline(inf,line,'\n');
        RemoveComment(line,commentTag);
        v = delimit(line,str);
        if(!v.empty()) { return true; }
    }
    return false;
}
// ============================================================================
//                                                      Get Line and Delimit
inline bool MyTools::GetLineAndDelimit(std::ifstream& inf, std::string str,
                                       std::vector<std::string>& v) {
    v.clear();
    while(!inf.eof()) {
        std::string line;
        getline(inf,line,'\n');
        v = delimit(line,str);
        if(!v.empty()) {
            return true;
        }
    }
    return false;
}
// ============================================================================
//                                                      Throw Error
template <typename T>
inline void MyTools::ThrowError(T expression) {
    std::cout << "Error! " << expression
    << "\nExiting.." << std::endl;
    exit(1);
}
// ============================================================================
//                                                      Throw Error
template <typename T,typename U>
inline void MyTools::ThrowError(T expression, U title) {
    std::cerr << "Error (" << expression << ") with " << title
    << "\nExiting..." << std::endl;
    exit(1);
}
// ============================================================================
//                                                      Throw Error
template <typename T,typename U,typename V>
inline void MyTools::ThrowError(T expression, U section, V title) {
    std::cerr << "Error (" << expression << ") with " << section
    << " of " << title << "\nExiting..." << std::endl;
    exit(1);
}
// ============================================================================
//                                                      XML Open
inline std::string MyTools::xmlopen(std::string str,std::string space)
{
    return space + "<" + str + ">";
}
// ============================================================================
//                                                      XML Open
inline std::string MyTools::xmlopen(std::string str)
{
    return "<" + str + ">";
}
// ============================================================================
//                                                      XML Close
inline std::string MyTools::xmlclose(std::string str,std::string space)
{
    return space + "</" + str + ">";
}
// ============================================================================
//                                                      XML Close
inline std::string MyTools::xmlclose(std::string str)
{
    return "</" + str + ">";
}
// ============================================================================
//                                                      Tab level
inline std::string MyTools::tablevel(int level) {
    std::string str = "";
    if(level < 0) { level = 0; }
    for(int i = 0; i < level; ++i) {
        str += "\t";
    }
    return str;
}
// ============================================================================
//                                                      Determine Close
inline int MyTools::DetermineClose(int pos,std::vector<int> vec)
{
    std::vector<int>::iterator ite = std::find(vec.begin(),vec.end(),pos);
    if(ite == vec.end()) { return -1; }
    else { return std::distance(vec.begin(),ite); }
}
// ============================================================================
//                                                      Add value
template <typename T>
inline void MyTools::AddValue(std::vector<T>& v, T val)
{
    for(uint32 i = 0; i < v.size(); ++i) {
        v[i] += val;
    }
}
// ============================================================================
//                                                      Get Min
template <typename T>
T MyTools::GetMin(const std::vector<T>& v) {
    return *std::min_element(v.begin(),v.end());
}
// ============================================================================
//                                                      Get Max
template <typename T>
T MyTools::GetMax(const std::vector<T>& v) {
    return *std::max_element(v.begin(),v.end());
}
// ============================================================================
//                                                      Get Sum
template <typename T>
T MyTools::Sum(const std::vector<T>& v) {
    T value;
    return std::accumulate(v.begin(),v.end(),value);
}
// ============================================================================
//                                                      Generate XML block
template <typename T, typename U>
inline void MyTools::GenerateXMLblock(std::ostream& out, const std::vector<T>& tags,
                                      std::vector<int> levels,
                                      const std::vector<int>& subsets,
                                      const std::vector<U>& values)
{
    if(tags.size() != values.size() ||
       tags.size() != levels.size() ||
       tags.size() != subsets.size() ) {
        ThrowError("incompatible sizes for GenerateXMLblock"); }
    
    std::vector<T> closeTags;
    std::vector<int> closeTagsLoc;
        
    std::vector<int>::iterator upper;
    upper = std::upper_bound(levels.begin(),levels.end(),0);
    
    if(upper == levels.end()) {
        ThrowError("no master xml tag");
    } else {
        closeTags.push_back(tags[std::distance(levels.begin(),upper)]);
        closeTagsLoc.push_back(tags.size()+1);
    }
    
    for(uint32 i = 0; i < tags.size(); ++i) {
        if(DetermineClose(i,closeTagsLoc) > 0) {
            int index = DetermineClose(i,closeTagsLoc);
            out << tablevel(levels[index]) << xmlclose(closeTags[index])
            << std::endl << std::endl;
            closeTags.pop_back();
            closeTagsLoc.pop_back();
        }
        if(subsets[i] == 0) {
            out << tablevel(levels[i]) << xmlopen(tags[i])
            << values[i] << xmlclose(tags[i]) << std::endl;
        } else if( i != 0 && levels[i] != 1 ) {
            closeTags.push_back(tags[i]);
            int newLoc = i+subsets[i]+1;
            if(newLoc == closeTagsLoc[i-1]) { newLoc -= 1; }
            closeTagsLoc.push_back(newLoc);
            out << tablevel(levels[i]) << xmlopen(tags[i]) << std::endl;
        }
        if ( i == 0 ) {
            out << tablevel(levels[i]) << xmlopen(tags[i]) << std::endl;
        }
    }
    
    if(closeTags.size() > 0) {
        for(uint32 i = tags.size(); i < 2*tags.size(); ++i) {
            if(closeTags.size() == 0) { break; }
            int index = DetermineClose(i,closeTagsLoc);
            if(index < 0) { continue; }
            out << tablevel(levels[index]) << xmlclose(closeTags[index])
            << std::endl << std::endl;
            closeTags.pop_back();
            closeTagsLoc.pop_back();
        }
    }
    
}
// ============================================================================
//                                                      Add Line
inline void MyTools::AddLine(std::ifstream& in, std::string& prev, std::string tag)
{
    if(prev.find(tag) == std::string::npos) { return; }
    else {
        while(true) {
            prev.erase(prev.find(tag));
            std::string line;
            getline(in,line,'\n');
            prev += line;
            if(lower(prev).find(tag) == std::string::npos) { break; }
        }
    }
}
//===================================================================================
//                                                      Line XML
template <typename T>
inline void MyTools::lineXML(std::ostream& out, std::string tag, T value, std::string delim)
{
    out << delim << xmlopen(tag) << value << xmlclose(tag) << std::endl;
}
//===================================================================================
//                                                      Line XML
template <typename T>
inline void MyTools::lineXML(std::ostream& out, std::string tag, T value, std::string delim, int prec)
{
    out << delim << xmlopen(tag) << std::setprecision(prec) << value << xmlclose(tag) << std::endl;
}
// ============================================================================
//                                                      Get Bool From String
inline bool MyTools::GetBoolFromString(std::string str)
{
    MakeLower(str);
    if(str.length() == 1) {
        if(str == "t" || str == "y") { return true; }
        else { return false; }
    } else {
        if(str.find("true") != std::string::npos ||
           str.find("yes") != std::string::npos ) { return true; }
        else { return false; }
    }
}
// ============================================================================
//                                                      Get Bool From String
inline std::vector<bool> MyTools::GetBoolFromString(const std::vector<std::string>& str)
{
    std::vector<bool> boolstr(str.size(),false);
    for(uint32 i = 0; i < str.size(); ++i) {
        boolstr.at(i) = GetBoolFromString(str.at(i));
    }
    return boolstr;
}
// ============================================================================
//                                                      Require size
template <typename T>
inline void MyTools::RequireSize(const std::vector<T>& v, int i, std::string message) {
    if(i == v.size()) { return; }
    else {
        int diff = i - v.size();
        std::cout << "Vector ";
        for(uint32 ii = 0; ii < v.size(); ++ii) {
            std::cout << v[ii] << " ";
        }
        std::cout << std::endl;
        std::string cond = "too ";
        cond += (diff > 0) ? "small " : "large ";
        cond += "by " + NumberToString(abs(diff));
        ThrowError("bad number of parameters",message,cond);
    }
}
// ============================================================================
//                                                      Check Size
template <typename T>
inline bool MyTools::CheckSize(const std::vector<T>& v, int i, std::string message) {
    if(int(v.size()) >= i) { return true; }
    else {
        int diff = i - v.size();
        std::cout << "Vector ";
        for(uint32 ii = 0; ii < v.size(); ++ii) {
            std::cout << v[ii] << " ";
        }
        std::cout << std::endl;
        std::string cond = "too small by " + NumberToString(diff);
        ThrowError("bad number of parameters",message,cond);
    }
    return false;
}
// ============================================================================
//                                                      Get Param
inline std::string MyTools::GetParam(const std::vector<std::string>& v, std::string param)
{
    for(uint32 i = 0; i < v.size(); ++i) {
        if(strcmpic(v[i],param)) {
            CheckSize(v,i+1,"looking for "+param);
            return v[i+1];
        }
    }
    ThrowError("could not find requested parameter",param,undelimit(v," "));
    return "";      // gets rid of warnings
}
// ============================================================================
//                                                      Check For Param
template <typename T>
inline std::pair<bool,uint32> MyTools::CheckForParam(const std::vector<T>& v, T param)
{
    for(uint32 i = 0; i < v.size(); ++i) {
        if(v[i] == param) { return std::pair<bool,uint32>(true,i); }
    }
    return std::pair<bool,uint32>(false,0);
}
// ============================================================================
//                                                      Check For Param
inline std::pair<bool,uint32> MyTools::CheckForParam(const std::vector<std::string>& v,
                                                     std::string param)
{
    for(uint32 i = 0; i < v.size(); ++i) {
        if(strcmpic(v[i],param)) {
            return std::pair<bool,uint32>(true,i);
        }
    }
    return std::pair<bool,uint32>(false,0);
}
// ============================================================================
//                                                      Convert Vector
// ============================================================================
template<typename T, typename U>
inline std::vector<T> MyTools::ConvertVector(const std::vector<U>& v)
{
    std::vector<T> tv;
    for(uint32 i = 0; i < v.size(); ++i) {
        tv.push_back(static_cast<T>(v.at(i)));
    }
    return tv;
}
// ============================================================================
//                                          Convert String Vector to Numbers
// ============================================================================
template<typename T>
inline std::vector<T> MyTools::ConvertStringVector(const std::vector<std::string>& v) {
    std::vector<T> tv;
    for(uint32 i = 0; i < v.size(); ++i) {
        tv.push_back(s2n<T>(v.at(i)));
    }
    return tv;
}
// ============================================================================
//                                      Truncate after certain decimal place
// ============================================================================
inline void MyTools::SetDoublePrecision(double& value, int precision, double units)
{
    uint16 fac = (value >= 0) ? 1 : -1;
    value = fabs(value);
    uint64 tmp = static_cast<uint64>((value/units) * pow(10,precision));
    value = (static_cast<double>(tmp) * pow(10.,-precision))*units;
    value *= static_cast<double>(fac);
}
// ============================================================================
//                                                  Get Iterator at Position
// ============================================================================
template <typename T>
inline typename std::vector<T>::iterator
MyTools::GetIteratorPosition(std::vector<T> v, T criteria)
{
    typename std::vector<T>::iterator it;
    for(it = v.begin(); it != v.end(); ++it) {
        if( *it == criteria ) { return it; }
    }
    ThrowError("GetIteratorPosition - could not find component matching critera",
               criteria);
    return v.end();
}
// ============================================================================
//                                                  Get Iterator Index
// ============================================================================
template <typename T>
inline uint32 MyTools::GetIteratorIndex(std::vector<T> v, T criteria)
{
    typename std::vector<T>::iterator it;
    for(it = v.begin(); it != v.end(); ++it) {
        if( *it == criteria ) { return std::distance(v.begin(),it); }
    }
    ThrowError("GetIteratorIndex - could not find component matching critera",
               criteria);
    return v.size();
}
// ============================================================================
//                                                  Add directory
// ============================================================================
inline bool MyTools::AddDirectory(std::string str)
{
    int r = 0;
    std::string cmd = "if [ ! -d " + str + " ]; then mkdir -p " + str + "; fi";
    if(system(NULL)) {
        r = system(cmd.c_str());
        if(r != 0) {
            std::cout << "\n\tWARNING!!! Return value " << r << " != 0 typically \
            indicates " << cmd << " was not successful.\n" << std::endl;
            sleep(2);
            return false;
        } else {
            return true;
        }
    } else {
        std::cout << "\n\tWARNING!!! Command " << cmd << " was not executed because \
        command processor \"system\" was not available.\n" << std::endl;
        sleep(2);
        return false;
    }
    return false;
}
// ============================================================================
//                                                  Make directory string
// ============================================================================
inline std::string MyTools::MakeDirectoryString(std::string& str)
{
    if(str.find("/") != str.length()-1) { return str += "/"; }
    else { return str; }
}
// ============================================================================
//                                                  Make directory string
// ============================================================================
inline std::string MyTools::MakeDirectoryString(const char cstr[])
{
    std::string str = std::string(cstr);
    return MakeDirectoryString(str);
}
// ============================================================================
//                                                  Interpolate for new value
// ============================================================================
template <typename T>
inline T MyTools::interpolate (T x_new, T x1, T x2, T y1, T y2)
{
    if(!(x_new > x1 && x_new < x2 )) {
        std::cout << "x_new = " << x_new << ", between " << x1 << " and " << x2 << std::endl;
        Exception("Interpolate","Interpolating for value not between bounds",n2s(x_new),Warning);
    }
    return y2 - ((y2-y1)/(x2-x1))*(x2 - x_new);
}
// ============================================================================
//                                                  Interpolate for new value
// ============================================================================
template <typename T>
inline T MyTools::loginterpolate (T x_new, T x1, T x2, T y1, T y2)
{
    if(!(x_new > x1 && x_new < x2 )) {
        std::cout << "x_new = " << x_new << ", between " << x1 << " and " << x2 << std::endl;
        Exception("LogInterpolate","Interpolating for value not between bounds",n2s(x_new),Warning);
    }
    return log(y2) - ((log(y2)-log(y1))/(log(x2)-log(x1)))*(log(x2) - log(x_new));
}
// ============================================================================
//                                                  Explicit Euler Integration
// ============================================================================
inline double MyTools::ExplicitEuler(double x_, double dydx_, double h_)
{
    return x_ + h_*dydx_;
}
// ============================================================================
//                                                  Implicit Euler Integration
// ============================================================================
inline double MyTools::ImplicitEuler(double x_, double dydx_, double h_)
{
    return x_ + h_*(dydx_ + 0.5*h_*dydx_);
}
// ============================================================================
//                                                  2nd order Runge-Kutta
// ============================================================================
inline double MyTools::RungeKutta2nd(double x_, double dydx_, double h_)
{
    // Easier to read version:
    // k1 = h_*dydx_
    // k2 = h_*(dydx_ + 0.5*k1)
    // x_new = x_ + (k1 + k2)*0.5
    //
    // Version with optimization of temporaries removal
    return x_ + h_*(dydx_ + (dydx_ + 0.5*h_*dydx_))*0.5;
}
// ============================================================================
//                                                  4th order Runge-Kutta
// ============================================================================
inline double MyTools::RungeKutta4th(double x_, double dydx_, double h_)
{
    // Easier to read version:
    // auto k1 = h_*dydx_;
    // auto k2 = h_*(dydx_ + 0.5*k1);
    // auto k3 = h_*(dydx_ + 0.5*k2);
    // auto k4 = h_*(dydx_ + k3);
    // x_new = x_ + (k1 + 2*k2 + 2*k3 + k4)/6
    //
    // Version with optimization of temporaries removal
    return x_ + ( (h_*dydx_) + 2*(h_*(dydx_ + 0.5*(h_*dydx_))) +
                 2*(h_*(dydx_ + 0.5*(h_*(dydx_ + 0.5*(h_*dydx_))))) +
                 (h_*(dydx_ + (h_*(dydx_ + 0.5*(h_*(dydx_ + 0.5*(h_*dydx_))))))))/6;
}
// ============================================================================
//                                                  4th order Runge-Kutta
// ============================================================================
inline double MyTools::RungeKutta4thUnits(double x_, double dydx_, double h_)
{
    // Easier to read version:
    // k1 = h_*dydx_
    // k2 = h_*(dydx_ + 0.5*k1)
    // k3 = h_*(dydx_ + 0.5*k2)
    // k4 = h_*(dydx_ + k3)
    // x_new = x_ + (k1 + 2*k2 + 2*k3 + k4)/6
    //
    // Version with optimization of temporaries removal
    
    // TEST VERSION (UNITS)
    
    
    auto k1 = h_*dydx_;
    auto k2 = 0.5*(k1 + h_*dydx_);
    auto k3 = 0.5*(k2 + h_*dydx_);
    auto k4 = k1 + h_*dydx_;
    
    return x_ + (k1 + 2*k2 + 2*k3 + k4)/6;
    
    /*
    return x_ + ( (h_*dydx_) + 2*(h_*(dydx_ + 0.5*(dydx_))) +
                 2*(h_*(dydx_ + 0.5*((dydx_ + 0.5*(dydx_))))) +
                 (h_*(dydx_ + ((dydx_ + 0.5*((dydx_ + 0.5*(dydx_))))))))/6;
     */

}
// ============================================================================
//                                                  Factorial
// ============================================================================
/*template <typename T>
inline T MyTools::factorial(uint32 number)
{
    //T factor = static_cast<T>(U);
    //while (factor > 0. ) {
    //    val *= (--factor);
    //}
    //return val;
    
    T temp;
    
    if(number <= 1) return 1;
    
    temp = number * factorial(number - 1);
    return temp;    

}*/

#endif







