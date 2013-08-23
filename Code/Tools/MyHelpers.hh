//
//
//
//
//
//


#ifndef myhelpers_hh_
#define myhelpers_hh_

#include <stdexcept>
#include <string>

//================================================================

enum MyExceptionSeverity
{
    Fatal,
    FATAL,
    Warning,
    WARNING
};

//================================================================

class Exception
{
public:
    Exception(std::string ss, MyExceptionSeverity severity)
    {
        std::cerr << "Exception: " << ss;
        switch(severity) {
            case Fatal:
            case FATAL:
                std::cerr << " - FATAL ERROR\n" << std::endl;
                abort();
                break;
            case Warning:
            case WARNING:
                std::cerr << " - WARNING" << std::endl;
                break;
        }
    }
    
    Exception(std::string ss1, std::string ss2, MyExceptionSeverity severity)
    {
        std::cerr << "Exception: " << ss2 << " **[ Issued by " << ss1 << " ]**";
        switch(severity) {
            case Fatal:
            case FATAL:
                std::cerr << " - FATAL ERROR\n" << std::endl;
                abort();
                break;
            case Warning:
            case WARNING:
               std::cerr << " - WARNING" << std::endl;
                break;
        }
    }
    
    template <typename T>
    Exception(std::string ss1, std::string ss2, T ss3, MyExceptionSeverity severity)
    {
        std::cerr << "Exception: " << ss2 << " for " << ss3 << " **[ Issued by " << ss1 << " ]**";
        switch(severity) {
            case Fatal:
            case FATAL:
                std::cerr << " - FATAL ERROR\n" << std::endl;
                abort();
                break;
            case Warning:
            case WARNING:
                std::cerr << " - WARNING" << std::endl;
                break;
        }
    }
 
};

//================================================================

class CheckAlpha
{
public:
    bool operator()(char c) {
        if(isalpha(c) && c != 'e' && c != 'E') {
            std::string ss(&c);
            Exception("Invalid conversion of alpha character ("+ss+") to number",Fatal);
        }
        return isalpha(c);
    }
};
/*
//================================================================

template <typename T>
const T operator+(const T& lhs, const T& rhs) {
    return T(lhs) += rhs;
}

//================================================================

template <typename T>
const T operator-(const T& lhs, const T& rhs) {
    return T(lhs) -= rhs;
}

//================================================================
*/
    
    


#endif