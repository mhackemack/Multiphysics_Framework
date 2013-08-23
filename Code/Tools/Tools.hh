//
//
//
//
//
//
//

#ifndef tools_hh_
#define tools_hh_ 1

#include "MyTools.hh"
#include "MyHelpers.hh"

//==================================================================================

namespace Dimensions {
    // for namespace "using" declarations
    enum dimensions { dx = 0, dy = 1, dz = 2 };
    
    // for class inheritance
    class Dimensions_ {
    public:
        Dimensions_()
        : dx(0),dy(1),dz(2) { }
        
        ~Dimensions_() { }
        
    protected:
        uint32 dx, dy, dz;
    };
};

//==================================================================================

namespace CoordinateSystem
{
    enum CoordSystem { CARTESIAN, POLAR, SPHERICAL };
    
    inline std::string CoordSystemName(CoordSystem coordSystem)
    {
        switch (coordSystem) {
            case CARTESIAN:
                return std::string("Cartesian");
                break;
            case POLAR:
                return std::string("Polar");
                break;
            case SPHERICAL:
                return std::string("Spherical");
                break;
            default:
                break;
        }
        Exception("Coordinate System Name","Undefined coordinate system",Warning);
        return std::string("Undefined");
    }
};

//==================================================================================

template <typename T>
class RealNumber
{
public:
    //RealNumber() { }
    //~RealNumber() { }
    
    void operator()(T* value) {
        if(isinf(value) || isnan(value)) { *value = static_cast<T>(0); }
    }
    
private:
    inline bool isinf(T* value) {
        return (std::numeric_limits<T>::has_infinity &&
                *value == std::numeric_limits<T>::infinity()) ? true : false;
    }
    
    inline bool isnan(T* value) { return (*value != *value) ? true : false; }
    
};

//==================================================================================

using CoordinateSystem::CARTESIAN;
using CoordinateSystem::POLAR;
using CoordinateSystem::SPHERICAL;
using CoordinateSystem::CoordSystemName;
typedef CoordinateSystem::CoordSystem CoordSystem;

//==================================================================================
//      For 3-dimensional vectors (position,velocity,acceleration,etc.)
template <typename T>
class ThreeVector : public Dimensions::Dimensions_
{
public:
    enum CoordSystemChange {
        CARTESIAN_TO_POLAR, CARTESIAN_TO_SPERICAL,
        POLAR_TO_CARTESIAN, POLAR_TO_SPERICAL,
        SPHERICAL_TO_CARTESIAN, SPHERICAL_TO_POLAR
    };
    
    // ==============================================================================
    // Constructors
    // ==============================================================================    
    // No parameters
    inline
    ThreeVector():x_(0),y_(0),z_(0),coordSystem(CARTESIAN) {}
    
    // 1 parameter
    inline
    ThreeVector(T num)
    :   x_((num!=num)?0:num),y_((num!=num)?0:num),z_((num!=num)?0:num),coordSystem(CARTESIAN) {}
  
    // 3 parameters
    inline
    ThreeVector(T xn, T yn, T zn)
    :   x_(xn),y_(yn),z_(zn),coordSystem(CARTESIAN) {}
    
    // Direction vector
    inline
    ThreeVector(uint32 tmpdim, T value)
    :   x_((tmpdim==dx)?value:0),y_((tmpdim==dy)?value:0),
        z_((tmpdim==dz)?value:0),coordSystem(CARTESIAN) {}
    
    // Projected Vector -> for (x,y,z) vector projected into 2 dimensions
    inline
    ThreeVector(uint32 tmpdim, ThreeVector tmp)
    :   x_((tmpdim==dx)?0:tmp.x()),y_((tmpdim==dy)?0:tmp.y()),
        z_((tmpdim==dz)?0:tmp.z()),coordSystem(CARTESIAN) {}
    
    // 3 parameters + Coordinate System
    inline
    ThreeVector(T xn, T yn, T zn, CoordSystem coordSys)
    :   x_(xn),y_(yn),z_(zn),coordSystem(coordSys) {}
    
    // std::vector parameters
    ThreeVector(std::vector<T> vec)
    : x_(0),y_(0),z_(0),coordSystem(CARTESIAN)
    {
        switch (vec.size()) {
            case 0:
                break;
            case 1:
                x_ = vec.at(0);
                y_ = vec.at(0);
                z_ = vec.at(0);
                break;
            case 2:
                Exception(std::string(typeid(*this).name()),
                          "Initialization of ThreeVector with vector of size equal to 2",Fatal);
                break;
            case 3:
                x_ = vec.at(0);
                y_ = vec.at(1);
                z_ = vec.at(2);
                break;
            default:
                std::cout << "\n" << std::endl;
                for(uint32 i = 3; i < vec.size(); ++i) {
                    std::cout << "\tThreeVector<" << std::string(typeid(T).name())
                    << "> :: unused parameter " << i << ": "
                    << vec.at(i) << std::endl;
                }
                std::cout << "\n" << std::endl;
                Exception("ThreeVector",
                          "Excess parameters in initialization of ThreeVector with vector",Fatal);
                if(vec.size() <= 3) { break; }
                break;
        }
    }
    
    ThreeVector(std::vector<std::string> svec)
    :   x_(0),y_(0),z_(0),coordSystem(CARTESIAN)
    {
        *this = ThreeVector(MyTools::s2n<T>(svec));
    }
    
    /*ThreeVector(const ThreeVector& rhs)
    :   x_(rhs.x_),y_(rhs.y_),z_(rhs.z_),coordSystem(rhs.coordSystem)
    { }*/
    
    
    // ==============================================================================
    // Operators
    // ==============================================================================
    // Assignment
    inline ThreeVector& operator=(const ThreeVector &rhs) {
        if (this != &rhs ) {
            x_ = rhs.x_; y_ = rhs.y_; z_ = rhs.z_; coordSystem = rhs.coordSystem;
        }
        return *this;
    }
    // ==============================================================================
    // Equivalence
    inline bool operator==(const ThreeVector &rhs) const 
    {   if(this->x_ == rhs.x_ &&
           this->y_ == rhs.y_ &&
           this->z_ == rhs.z_ &&
           this->coordSystem == rhs.coordSystem) { return true; }
        else { return false; }
    }
    
    inline bool operator!=(const ThreeVector &rhs) const
    { return !(*this == rhs); }
    
    inline bool operator<=(const ThreeVector& rhs) {
        return (x_ <= rhs.x_ && y_ <= rhs.y_ && z_ <= rhs.z_) ? true : false;
    }

    inline bool operator<(const ThreeVector& rhs) {
        return (x_ < rhs.x_ && y_ < rhs.y_ && z_ < rhs.z_) ? true : false;
    }
    
    inline bool operator>=(const ThreeVector& rhs) {
        return (x_ >= rhs.x_ && y_ >= rhs.y_ && z_ >= rhs.z_) ? true : false;
    }
    
    inline bool operator>(const ThreeVector& rhs) {
        return (x_ > rhs.x_ && y_ > rhs.y_ && z_ > rhs.z_) ? true : false;
    }
    
    // ==================== ThreeVector Operators ====================
    // Addition of ThreeVector
    ThreeVector& operator+=(const ThreeVector& n)
    {
        std::auto_ptr<ThreeVector<T>> real(new ThreeVector<T>(this->convert(CARTESIAN)));
        std::auto_ptr<ThreeVector<T>> temp(new ThreeVector<T>((&const_cast<ThreeVector&>(n))->convert(CARTESIAN)));
        real->x_ += temp->x();
        real->y_ += temp->y();
        real->z_ += temp->z();
        return *this = real->convert(this->coord_system());
    }
    
    const ThreeVector operator+(const ThreeVector &rhs) const
    { return ThreeVector(*this) += rhs; }
    // ===============================================
    // Multiplication of ThreeVector
    ThreeVector& operator*=(const ThreeVector& n)
    {
        std::auto_ptr<ThreeVector<T>> real(new ThreeVector<T>(this->convert(CARTESIAN)));
        std::auto_ptr<ThreeVector<T>> temp(new ThreeVector<T>((&const_cast<ThreeVector&>(n))->convert(CARTESIAN)));
        real->x_ *= temp->x();
        real->y_ *= temp->y();
        real->z_ *= temp->z();
        return *this = real->convert(this->coord_system());
    }
    
    const ThreeVector operator*(const ThreeVector& rhs) const
    { return ThreeVector(*this) *= rhs; }
    // ===============================================
    // Subtraction of ThreeVector
    ThreeVector& operator-=(const ThreeVector& n)
    {
        std::auto_ptr<ThreeVector<T>> real(new ThreeVector<T>(this->convert(CARTESIAN)));
        std::auto_ptr<ThreeVector<T>> temp(new ThreeVector<T>((&const_cast<ThreeVector&>(n))->convert(CARTESIAN)));
        real->x_ -= temp->x();
        real->y_ -= temp->y();
        real->z_ -= temp->z();
        //*this = real->convert(this->coord_system());
        return *this = real->convert(this->coord_system());
    }
    
    const ThreeVector operator-(const ThreeVector& rhs) const
    { return ThreeVector(*this) -= rhs; }
    // ===============================================
    // Division of ThreeVector
    ThreeVector& operator/=(ThreeVector n)
    {
        std::auto_ptr<ThreeVector<T>> real(new ThreeVector<T>(this->convert(CARTESIAN)));
        std::auto_ptr<ThreeVector<T>> temp(new ThreeVector<T>((&const_cast<ThreeVector&>(n))->convert(CARTESIAN)));
        real->x_ /= temp->x();
        real->y_ /= temp->y();
        real->z_ /= temp->z();
        return *this = real->convert(this->coord_system());
    }
    
    const ThreeVector operator/(const ThreeVector &rhs) const
    { return ThreeVector(*this) /= rhs; }

    T& operator()(uint32 d_)
    {
        if(d_==0) { return x_; }
        else if(d_==1) { return y_; }
        else if(d_==2) { return z_; }

        Exception("ThreeVector","Operator() - out of range (>= 3)",MyTools::n2s(d_),Fatal);
        return x_;
    }

    // ===============================================
    // Output Stream operator
    friend std::ostream& operator<<(std::ostream& out,
                                    ThreeVector tv)
    { tv.print(", ",out); return out; }    

private:
    inline void check(T*);
    inline void check(const T*);
    
public:
    // Get functions
    inline T& x() { check(&x_); return x_; }
    inline T& y() { check(&y_); return y_; }
    inline T& z() { check(&z_); return z_; }
    // Get functions
    inline T cx() const { /*check(const_cast<T*>(&x_));*/ return x_; }
    inline T cy() const { /*check(const_cast<T*>(&y_));*/ return y_; }
    inline T cz() const { /*check(const_cast<T*>(&z_));*/ return z_; }
    // Get functions
    inline T tx() { check(&x_); return x_; }
    inline T ty() { check(&y_); return y_; }
    inline T tz() { check(&z_); return z_; }
    
    inline T magnitude() { return sqrt( x_*x() + y_*y() + z_*z() ); }
    inline T cmagnitude() const { return sqrt( x_*cx() + y_*cy() + z_*cz() ); }
    
    // Get value for dimension
    inline const T& dim(int d) const
    {
        switch(d) {
            case 0:
                return x_;
                break;
            case 1:
                return y_;
                break;
            case 2:
                return z_;
                break;
        }
        std::cout << "Passing dimension " << d << std::endl;
        Exception("ThreeVector","Bad dimension value passed to dim(int)",MyTools::n2s(d),Fatal);
        return x_;
    }
    
    // for printing as std::string
    inline std::string asString(T val) {
        check(&val);
        std::string delim = "";
        if(val == (-0)) { val = 0; }
        if(magnitude() > 0 && (fabs(val)/magnitude()) < 1.e-10) { val = 0.; }
        if(!(val<0) && 0 != MyTools::n2s(val).find("-")) { delim += " "; }
        return delim + MyTools::n2s(val);
    }
    
    // Get the coordinate system
    inline CoordSystem coord_system() { return coordSystem; }
    
    // Getting value as string
    inline std::string sx() { return asString(x_); }
    inline std::string sy() { return asString(y_); }
    inline std::string sz() { return asString(z_); }

    // Normalize the vector to 1 (permanent)
    inline void normalize() {
        T mag = magnitude();
        bool d_by_0 = DivideByZero(mag);
        if(d_by_0) { return; }
        x_ /= mag;
        y_ /= mag;
        z_ /= mag;
    }
    
    // For projection
    inline ThreeVector nullify(uint32 tmpdim) {
        switch(tmpdim) {
            case 0:
                return ThreeVector(0.,this->y_,this->z_);
                break;
            case 1:
                return ThreeVector(this->x_,0.,this->z_);
                break;
            case 2:
                return ThreeVector(this->x_,this->y_,0);
                break;
            default:
                Exception("ThreeVector","Nullify - passed bad index",Fatal);
                break;
        }
    }
    
    // For direction
    inline ThreeVector direction(uint32 tmpdim) {
        return ThreeVector(tmpdim,(this->dim(tmpdim)!=0)?this->dim(tmpdim):1);
    }
    
    // For normalization without permanent change
    inline ThreeVector normal() {
        return (magnitude() != 0) ? ThreeVector(x_/magnitude(),y_/magnitude(),z_/magnitude()) : ThreeVector(0.);
    }
    
    // For normalization without permanent change
    inline const ThreeVector cnormal() const {
        return (cmagnitude() != 0) ? ThreeVector(cx()/cmagnitude(),cy()/cmagnitude(),cz()/cmagnitude()) : ThreeVector(0.);
    }
    
    // Get absolute
    inline ThreeVector absolute() {
        return ThreeVector(((x_<0)?-x_:x_),
                           ((y_<0)?-y_:y_),
                           ((z_<0)?-z_:z_));
    }
    
    // Get an inverse
    inline ThreeVector inverse() {
        return ThreeVector(-1.*x_,-1.*y_,-1.*z_,coordSystem);
    }
    
    // Get an inverse around a value
    inline ThreeVector inverse(const T& val)
    {
        return ThreeVector( (x_ >= 0.) ? (x_ - val) : ( val - x_ ),
                            (y_ >= 0.) ? (y_ - val) : ( val - y_ ),
                            (z_ >= 0.) ? (z_ - val) : ( val - z_ ) );
    }
    
    // Invert
    inline void invert() { x_ *= -1.0; y_ *= -1.0; z_ *= -1.0; }
    // Invert around a value
    inline void invert(const T& val)
    {
        x_ = (x_ >= 0.) ? (x_ - val) : ( val - x_ );
        y_ = (y_ >= 0.) ? (y_ - val) : ( val - y_ );
        z_ = (z_ >= 0.) ? (z_ - val) : ( val - z_ );
    }
    
    // Flip the three vector
    inline void flip() { std::swap(x_,z_); }
    // Flip three vector
    inline ThreeVector flipped() { return ThreeVector(z_,y_,x_,coordSystem); }
    
    // Check if divide by zero
    inline bool DivideByZero(T n) {
        return (0 == n) ? true : false;
    }
    
    // Get Dimension of value
    inline uint32 GetDimension(T val) {
        if( fabs(val-x_) < 1.0e-6 ) {
            return 0;
        } else if( fabs(val-y_) < 1.0e-6 ) {
            return 1;
        } else if( fabs(val-z_) < 1.0e-6 ) {
            return 2;
        } else {
            Exception("ThreeVector","Could not find dimension for requested value",Fatal);
        }
    }
    
    inline ThreeVector units(T factor)
    {
        return ThreeVector(x_/factor,y_/factor,z_/factor);
    }
    
    // Convert to another coordinate system
    inline ThreeVector convert(CoordSystem);
    
    inline ThreeVector as_cartesian() { return convert(CARTESIAN); }
    inline ThreeVector as_polar() { return convert(POLAR); }
    inline ThreeVector as_spherical() { return convert(SPHERICAL); }
    
    // Set functions
    inline void SetX(T v) { x_ = v; }
    inline void SetY(T v) { y_ = v; }
    inline void SetZ(T v) { z_ = v; }
    
    inline void SetDim(int d, T v) {
        if(d==dx) { x_ = v;; }
        else if(d==dy) { y_ = v; }
        else if(d==dz) { z_ = v; }
        else {
            Exception(std::string(typeid(*this).name()),
                      "Invalid dimension - "+MyTools::n2s(d),Warning);
        }
    }
    
    inline void set(int d, T v) {
        if(d==dx) { x_ = v;; }
        else if(d==dy) { y_ = v; }
        else if(d==dz) { z_ = v; }
        else {
            Exception(std::string(typeid(*this).name()),
                      "Invalid dimension - "+MyTools::n2s(d),Warning);
        }
    }
    
    inline void add(int d, T v) {
        if(d == dx || d == 0) { x_ += v; }
        else if(d == dy || d == 1) { y_ += v; }
        else if(d == dz || d == 2) { z_ += v; }
        else {
            Exception(std::string(typeid(*this).name()),
                      "Invalid dimension - "+MyTools::n2s(d),Warning);
        }
    }
    
    inline void multiply(int d, T v) {
        if(d == dx || d == 0) { x_ *= v; }
        else if(d == dy || d == 1) { y_ *= v; }
        else if(d == dz || d == 2) { z_ *= v; }
        else {
            Exception(std::string(typeid(*this).name()),
                      "Invalid dimension - "+MyTools::n2s(d),Warning);
        }
    }

    inline ThreeVector add(T v) {
        return ThreeVector(x_+v,y_+v,z_+v);
    }
    
    inline ThreeVector multiply(T v) {
        return ThreeVector(x_*v,y_*v,z_*v);
    }
    
    inline ThreeVector radians() {
        return *this;
    }
    
    ThreeVector degrees() {
        auto radtodeg = [](T value) { return value * 180./M_PI; };
        switch (coordSystem) {
            case CARTESIAN:
                return ThreeVector(radtodeg(this->x()),radtodeg(this->y()),
                                   radtodeg(this->z()));
                break;
            case POLAR:
                return ThreeVector(this->x(),radtodeg(this->y()),this->z());
                break;
            case SPHERICAL:
                return ThreeVector(this->x(),radtodeg(this->y()),radtodeg(this->z()));
                break;
            default:
                break;
        }
        return *this;
    }
    
    void Set2Radians()
    {
        auto deg2rad = [](T value) { return value * M_PI/180.; };
        x_ = deg2rad(x_);
        y_ = deg2rad(y_);
        z_ = deg2rad(z_);
    }
    
    inline T min() { return (x_ < y_) ? ((x_ < z_)?x_:z_):((y_ < z_)?y_:z_); }
    inline T max() { return (x_ > y_) ? ((x_ > z_)?x_:z_):((y_ > z_)?y_:z_); }
    
    inline ThreeVector GetMaximums(ThreeVector compare) {
        return ThreeVector((x_ > compare.x_) ? x_ : compare.x_,
                           (y_ > compare.y_) ? y_ : compare.y_,
                           (z_ > compare.z_) ? z_ : compare.z_);
    }
    
    inline ThreeVector GetMinimums(ThreeVector compare) {
        return ThreeVector((x_ < compare.x_) ? x_ : compare.x_,
                           (y_ < compare.y_) ? y_ : compare.y_,
                           (z_ < compare.z_) ? z_ : compare.z_);
    }
    
    void print(std::ostream& o) {
        o << "\t" << x() << " " << y() << " " << z() << std::endl;
    }
    
    void print(std::ostream& o, std::string delim) {
        o << " <" << x() << delim << y() << delim << z() << "> " << std::endl;
    }

    void print(std::string delim, std::ostream& o) {
        std::stringstream ss;
        ss << "<" << sx() << delim << sy() << delim << sz() << ">";
        o << ss.str();
    }
    
    inline std::string as_string() {
        std::stringstream ss;
        ss << "<" << sx() << ", " << sy() << ", " << sz() << ">";
        return ss.str();
    }
    
    inline void SetCoordinateSystem(CoordSystem coordSys) { coordSystem = coordSys; }
    inline CoordSystem GetCoordinateSystem() { return coordSystem; }
    
private:
    T x_;
    T y_;
    T z_;
    CoordSystem coordSystem;
    RealNumber<T> realNumberCheck;

};

//=================================================================================================================================
//
//=================================================================================================================================
template <>
class RealNumber<ThreeVector<double>>
{
public:
    
    void operator()(ThreeVector<double>* v3) {
        if(isinf(v3->x()) || isnan(v3->x()) ) { v3->SetX(static_cast<double>(0.)); }
        if(isinf(v3->y()) || isnan(v3->y()) ) { v3->SetY(static_cast<double>(0.)); }
        if(isinf(v3->z()) || isnan(v3->z()) ) { v3->SetZ(static_cast<double>(0.)); }
    }
    
private:
    inline bool isinf(double value) {
        return (std::numeric_limits<double>::has_infinity &&
                fabs(value) == std::numeric_limits<double>::infinity()) ? true : false;
    }
    
    inline bool isnan(double value) { return (value != value) ? true : false; }
};

//=================================================================================================================================
//
//=================================================================================================================================
template <>
class RealNumber<ThreeVector<float>>
{
public:
    
    void operator()(ThreeVector<float>* v3) {
        if(isinf(v3->x()) || isnan(v3->x()) ) { v3->SetX(static_cast<float>(0.)); }
        if(isinf(v3->y()) || isnan(v3->y()) ) { v3->SetY(static_cast<float>(0.)); }
        if(isinf(v3->z()) || isnan(v3->z()) ) { v3->SetZ(static_cast<float>(0.)); }
        
    }
private:
    inline bool isinf(float value) {
        return (std::numeric_limits<float>::has_infinity &&
                fabs(value) == std::numeric_limits<float>::infinity()) ? true : false;
    }
    
    inline bool isnan(float value) { return (value != value) ? true : false; }
};

//=================================================================================================================================

template <typename T>
void ThreeVector<T>::check(T* val) {
    realNumberCheck(val);
}

template <typename T>
void ThreeVector<T>::check(const T* val) {
    realNumberCheck(const_cast<T*>(val));
}
//=================================================================================================================================

template <typename T>
ThreeVector<T> ThreeVector<T>::convert(CoordSystem newSystem)
{
    ThreeVector<T> tmp = *this;
    
    if(!(newSystem == coordSystem)) {
        int systemcase = 0;
        switch (newSystem) {
            case CARTESIAN:
                systemcase = (coordSystem == POLAR) ? POLAR_TO_CARTESIAN : SPHERICAL_TO_CARTESIAN;
                break;
            case SPHERICAL:
                systemcase = (coordSystem == POLAR) ? POLAR_TO_SPERICAL : CARTESIAN_TO_SPERICAL;
                break;
            case POLAR:
                systemcase = (coordSystem == CARTESIAN) ? CARTESIAN_TO_POLAR : SPHERICAL_TO_POLAR;
            default:
                break;
        }
        
        tmp.SetCoordinateSystem(newSystem);
        
        auto magnitude = [](T lx, T ly) {
            return sqrt(lx*lx + ly*ly);
        };
        switch (systemcase) {
            case POLAR_TO_CARTESIAN:
                return ThreeVector<T>(x_*cos(y_),x_*sin(y_),z_,newSystem);
                break;
            case SPHERICAL_TO_CARTESIAN:
                return ThreeVector<T>(x_*cos(y_)*sin(z_),x_*sin(y_)*sin(z_),x_*cos(z_),newSystem);
                break;
            case POLAR_TO_SPERICAL:
                return ThreeVector<T>(magnitude(x_,z_),y_,atan2(x_,z_),newSystem);
                break;
            case CARTESIAN_TO_SPERICAL:
                return ThreeVector<T>(this->magnitude(),atan2(y_,x_),acos(z_/this->magnitude()),newSystem);
                break;
            case CARTESIAN_TO_POLAR:
                return ThreeVector<T>(magnitude(x_,y_),atan2(y_,x_),z_,newSystem);
                break;
            case SPHERICAL_TO_POLAR:
                return ThreeVector<T>(x_*sin(z_),y_,x_*cos(z_),newSystem);
                break;
            default:
                break;
        }
    }
    
    return *this;
}

const ThreeVector<double> funitVector = ThreeVector<double>(1.0);

namespace ThreeVectorLinearAlgebra
{
    using Dimensions::dx;
    using Dimensions::dy;
    using Dimensions::dz;
    
    //===============================================================================
    template <typename T>
    inline T Magnitude(const ThreeVector<T>& lhs, const ThreeVector<T>& rhs) {
        return sqrt(pow(lhs.cx()-rhs.cx(),2)+pow(lhs.cy()-rhs.cy(),2)+pow(lhs.cz()-rhs.cz(),2));
    }
    //===============================================================================
    template <typename T>
    inline T Magnitude(const ThreeVector<T>& lhs) {
        return sqrt(pow(lhs.cx(),2)+pow(lhs.cy(),2)+pow(lhs.cz(),2));
    }
    //===============================================================================
    template <typename T>
    inline T MagnitudeNC(ThreeVector<T> lhs, ThreeVector<T> rhs) {
        return sqrt(pow(lhs.x()-rhs.x(),2)+pow(lhs.y()-rhs.y(),2)+pow(lhs.z()-rhs.z(),2));
    }
    //===============================================================================
    template <typename T>
    inline T dotproduct(const ThreeVector<T>& left, const ThreeVector<T>& right) {
        return ((left.cx()*right.cx()) + (left.cy()*right.cy()) + (left.cz()*right.cz()));
    }
    //===============================================================================
    template <typename T>
    inline ThreeVector<T> crossproduct(const ThreeVector<T>& top, const ThreeVector<T>& bot) {
        return ThreeVector<T>(  (top.cy()*bot.cz()) - (bot.cy()*top.cz()),
                              -((top.cx()*bot.cz()) - (bot.cx()*top.cz())),
                                (top.cx()*bot.cy()) - (bot.cx()*top.cy()) );
    }
    //===============================================================================
    template <typename T>
    inline ThreeVector<T> projection(const ThreeVector<T>& A, const ThreeVector<T>& ontoB) {
        auto proj = [] (const ThreeVector<T>& a_, const ThreeVector<T>& b_) {
            return b_*(dotproduct(a_,b_)/pow(b_.magnitude(),2));
        };
        return ThreeVector<T>(proj(A,ontoB));
    }
    //===============================================================================
    template <typename T>
    inline double GetAngle(const ThreeVector<T>& v1, const ThreeVector<T>& v2) {
        return acos( (dotproduct(v1,v2) / (Magnitude(v1)*Magnitude(v2)) ) );
    }
    //===============================================================================
    template <typename T>
    inline ThreeVector<T> rejection(ThreeVector<T> A, ThreeVector<T> ontoB) {
        auto proj = [] (ThreeVector<T>& a_, ThreeVector<T>& b_) {
            return b_*(dotproduct(a_,b_)/pow(b_.magnitude(),2));
        };
        ThreeVector<T> tmp(proj(A,ontoB));
        return ThreeVector<T>(A.x()-tmp.x(),A.y()-tmp.y(),A.z()-tmp.z());
    }
    //===============================================================================
    template <typename T>
    inline ThreeVector<T> ComputeEulerAngles(ThreeVector<T> moved) {
        
        //if(moved == ThreeVector<T>(0.)) { moved = ThreeVector<T>(0.,0.,1.); }
        
        /*
        std::cout << "\n\tA - Projection onto xy-plane is using " << moved.direction(dx) << " and " << moved.nulliy_(dz) << std::endl;
        std::cout << "\tB - Projection onto yz-plane is using " << moved.nulliy_(dz) << " and " << moved << std::endl;
        std::cout << "\tG - Projection onto xz-plane is using " << moved.direction(dz) << " and " << moved.nulliy_(dy) << std::endl;
        std::cout << "\n" << std::endl;
        */
        
        auto GreaterThanPi = [] (T a, T b) -> T {
            uint32 fac = 1;
            // Distinguish between quadrants 1,3 or 2,4
            if(a*-b - b*-a > 0) { fac = 2; } // if {} -> quadrant 2 or 4
                                             // else { fac = 1; } -> quadrant 1 or 3
                                             
            // Distinguish between quadrant 1,2 or 3,4
            if(fac%2 == 1) {
                if(a+b < 0) { fac += 2; }   // if {} -> quadrant = 3
                                            // else {} -> quadrant = 1
            } else {
                if(b-a > 0) { fac += 2; }   // if {} -> quadrant = 4
                                            // else {} -> quadrant = 2
            }
            
            // Return appropriate value based on quadrant
            switch(fac) {
                case 1:
                case 2:
                    return 0;
                    break;
                case 3:
                case 4:
                    return -2*M_PI;
                    break;
            }
        };
        
        return ThreeVector<T>(acos(projection(moved.nullify(dz),moved.direction(dx)).x()/moved.nullify(dz).magnitude())+
                                GreaterThanPi(moved.x(),moved.y()),
                              acos(projection(moved,moved.nullify(dz)).y()/moved.nullify(dz).magnitude()),
                              //0.,
                              acos(projection(moved.nullify(dy),moved.direction(dz).absolute()).z()/moved.nullify(dy).magnitude())+
                                GreaterThanPi(moved.z(),moved.y()));
    }
    //===============================================================================
    template <typename T>
    inline ThreeVector<T> ComputeEulerAngles(const ThreeVector<T>& a, const ThreeVector<T>& b) {
        ThreeVector<T> a_normal = const_cast<ThreeVector<T>&>(a).normal();
        ThreeVector<T> b_normal = const_cast<ThreeVector<T>&>(b).normal();
        std::cout << "A normal = \t" << a_normal << std::endl;
        std::cout << "B normal = \t" << b_normal << std::endl;
        auto delta_x = b_normal.x() - a_normal.x();
        auto delta_y = b_normal.y() - a_normal.y();
        auto delta_z = b_normal.z() - a_normal.z();
        std::cout << "Delta X = " << delta_x << std::endl;
        std::cout << "Delta Y = " << delta_y << std::endl;
        std::cout << "Delta Z = " << delta_z << std::endl;
        //return ThreeVector<T>(atan2(delta_z,1.),
        //                      atan2(delta_z,1.),
        //                      atan2(delta_x,delta_y));
        return ThreeVector<T>(asin(delta_z),
                              asin(delta_x),
                              asin(delta_y));
        //return ThreeVector<T>(asin(b_normal.z()-a_normal.z()),
        //                      asin(b_normal.z()-a_normal.z()),
        //                      asin(b_normal.x)
        //auto aEuler = ComputeEulerAngles(a);
        //auto bEuler = ComputeEulerAngles(b);
        //return bEuler - aEuler;
    }
    //===============================================================================

};

using ThreeVectorLinearAlgebra::Magnitude;
using ThreeVectorLinearAlgebra::MagnitudeNC;
using ThreeVectorLinearAlgebra::dotproduct;
using ThreeVectorLinearAlgebra::crossproduct;
using ThreeVectorLinearAlgebra::projection;
using ThreeVectorLinearAlgebra::rejection;
using ThreeVectorLinearAlgebra::ComputeEulerAngles;
using ThreeVectorLinearAlgebra::GetAngle;

typedef ThreeVector<double> dThreeVector;
typedef ThreeVector<double> d3Vector;
typedef ThreeVector<float> fThreeVector;
typedef ThreeVector<float> f3Vector;


//==================================================================================

namespace tools
{
    //using namespace MyTools;
    
    using MyTools::MakeLower;
    using MyTools::MakeUpper;
    using MyTools::lower;
    using MyTools::upper;
    using MyTools::checkpos;
    using MyTools::checkposic;
    using MyTools::s2i;
    using MyTools::s2d;
    using MyTools::s2n;
    using MyTools::n2s;
    using MyTools::RemoveElement;
    using MyTools::Print;
    using MyTools::subvector;
    using MyTools::checkforID;
    using MyTools::removepos;
    using MyTools::Purge;
    using MyTools::RemoveDuplicates;
    using MyTools::count;
    using MyTools::IsDuplicate;
    using MyTools::Scale;
    using MyTools::mean;
    using MyTools::StdDev;
    using MyTools::RelativeError;
    using MyTools::ExitWithoutParams;
    using MyTools::OpenFile;
    using MyTools::RemoveComment;
    using MyTools::GetLineAndDelimit;
    using MyTools::ThrowError;
    using MyTools::delimit;
    using MyTools::undelimit;
    using MyTools::Capitalize;
    using MyTools::Uncapitalize;
    using MyTools::RemoveWhitespace;
    using MyTools::strcmpic;
    using MyTools::xmlopen;
    using MyTools::xmlclose;
    using MyTools::tablevel;
    using MyTools::AddValue;
    using MyTools::GetMin;
    using MyTools::GetMax;
    using MyTools::Sum;
    using MyTools::GenerateXMLblock;
    using MyTools::AddLine;
    using MyTools::lineXML;
    using MyTools::GetBoolFromString;
    using MyTools::RequireSize;
    using MyTools::CheckSize;
    using MyTools::GetParam;
    using MyTools::CheckForParam;
    using MyTools::ConvertVector;
    using MyTools::ConvertStringVector;
    using MyTools::SetDoublePrecision;
    using MyTools::GetIteratorPosition;
    using MyTools::GetIteratorIndex;
    using MyTools::AddDirectory;
    using MyTools::MakeDirectoryString;
    using MyTools::interpolate;
    using MyTools::loginterpolate;
    using MyTools::ExplicitEuler;
    using MyTools::ImplicitEuler;
    using MyTools::RungeKutta2nd;
    using MyTools::RungeKutta4th;
    using MyTools::RungeKutta4thUnits;
    //using MyTools::factorial;
    
};

namespace tools
{
    inline void print_line_break(std::ostream&,uint32);
    inline void print_line_break(std::stringstream&,uint32);
    
}

inline void tools::print_line_break(std::ostream& out, uint32 ncount)
{
    for(auto i = 0; i < ncount; ++i) { out << std::endl; }
}

inline void tools::print_line_break(std::stringstream& out, uint32 ncount)
{
    for(auto i = 0; i < ncount; ++i) { out << std::endl; }
}

#endif
