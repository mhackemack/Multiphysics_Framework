//
//
//
//
//
//

#ifndef rotationmatrix_hh_
#define rotationmatrix_hh_

#include "Tools.hh"

namespace GLOBAL {

    
    template <typename T> ThreeVector<ThreeVector<T>> identity_matrix()
    {
        return ThreeVector<ThreeVector<T>>(ThreeVector<T>(1.,0.,0.),
                                           ThreeVector<T>(0.,1.,0.),
                                           ThreeVector<T>(0.,0.,1.));
    }
    //typedef ThreeVector<ThreeVector<T>> RotateMatrix;
    
    const ThreeVector<double> Zero3Vector(0.);
}

using namespace GLOBAL;

template <typename T>
class RotationMatrix
{
public:
    typedef ThreeVector<ThreeVector<T>> RotateMatrix;

    // constructors
    inline RotationMatrix();
    inline RotationMatrix(std::vector<std::vector<T>>);
    //RotationMatrix(T,T);
    inline RotationMatrix(RotateMatrix,ThreeVector<T>);
    
    inline T operator()(uint32 i, uint32 j) { return current_matrix_(i)(j); }
    
    // destructor
    ~RotationMatrix() { }
    
    // given axis, angle, and ThreeVector to rotate
    inline RotateMatrix get_origin_rotation_matrix(const ThreeVector<T>&);
    
    // =================================================================================================
    // rotation with euler angles
    inline ThreeVector<T> get_rotation(const ThreeVector<T>&, const ThreeVector<T>&,
                                            const ThreeVector<T>& = ThreeVector<T>(0,0,0));
    inline ThreeVector<T> get_rotation(T, T, const ThreeVector<T>&);
    inline ThreeVector<T> get_inverse_rotation(const ThreeVector<T>&, const ThreeVector<T>&,
                                                    const ThreeVector<T>& = ThreeVector<T>(0,0,0));
    inline ThreeVector<T> get_reverse_rotation(ThreeVector<T>&);
    inline ThreeVector<T> get_axis_rotation(T, const ThreeVector<T>&,
                                                 const ThreeVector<T>&,
                                                 ThreeVector<T> position =
                                                 const_cast<ThreeVector<T>&>(Zero3Vector));
    inline RotateMatrix get_unique_axis_rotation(T theta, const ThreeVector<T>& axis_);
    inline ThreeVector<T> get_rotation(const RotateMatrix&, const ThreeVector<T>&, const ThreeVector<T>& position =
                                       ThreeVector<T>(0.));
    // =================================================================================================

    // =================================================================================================
    // rotation with vectors
    inline ThreeVector<T> get_vector_rotation(ThreeVector<T>&, ThreeVector<T>&);
    inline RotateMatrix MakeRotationMatrix(ThreeVector<T>,ThreeVector<T>);
    inline void convert2angles(ThreeVector<T>&);
    // =================================================================================================
    
    // =================================================================================================
    // multiply 3x3 and 3x1 matrix
    inline ThreeVector<T> multiply_matrix(const RotateMatrix&, const ThreeVector<T>&);
    // multiply 3x1 and 1x3 matrix
    inline RotateMatrix multiply_matrix(const ThreeVector<T>&, const ThreeVector<T>&);
    // get inverse of ThreeVector<T>
    inline RotateMatrix GetInverse(const RotateMatrix&);
    // get transpose of ThreeVector<T>
    inline RotateMatrix GetTranspose(const RotateMatrix&);
    // get current rotation matrix
    inline RotateMatrix& GetCurrentMatrix() { return current_matrix_; }
    // set the current matrix
    inline void SetCurrentMatrix(const RotateMatrix& rot) { current_matrix_ = rot; }
    // get inverse of current rotation matrix    
    inline RotateMatrix GetCurrentMatrixInverse() { return GetInverse(current_matrix_); }
    // get transpose of current rotation matrix
    inline RotateMatrix GetCurrentMatrixTranspose() { return GetTranspose(current_matrix_); }
    // get the current position
    inline ThreeVector<T> GetCurrentPosition() { return current_position_; }
    // set the current position
    inline void SetCurrentPosition(const ThreeVector<T>& val) { current_position_ = val; }
    // get rotation matrix history (not including current)
    inline std::list<RotateMatrix>& GetRotationHistory() { return rotation_history_; }
    // get current angles
    inline ThreeVector<T>& GetCurrentAngles() { return current_angle_; }
    // get current rotation point
    inline ThreeVector<T>& GetCurrentPoint() { return current_position_; }
    // =================================================================================================


    // =================================================================================================
    // Retrieval w/o modifying history
    inline ThreeVector<T> retrieve_rotation(const ThreeVector<T>&, const ThreeVector<T>&);
    inline void apply_rotation(ThreeVector<T>&, const RotateMatrix&, const ThreeVector<T>&);
    // =================================================================================================

    
    // =================================================================================================
    // Histories for rotation and euler angles
    inline void pop_rotation_history();
    inline void push_rotation_history(const RotateMatrix&);
    inline void pop_angle_history();
    inline void push_angle_history(const ThreeVector<T>&);
    inline void pop_position_history();
    inline void push_position_history(const ThreeVector<T>&);
    inline ThreeVector<T> GetEulerAngles(const RotateMatrix&);
    inline ThreeVector<T> GetEulerAngles() { return GetEulerAngles(current_matrix_); }
    // =================================================================================================

    
    // =================================================================================================
    inline void print_rotation_history(std::ostream&);
    inline void print_rotation(const RotateMatrix&,std::ostream& out = std::cout);
    // =================================================================================================
    

    // =================================================================================================
    RotationMatrix& operator=(RotationMatrix& rhs) {
        if(this != &rhs) {
            current_matrix_ = rhs.current_matrix_;
            current_angle_ = rhs.current_angle_;
            current_position_ = rhs.current_position_;
            rotation_history_ = rhs.rotation_history_;
            angle_history_ = rhs.angle_history_;
            position_history_ = rhs.position_history_;
            verbose = rhs.verbose;
        }
        return *this;
    }
    // =================================================================================================
    friend std::ostream& operator<<(std::ostream& out, const RotationMatrix& rm) {
        rm.print_rotation(rm.current_matrix_,out);
        return out;
    }
    // =================================================================================================

private:
    RotateMatrix current_matrix_;
    ThreeVector<T> current_angle_;
    ThreeVector<T> current_position_;
    std::list<RotateMatrix> rotation_history_;
    std::list<ThreeVector<T>> angle_history_;
    std::list<ThreeVector<T>> position_history_;
    uint32 verbose;
    RotateMatrix cofactor;
    
public:
    RotationMatrix& operator=(const RotationMatrix& rhs)
    {
        if(this != &rhs) {
            current_matrix_ = rhs.current_matrix_;
            current_angle_ = rhs.current_angle_;
            current_position_ = rhs.current_position_;
            rotation_history_ = rhs.rotation_history_;
            angle_history_ = rhs.angle_history_;
            position_history_ = rhs.position_history_;
            verbose = rhs.verbose;
            cofactor = rhs.cofactor;
        }
        return *this;
    }
    
    bool operator==(const RotationMatrix& rhs)
    {
        return (current_matrix_ == rhs.current_matrix_ && current_angle_ == rhs.current_angle_ &&
                current_position_ == rhs.current_position_ && rotation_history_ == rhs.rotation_history_ &&
                angle_history_ == rhs.angle_history_ && position_history_ == rhs.position_history_) ? true : false;
    }
    
    bool operator!=(const RotationMatrix& rhs)
    {
        return !(*this == rhs);
    }
};


//=============================================================================================================================
//=============================================================================================================================
//=============================================================================================================================
//
//                  Implementations
//
//=============================================================================================================================
//=============================================================================================================================
//=============================================================================================================================
template <typename T>
inline
RotationMatrix<T>::RotationMatrix()
:   current_matrix_(identity_matrix<T>()),current_angle_(ThreeVector<T>(static_cast<T>(0))),
    current_position_(ThreeVector<T>(static_cast<T>(0))),verbose(0),
    cofactor(ThreeVector<ThreeVector<T>>(ThreeVector<T>(1.,-1.,1),
                                         ThreeVector<T>(-1.,1.,-1),
                                         ThreeVector<T>(1.,-1.,1)))
{
    rotation_history_.resize(1,identity_matrix<T>());
    angle_history_.resize(1);
    position_history_.resize(1);
}

//=============================================================================================================================
template <typename T>
inline
typename RotationMatrix<T>::RotateMatrix RotationMatrix<T>::get_origin_rotation_matrix(const ThreeVector<T>& angle_)
{    
    return RotateMatrix(ThreeVector<T>(cos(angle_.cy())*cos(angle_.cz()),
                                            (-cos(angle_.cx())*sin(angle_.cz())) + (sin(angle_.cx())*sin(angle_.cy())*cos(angle_.cz())),
                                            (sin(angle_.cx())*sin(angle_.cz())) + (cos(angle_.cx())*sin(angle_.cy())*cos(angle_.cz()))),
                        ThreeVector<T>(cos(angle_.cy())*sin(angle_.cz()),
                                            (cos(angle_.cx())*cos(angle_.cz())) + (sin(angle_.cx())*sin(angle_.cy())*sin(angle_.cz())),
                                            (-sin(angle_.cx())*cos(angle_.cz())) + (cos(angle_.cx())*sin(angle_.cy())*sin(angle_.cz()))),
                        ThreeVector<T>(-sin(angle_.cy()),
                                            sin(angle_.cx())*cos(angle_.cy()),
                                            cos(angle_.cx())*cos(angle_.cy())));
}

//=============================================================================================================================
template <typename T>
inline
typename RotationMatrix<T>::RotateMatrix RotationMatrix<T>::get_unique_axis_rotation(T theta, const ThreeVector<T>& axis_)
{
    ThreeVector<T> row1 = ThreeVector<T>(cos(theta)+(axis_.cx()*axis_.cx()*(1.-cos(theta))),
                                                   (axis_.cx()*axis_.cy()*(1.-cos(theta))) - (axis_.cz()*sin(theta)),
                                                   (axis_.cx()*axis_.cz()*(1.-cos(theta))) + (axis_.cy()*sin(theta)));
    
    ThreeVector<T> row2 = ThreeVector<T>((axis_.cy()*axis_.cx()*(1.-cos(theta))) + (axis_.cz()*sin(theta)),
                                                   cos(theta) + (axis_.cy()*axis_.cy()*(1.-cos(theta))),
                                                   (axis_.cy()*axis_.cz()*(1.-cos(theta))) - (axis_.cx()*sin(theta)));
    
    ThreeVector<T> row3 = ThreeVector<T>((axis_.cz()*axis_.cx()*(1.-cos(theta))) - (axis_.cy()*sin(theta)),
                                                   (axis_.cz()*axis_.cy()*(1.-cos(theta))) + (axis_.cx()*sin(theta)),
                                                   cos(theta) + (axis_.cz()*axis_.cz()*(1.-cos(theta))));
    return RotateMatrix(row1,row2,row3);
}

//=============================================================================================================================
template <typename T>
inline
ThreeVector<T> RotationMatrix<T>::get_rotation(const ThreeVector<T>& angle_,
                                               const ThreeVector<T>& subject,
                                               const ThreeVector<T>& position)
{
    ThreeVector<T> subject_origin = subject - position;
    RotateMatrix rotation = get_origin_rotation_matrix(angle_);
    
    push_rotation_history(rotation);
    push_angle_history(angle_);
    push_position_history(position);
    
    return multiply_matrix(rotation,subject_origin) + position;
}
//=============================================================================================================================
template <typename T>
inline
ThreeVector<T> RotationMatrix<T>::get_rotation(T phi, T theta, const ThreeVector<T>& momentum_direction)
{
    // phi is rotation about z
    // theta is incline angle
    ThreeVector<double> angles(theta,0.,phi);
    push_rotation_history(get_origin_rotation_matrix(angles));
    push_angle_history(angles);
    push_position_history(ThreeVector<T>(0.));
    
    if(momentum_direction.magnitude() > 1) {
        Exception("RotationMatrix","Getting angle rotation for unnormalized vector",WARNING);
    }
    
    return ThreeVector<T>((fabs(momentum_direction.z()) >= 1.) ? ( sin(theta)*cos(phi) ) :
                          ((sin(theta)*(momentum_direction.x()*momentum_direction.z()*cos(phi) -
                                        momentum_direction.y()*sin(phi)))/std::sqrt(1.-pow(momentum_direction.z(),2)) +
                           momentum_direction.x()*cos(theta)),
                          (fabs(momentum_direction.z()) >= 1.) ? ( sin(theta)*sin(phi) * ((momentum_direction.z() > 0.)
                                                                                          ? 1. : -1.) ) :
                          ((sin(theta)*(momentum_direction.y()*momentum_direction.z()*cos(phi) -
                                        momentum_direction.x()*sin(phi)))/std::sqrt(1.-pow(momentum_direction.z(),2)) +
                           momentum_direction.y()*cos(theta)),
                          (fabs(momentum_direction.z()) >= 1.) ? ( cos(theta) * ((momentum_direction.z() > 0.)
                                                                                 ? 1. : -1.) ) :
                          -(sqrt(1.-pow(momentum_direction.z(),2))*sin(theta)*cos(phi)) +
                          momentum_direction.z()*cos(theta)).normal();
}

//=============================================================================================================================
template <typename T>
inline ThreeVector<T> RotationMatrix<T>::get_rotation(const RotateMatrix& rotation, const ThreeVector<T>& subject,
                                                      const ThreeVector<T>& position)
{
    push_rotation_history(rotation);
    push_angle_history(GetEulerAngles(rotation));
    push_position_history(position);
    
    return multiply_matrix(rotation,subject-position) + position;
}
//=============================================================================================================================
template <typename T>
inline void RotationMatrix<T>::apply_rotation(ThreeVector<T>& subject, const RotateMatrix& rot, const ThreeVector<T>& position)
{
    subject = multiply_matrix(rot,subject-position) + position;
}
//=============================================================================================================================
template <typename T>
inline
ThreeVector<T> RotationMatrix<T>::get_axis_rotation(T angle_, const ThreeVector<T>& subject,
                                                      const ThreeVector<T>& axis_, ThreeVector<T> position)
{
    ThreeVector<T> subject_origin = subject - position;
    RotateMatrix rotation = get_unique_axis_rotation(angle_,axis_.cnormal());
    
    push_rotation_history(rotation);
    push_angle_history(angle_);
    push_position_history(position);
    
    return multiply_matrix(rotation,subject_origin) + position;
}
//=============================================================================================================================
template <typename T>
inline
ThreeVector<T> RotationMatrix<T>::get_inverse_rotation(const ThreeVector<T>& angle_,
                                                    const ThreeVector<T>& subject,
                                                    const ThreeVector<T>& position)
{
    ThreeVector<T> subject_origin = subject - position;
    RotateMatrix inv_rotation = GetInverse(get_origin_rotation_matrix(angle_));
    
    //print_rotation(inv_rotation);
    
    push_rotation_history(inv_rotation);
    push_angle_history(ThreeVector<T>(static_cast<T>(angle_.cz()),
                                      static_cast<T>(angle_.cy()),
                                      static_cast<T>(angle_.cx())));
    push_position_history(position);
    
    return multiply_matrix(inv_rotation,subject_origin) + position;
}

//=============================================================================================================================
template <typename T>
inline
ThreeVector<T> RotationMatrix<T>::retrieve_rotation(const ThreeVector<T>& subject, const ThreeVector<T>& position)
{
    return multiply_matrix(current_matrix_,subject-current_position_) + position;
}

//=============================================================================================================================
template <typename T>
inline
ThreeVector<T> RotationMatrix<T>::get_vector_rotation(ThreeVector<T>& v1,
                                                        ThreeVector<T>& v2)
{
    RotateMatrix rotation = MakeRotationMatrix(v1,v2);
    
    push_rotation_history(rotation);
    push_angle_history(multiply_matrix(rotation,v1));
    push_position_history(ThreeVector<T>(0));
    
    print_rotation(rotation);
    std::cout << " Angle is " << current_angle_ << std::endl;
    std::cout << " Current Position is " << current_position_ << std::endl;
    
    return multiply_matrix(rotation,v1);
}

//=============================================================================================================================
template <typename T>
inline
typename RotationMatrix<T>::RotateMatrix RotationMatrix<T>::MakeRotationMatrix(ThreeVector<T> X, ThreeVector<T> Y)
{
    X.normalize();
    Y.normalize();
    // make sure that we actually have two unique vectors.
    assert( X != Y );
        
    RotateMatrix M;
    M.x() = X;
    M.z() = crossproduct(X,Y).normal();
    M.y() = crossproduct(M.z(),X).normal();
    
    convert2angles(M(0));
    convert2angles(M(1));
    convert2angles(M(2));
    
    return M;
}

//=============================================================================================================================
template <typename T>
inline
void RotationMatrix<T>::convert2angles(ThreeVector<T>& v3)
{
    std::cout << std::endl;
    auto mag = [] (T x_, T y_) -> decltype(x_) {
        T answer = sqrt(pow(x_,2) + pow(y_,2));
        std::cout << "Magnitude of " << x_ << " and " << y_ << " is " << answer << std::endl;
        return ((answer == 0) ? 1 : answer);
    };
    
    T ax = v3.y()/mag(v3.y(),v3.z());
    T ay = v3.z()/mag(v3.x(),v3.z());
    T az = v3.x()/mag(v3.x(),v3.y());

    ThreeVector<T> angles(-acos(ax),acos(ay),-acos(az));
    
    std::cout << "Angles is " << angles << std::endl;
    std::cout << v3 << " converted to angles is " << std::flush;
    v3 = angles;
    std::cout << v3 << std::endl;
}

//=============================================================================================================================
template <typename T>
inline
ThreeVector<T> RotationMatrix<T>::get_reverse_rotation(ThreeVector<T>& subject)
{
    ThreeVector<T> subject_origin = subject - current_position_;
    RotateMatrix trans_rotation = GetTranspose(current_matrix_);
    pop_rotation_history();
    pop_angle_history();
    
    ThreeVector<T> rt = multiply_matrix(trans_rotation,subject_origin) + current_position_;
    pop_position_history();
    return rt;
}

//=============================================================================================================================
template <typename T>
inline
ThreeVector<T> RotationMatrix<T>::multiply_matrix(const RotateMatrix& mtx, const ThreeVector<T>& rhs)
{
    return ThreeVector<T>(mtx.dim(0).dim(0)*rhs.dim(0) + mtx.dim(0).dim(1)*rhs.dim(1) + mtx.dim(0).dim(2)*rhs.dim(2),
                               mtx.dim(1).dim(0)*rhs.dim(0) + mtx.dim(1).dim(1)*rhs.dim(1) + mtx.dim(1).dim(2)*rhs.dim(2),
                               mtx.dim(2).dim(0)*rhs.dim(0) + mtx.dim(2).dim(1)*rhs.dim(1) + mtx.dim(2).dim(2)*rhs.dim(2));
}

//=============================================================================================================================
template <typename T>
inline
typename RotationMatrix<T>::RotateMatrix RotationMatrix<T>::multiply_matrix(const ThreeVector<T>& a, const ThreeVector<T>& b)
{
    return ThreeVector<ThreeVector<T>>(ThreeVector<T>(a.dim(0)*b.dim(0),a.dim(0)*b.dim(1),a.dim(0)*b.dim(2)),
                                            ThreeVector<T>(a.dim(1)*b.dim(0),a.dim(1)*b.dim(1),a.dim(1)*b.dim(2)),
                                            ThreeVector<T>(a.dim(2)*b.dim(0),a.dim(2)*b.dim(1),a.dim(2)*b.dim(2)));
}

//=============================================================================================================================
template <typename T>
inline
typename RotationMatrix<T>::RotateMatrix RotationMatrix<T>::GetInverse(const RotateMatrix& mtx)
{
    T determinate_ = 0;
    
    auto MinorDeterminant = [&](uint32 row, uint32 col) {
        std::vector<std::vector<T>> elements(2);
        uint32 index = 0;
        for(auto j = 0; j < 3; ++j) {
            if(j == row) { continue; }
            for(auto i = 0; i < 3; ++i) {
                if(i == col) { continue; }
                elements.at(index).push_back(mtx.dim(j).dim(i));
            }
            if(j != row) { ++index; }
        }
        return elements.at(0).at(0)*elements.at(1).at(1) - elements.at(0).at(1)*elements.at(1).at(0);
    };
    
    for(auto i = 0; i < 3; ++i) {
        determinate_ += cofactor(i)(0)*mtx.dim(i).dim(0)*MinorDeterminant(i,0);
    }
    
    auto cp = [] ( T a, T b, T c, T d) { return a*d - b*c; };
    double f = 1./determinate_;
    ThreeVector<ThreeVector<T>> results(ThreeVector<T>(f*cp(mtx.dim(1).dim(1),mtx.dim(1).dim(2),mtx.dim(2).dim(1),mtx.dim(2).dim(2)),
                                                       f*cp(mtx.dim(0).dim(2),mtx.dim(0).dim(1),mtx.dim(2).dim(2),mtx.dim(2).dim(1)),
                                                       f*cp(mtx.dim(0).dim(1),mtx.dim(0).dim(2),mtx.dim(1).dim(1),mtx.dim(1).dim(2))),
                                        
                                        ThreeVector<T>(f*cp(mtx.dim(1).dim(2),mtx.dim(1).dim(0),mtx.dim(2).dim(2),mtx.dim(2).dim(0)),
                                                       f*cp(mtx.dim(0).dim(0),mtx.dim(0).dim(2),mtx.dim(2).dim(0),mtx.dim(2).dim(2)),
                                                       f*cp(mtx.dim(0).dim(2),mtx.dim(0).dim(0),mtx.dim(1).dim(2),mtx.dim(1).dim(0))),
                                        
                                        ThreeVector<T>(f*cp(mtx.dim(1).dim(0),mtx.dim(1).dim(1),mtx.dim(2).dim(0),mtx.dim(2).dim(1)),
                                                       f*cp(mtx.dim(0).dim(1),mtx.dim(0).dim(0),mtx.dim(2).dim(1),mtx.dim(2).dim(0)),
                                                       f*cp(mtx.dim(0).dim(0),mtx.dim(0).dim(1),mtx.dim(1).dim(0),mtx.dim(1).dim(1))));
    
    return results;
}

//=============================================================================================================================
template <typename T>
inline
typename RotationMatrix<T>::RotateMatrix RotationMatrix<T>::GetTranspose(const RotateMatrix& mtx)
{
    
    return RotateMatrix(ThreeVector<T>(mtx.cx().cx(),mtx.cy().cx(),mtx.cz().cx()),
                        ThreeVector<T>(mtx.cx().cy(),mtx.cy().cy(),mtx.cz().cy()),
                        ThreeVector<T>(mtx.cx().cz(),mtx.cy().cz(),mtx.cz().cz()));
}

//=============================================================================================================================

template <typename T>
inline ThreeVector<T> RotationMatrix<T>::GetEulerAngles(const RotateMatrix& rot)
{
    ThreeVector<double> angles(0.);
    
    // psi = x
    // theta = y
    // phi = z
    if(fabs(rot.dim(2).dim(0)) != 1.) {
        angles(1) = -asin(rot.dim(2).dim(0));
        angles(0) = atan2(rot.dim(2).dim(1)/cos(angles(1)),rot.dim(2).dim(2)/cos(angles(1)));
        angles(2) = atan2(rot.dim(1).dim(0)/cos(angles(1)),rot.dim(0).dim(0)/cos(angles(1)));
    } else {
        angles(2) = 0.;
        if(rot.dim(2).dim(0) == -1) {
            angles(1) = 0.5*M_PI;
            angles(0) = angles(2) + atan2(rot.dim(0).dim(1),rot.dim(0).dim(2));
        } else {
            angles(1) = -0.5*M_PI;
            angles(0) = -angles(2) + atan2(-rot.dim(0).dim(1),-rot.dim(0).dim(2));
        }
    }
    
    return angles;
}

//=============================================================================================================================





//=============================================================================================================================
//=============================================================================================================================
//=============================================================================================================================
//
//                  Rarely used constructors
//
//=============================================================================================================================
//=============================================================================================================================
//=============================================================================================================================
template <typename T>
inline
RotationMatrix<T>::RotationMatrix(std::vector<std::vector<T>> v)
:   current_matrix_(identity_matrix<T>()),current_angle_(ThreeVector<T>(static_cast<T>(0))),
    current_position_(ThreeVector<T>(static_cast<T>(0))),verbose(0),
    cofactor(ThreeVector<ThreeVector<T>>(ThreeVector<T>(1.,-1.,1),
                                         ThreeVector<T>(-1.,1.,-1),
                                         ThreeVector<T>(1.,-1.,1)))
{
    // first ThreeVector level of matrix is columns, each ThreeVector element in second level is row
    for(auto j = 0; j < 3; ++j) {
        current_matrix_.set(j,ThreeVector<T>(v.at(j)));
    }
    rotation_history_.resize(1,identity_matrix<T>());
    angle_history_.resize(1);
    position_history_.resize(1);
}

//=============================================================================================================================
template <typename T>
inline
RotationMatrix<T>::RotationMatrix(RotateMatrix v, ThreeVector<T> a)
:   current_matrix_(identity_matrix<T>()),current_angle_(ThreeVector<T>(static_cast<T>(0))),
    current_position_(ThreeVector<T>(static_cast<T>(0))),verbose(0),
    cofactor(ThreeVector<ThreeVector<T>>(ThreeVector<T>(1.,-1.,1),
                                         ThreeVector<T>(-1.,1.,-1),
                                         ThreeVector<T>(1.,-1.,1)))
{
    push_rotation_history(v);
    push_angle_history(a);
    position_history_.resize(1);
}

//=============================================================================================================================







//=============================================================================================================================
//=============================================================================================================================
//=============================================================================================================================
//
//                      Tedious Functions
//
//=============================================================================================================================
//=============================================================================================================================
//=============================================================================================================================
template <typename T>
inline
void RotationMatrix<T>::print_rotation(const RotateMatrix& rotation, std::ostream& out)
{
    out << std::endl;
    for(auto i = 0; i < 3; ++i) {
        out << "\t\t" << std::setw(50) << rotation.dim(i) << std::endl;
    }
    out << std::endl;
}

//=============================================================================================================================
template <typename T>
inline
void RotationMatrix<T>::pop_rotation_history()
{
    current_matrix_ = rotation_history_.back();
    rotation_history_.pop_back();
    
    if(rotation_history_.size() == 0) { Exception("Rotation Matrix","Rotation History has been popped too many times",Fatal); }
}

//=============================================================================================================================
template <typename T>
inline
void RotationMatrix<T>::push_rotation_history(const RotateMatrix& val)
{
    rotation_history_.push_back(current_matrix_);
    current_matrix_ = val;
    if(rotation_history_.size() > 500) { rotation_history_.pop_front(); }
}

//=============================================================================================================================
template <typename T>
inline
void RotationMatrix<T>::pop_angle_history()
{
    current_angle_ = angle_history_.back();
    angle_history_.pop_back();
    if(angle_history_.size() == 0) { Exception("Rotation Matrix","Angle History has been popped too many times",Fatal); }

}

//=============================================================================================================================
template <typename T>
inline
void RotationMatrix<T>::push_angle_history(const ThreeVector<T>& val)
{
    angle_history_.push_back(current_angle_);
    current_angle_ = val;
    if(angle_history_.size() > 500) { angle_history_.pop_front(); }

}

//=============================================================================================================================
template <typename T>
inline
void RotationMatrix<T>::pop_position_history()
{
    current_position_ = position_history_.back();
    position_history_.pop_back();
    if(position_history_.size() == 0) { Exception("Rotation Matrix","Position History has been popped too many times",Fatal); }
}

//=============================================================================================================================
template <typename T>
inline
void RotationMatrix<T>::push_position_history(const ThreeVector<T>& val)
{
    position_history_.push_back(current_position_);
    current_position_ = val;
    if(position_history_.size() > 500) { position_history_.pop_front(); }

}

//=============================================================================================================================
template <typename T>
inline
void RotationMatrix<T>::print_rotation_history(std::ostream& out)
{
    out << "\nROTATION HISTORY:" << std::endl;
    for(auto ite = rotation_history_.begin(); ite != rotation_history_.end(); ++ite) {
        out << "\t" << std::distance(rotation_history_.begin(),ite) << ":" << std::endl;
        for(auto j = 0; j < 3; ++j) {
            out << "\t\t" << (*ite)(j) << std::endl;
        }
    }
    out << "CURRENT ROTATION:" << std::endl;
    for(auto j = 0; j < 3; ++j) {
        out << "\t" << current_matrix_.dim(j) << std::endl;
    }
}
//=============================================================================================================================


#endif



