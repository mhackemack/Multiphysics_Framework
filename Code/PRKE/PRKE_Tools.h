#ifndef prke_tools_h_
#define prke_tools_h_

#include "../Global/GlobalHeaders.h"

namespace PRKE_Tools {

// Identity Matrix Routine
vector< vector<double>> create_identity_matrix(unsigned int dim) {
  
  vector< vector<double>> IMat (dim, vector<double> (dim, 0.0));
  
  for (unsigned int i=0; i < dim; ++i) {
    IMat.at(i).at(i) = 1.0;
  }
  
  return IMat;
};

// Matrix Multiplication Routine
template <typename T>
vector< vector<T>> matrix_multiplication(vector< vector<T>> AMat, vector< vector<T>> BMat) {
  
  vector< vector<T>> finalMat (AMat.size(), vector<T> (BMat.at(0).size(), 0.0));
  
  for (unsigned int i=0; i < AMat.size(); ++i) {
    for (unsigned int j=0; j < AMat.at(i).size(); ++j) {
      
    }
  }
  
  return finalMat;
};

// Matrix_Vector Multiplication Routine
template <typename T>
vector<T> matrix_vector_multiplication(vector< vector<T>> AMat, vector<T> BVect) {
   
  vector<T> finalVect (BVect.size(), 0.0);
  
  for (unsigned int i=0; i < AMat.size(); ++i) {
    for (unsigned int j=0; j < AMat.at(i).size(); ++j) {
      finalVect.at(i) += AMat.at(i).at(j) * BVect.at(j);
    }
  }
  
  return finalVect;
};

template <typename T>
void multiply_by_scalar(vector< vector<T>> *v, T x) {
  
  for (unsigned int i=0; i < v->size(); ++i) {
    for (unsigned int j=0; j < v->at(i).size(); ++j) {
      v->at(i).at(j) *= x;
    }
  }
  
  return;
};

template <typename T>
void multiply_by_scalar(vector<T> *v, T x) {
  
  for (unsigned int i=0; i < v->size(); ++i) {
    v->at(i) *= x;
  }
  
  return;
};

// Matrix Addition Routine
template <typename T>
vector< vector<T>> matrix_addition(vector< vector<T>> AMat, vector< vector<T>> BMat) {
  
  vector< vector<T>> finalMat (AMat.size(), vector<T> (AMat.at(0).size(), 0.0));
  
  for (unsigned int i=0; i < AMat.size(); ++i) {
    for (unsigned int j=0; j < AMat.at(i).size(); ++j) {
      finalMat.at(i).at(j) = AMat.at(i).at(j) + BMat.at(i).at(j);
    }
  }
  
  return finalMat;
};

}




#endif
