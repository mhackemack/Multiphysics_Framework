#ifndef prke_h_
#define prke_h_

#include "../Global/GlobalHeaders.h"
#include "PRKE_Tools.h"

class PRKE
{
  public:
    
    // Constructor
    PRKE (unsigned int);
    
    // Destructor
    ~PRKE () {}
    
	  void perform_prke_step();
	  
	  void set_prke_parameters
	  (
	    const vector<double>& time_reactivity_,
	    const vector<double>& reactivity_,
	    const vector<double>& beta_,
	    const vector<double>& lambda_,
	    const double& meanGenTime_,
	    const double& timeStep_,
	    const vector<double>& initial_solution_
	  );
	  
    // Accessors
    double get_sumBeta() {return sumBeta;}
    double get_meanGenTime() {return meanGenTime;}
    double get_timeStep() {return timeStep;}
    vector<double> get_lambda() {return lambda;}
    vector<double> get_beta() {return beta;}
    vector<double> get_time_reactivity() {return time_reactivity;}
    vector<double> get_reactivity() {return reactivity;}
    vector<double> get_solution() {return solution;}
    vector<double> get_old_solution() {return old_solution;}
    vector<double> get_initial_solution() {return initial_solution;}
    
    // Manipulators
    void set_beta(vector<double> beta_) {beta = beta_;}
    void set_sumBeta(double sumBeta_) {sumBeta = sumBeta_;}
    void set_lambda(vector<double> lambda_) {lambda = lambda_;}
    void set_meanGenTime(double meanGenTime_) {meanGenTime = meanGenTime_;}
    void set_timeStep(double timeStep_) {timeStep = timeStep_;}
    void set_smallTimeStep(double smallTimeStep_) {smallTimeStep = smallTimeStep_;}
    void set_time_reactivity(vector<double> time_reactivity_) {time_reactivity = time_reactivity_;}
    void set_reactivity(vector<double> reactivity_) {reactivity = reactivity_;}
    void set_solution(vector<double> solution_) {solution = solution_;}
    void set_old_solution(vector<double> old_solution_) {old_solution = old_solution_;}
    void set_initial_solution(vector<double> initial_solution_) {initial_solution = initial_solution_;}
	
  private:
    
    const unsigned int n_precursors;

	  vector<double> reactivity;
	  vector<double> time_reactivity;
	  vector<double> lambda;
	  vector<double> beta;

	  vector<double> solution;
	  vector<double> old_solution;
	  vector<double> initial_solution;
	  
	  double sumBeta;
	  double meanGenTime;
	  double timeStep;
	  double smallTimeStep;
	  
	  
	  inline double interpolate_reactivity(double time);
	  inline vector< vector<double>> create_A_matrix(unsigned int dim);
	  
	  inline void forward_euler_step();
	  inline void backward_euler_step();
	  inline void crank_nicholson_step();
	  
};


// Accessors


// Manipulators


//Interpolate Reactivity Linearly
double PRKE::interpolate_reactivity(double time) {

  unsigned int final = time_reactivity.size() - 1;

  if (reactivity.size() == 1) {
    return reactivity.front();
  }
  else {
    if (time < time_reactivity.front()) {
      return reactivity.front() - (reactivity.at(1) - reactivity.front())/(time_reactivity.at(1) - time_reactivity.front())*(time_reactivity.front() - time);
    }
    else if (time > time_reactivity.back()) {
      return reactivity.back() - (reactivity.back() - reactivity.at(final - 1))/(time_reactivity.back() - time_reactivity.at(final - 1))*(time - time_reactivity.at(final - 1));
    }
    else {
      for (unsigned int i=1; i<reactivity.size(); ++i) {
        if (time <= time_reactivity.at(i)) {
          return reactivity.at(i-1) + (reactivity.at(i) - reactivity.at(i-1))/(time_reactivity.at(i) - time_reactivity.at(i-1))*(time - time_reactivity.at(i-1));
        }
      }
    }
  }
  
};

//Create PRKE matrix
vector< vector<double>> PRKE::create_A_matrix(unsigned int dim) {
  
  vector< vector<double>> AMat(dim, vector<double>(dim,0.0));
  
  AMat.at(0).at(0) = (interpolate_reactivity(smallTimeStep) - sumBeta) / (meanGenTime);
  for (unsigned int i=1; i < n_precursors; ++i) {
    AMat.at(i).at(0) = beta.at(i-1) / meanGenTime;
    AMat.at(i).at(i) = -lambda.at(i-1);
    AMat.at(0).at(i) = lambda.at(i-1);
  }
  
  return AMat;
};


//Forward Euler Step
void PRKE::forward_euler_step() {
  
  vector< vector<double>> A = create_A_matrix(n_precursors + 1);
  vector< vector<double>> eye = PRKE_Tools::create_identity_matrix(n_precursors + 1);
  
  vector< vector<double>> tempMat = A;
  PRKE_Tools::multiply_by_scalar(&tempMat,smallTimeStep);
  
  tempMat = PRKE_Tools::matrix_addition(tempMat, eye);  
  
  solution = PRKE_Tools::matrix_vector_multiplication(tempMat, old_solution);
  
  return;
};




#endif
