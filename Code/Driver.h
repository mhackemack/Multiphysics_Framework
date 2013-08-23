//===========================================================================//
/*!
 * \class Driver
 * \brief
 *
 * Routine to print out heading and various outputs.
 */
//===========================================================================//

#ifndef DRIVER_H
#define DRIVER_H

// Global Headers
#include "Global/GlobalHeaders.h"
#include "Global/Deal2Headers.h"
#include "Global/Constants.h"

// Nuetronics Headers
#include "Diffusion/Neutronics.h"
#include "PRKE/PRKE.h"

// Thermal-Hydraulics Headers
#include "TH/HeatConduction.h"

// Input/Ouput Headers
#include "InputOutput/ParameterHandlerFunctions.h"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/regex.hpp>

template <int dim>
class Driver
{
  //Parameter Class Declaration
  // --------------------------
  public:
  class Parameters
  {
    public:
      Parameters ();

      void get_parameters (const std::string &proj_name);

      unsigned int phys_type;
      unsigned int phys_couple;
      unsigned int kin_couple;
      unsigned int temp_dep;
      unsigned int fe_degree;
      unsigned int number_synch_points;
      
      double tolerance;
      
      std::string folder_name;
      std::string project_name;
      std::string problem_name;
  };
  
  // Constructors and Functions
  // --------------------------
  public:
  // Constructor
  Driver (Parameters &parameters);

  // Destructor
  ~Driver ();
  
  // Local Parameters Class
  // Really lives in Main
  // --------------------
  Parameters  &parameters;
  
  // Routines
  // --------
  public:
  void execute_problem ();
  
  private:
  void initialize_problem ();
  void steady_state_solve ();
  void run_time_step ();
  void run_full_synchronous ();
  void run_two_step ();
  void run_JFNK_step ();
  void terminate_problem ();
  void neutronics_time_step ();
  void TH_time_step ();
  void dtstep ();
  void dtkin ();
  void dtTH ();
  
  bool check_transient_status ();
  bool check_next_synch_point ();
  void update_next_synch_point ();

  private:
  // Physics Classes
  // ---------------
  Neutronics<dim> neutronics;
  HeatConduction<dim> TH;
  
  // Variables
  // ---------
  public:
  

  private:
  double        tolerance;
  unsigned int  phys_type;
  unsigned int  phys_couple;
  unsigned int  kin_couple;
  unsigned int  temp_dep;
  unsigned int  fe_degree;
  
  // Transient Variables
  vector<double>            synch_times;
  vector<double>            dt_master_inp;
  vector<double>            dt_TH_inp;
  vector<double>            dt_master_kin_inp;
  vector<vector<double>>    dt_kin_inp;
  
  // Time-Step Controls
  vector<unsigned int>      master_time_step_selection;
  vector<unsigned int>      kin_time_step_selection;
  vector<unsigned int>      TH_time_step_selection;
  vector<unsigned int>      kin_time_step_estimator;
  vector<unsigned int>      TH_time_step_estimator;
  vector<double>            kin_err_tol;
  vector<double>            TH_err_tol;
  vector<double>            kin_min_time_step;
  vector<double>            TH_min_time_step;
  vector<double>            kin_max_time_step;
  vector<double>            TH_max_time_step;
  
  double                    final_time;
  double                    master_dt;
  double                    master_kin_dt;
  double                    master_TH_dt;
  
  unsigned int              synch_counter;
  unsigned int              TH_advance;
  vector<unsigned int>      kin_advance;
  
  vector<unsigned int>      phys_coupling;
  vector<unsigned int>      kin_coupling;
  vector<unsigned int>      master_step_selection;
  vector<unsigned int>      kin_step_selection;
  vector<unsigned int>      TH_step_selection;
  
  vector<bool>              change_kin_mesh;
  vector<bool>              change_TH_mesh;
  
  private:
  Triangulation<dim>    mesh;
  Triangulation<dim>    upload_new_mesh (std::string str_in);
};

//======================================================================================================================
// Driver Constructor
template <int dim>
Driver<dim>::Driver (Parameters &parameters) 
    :
    parameters (parameters),
    neutronics (parameters.project_name,parameters.fe_degree),
    TH (parameters.project_name,parameters.fe_degree)
{
    // Local Temporary Variables
    ParameterHandler kin_param;
    ParameterHandler TH_param;
    ParameterHandler trans_param;
    unsigned int     n_erg;
    
    synch_counter = 0;
    initialize_problem ();
}

//======================================================================================================================
// Driver Destructor
template <int dim>
Driver<dim>::~Driver ()
{
    //delete neutronics;
    //delete TH;
}

//======================================================================================================================
// Parameters Constructor
template <int dim>
Driver<dim>::Parameters::Parameters ()  {}

//======================================================================================================================
// Upload New Mesh into GridIn Class
template <int dim>
Triangulation<dim> Driver<dim>::upload_new_mesh (std::string str_in)
{
  
  // Local Variables
  GridIn<dim>   gridin;
  Triangulation<dim> mesh_;
  std::ifstream filename (parameters.folder_name + str_in);
  
  gridin.attach_triangulation (mesh_);
  gridin.read_msh (filename);
  
  return mesh_;
}

//======================================================================================================================
// Get Parameters Routine
template <int dim>
void Driver<dim>::Parameters::get_parameters (const std::string &proj_name)
{
    // Strings
    project_name = proj_name;
    folder_name = Constants::input_str + proj_name + "/";
    
    ParameterHandler proj_handler;
    ParameterHandler kin_handler;
    ParameterHandler TH_handler;
    ParameterHandler trans_handler;
    
    declare_project_parameters (proj_handler, folder_name);
    //declare_kinetics_parameters (kin_handler, folder_name);
    //declare_TH_parameters (TH_handler, folder_name);
    //declare_transient_parameters (trans_handler, folder_name);
    
    // Assign Project Parameters
    problem_name = proj_handler.get ("Problem Name");
    phys_type = proj_handler.get_integer ("Physics Type");
    phys_couple = proj_handler.get_integer ("Steady-State Physics Coupling");
    fe_degree = proj_handler.get_integer ("FE degree");
    tolerance = proj_handler.get_double ("tolerance");
    temp_dep = proj_handler.get_bool ("Temperature Dependence");
    number_synch_points = proj_handler.get_integer ("Number Synch Points");
    
    // Assign Physics Type
    switch (phys_type) 
    {
      case 0:
        phys_type = Constants::Physics_Neutronics_Only;
        temp_dep = false;
        break;
      case 1:
        phys_type = Constants::Physics_TH_Only;
        break;
      case 2:
        phys_type = Constants::Physics_Both;
        break;
    }
    
    // Assign Steady-State Physics Coupling
    // Only Matters if it is JFNK (NOT IMPLEMENTED YET)
    // Otherwise, gets updated every synch point
    // ------------------------------------------------
    switch (phys_couple)
    {
      case 0:
        phys_couple = Constants::synchronous;
        break;
      case 1:
        phys_couple = Constants::full_synchronous;
        break;
      case 2:
        phys_couple = Constants::asynchronous;
        break;
      case 3:
        phys_couple = Constants::full_asynchronous;
        break;
      case 4:
        phys_couple = Constants::TwoStep;
        break;
      case 5:
        phys_couple = Constants::JFNK;
        break;
    }
}

//======================================================================================================================
// Initialize All Classes Within Driver
template <int dim>
void Driver<dim>::initialize_problem ()
{
    // Local Variables
    ParameterHandler kin_handler;
    ParameterHandler TH_handler;
    ParameterHandler trans_handler;
    
    // Set Output Directories
    // ----------------------
    std::string out_folder1 = Constants::output_str + parameters.project_name + "/";
    std::string out_folder2 = Constants::output_str + parameters.project_name + "/" + parameters.problem_name + "/";
    std::string mk_comm1 = "exec mkdir " + out_folder1;
    std::string mk_comm2 = "exec mkdir " + out_folder2;
    std::string rm_comm1 = "exec rm -r " + out_folder1 + "*";
    std::string rm_comm2 = "exec rm -r " + out_folder2 + "*";
    struct stat sb1;
    struct stat sb2;
    
    char * sf1 = new char[out_folder1.length() + 1];
    char * sf2 = new char[out_folder2.length() + 1];
    char * sm1 = new char[mk_comm1.length() + 1];
    char * sm2 = new char[mk_comm2.length() + 1];
    char * sr1 = new char[rm_comm1.length() + 1];
    char * sr2 = new char[rm_comm2.length() + 1];
    
    std::strcpy (sf1,out_folder1.c_str());
    std::strcpy (sf2,out_folder2.c_str());
    std::strcpy (sm1,mk_comm1.c_str());
    std::strcpy (sm2,mk_comm2.c_str());
    std::strcpy (sr1,rm_comm1.c_str());
    std::strcpy (sr2,rm_comm2.c_str());
    
    if (stat(sf1, &sb1) == 0 && S_ISDIR(sb1.st_mode)) // Check if Project Name Folder exists
      if (stat(sf2, &sb2) == 0 && S_ISDIR(sb2.st_mode)) // Check if Problem Name Folder exists
      {
        cout << std::endl;
        system (sr2);
        cout << std::endl;
      }
      else
        system (sm2);
    else
    {
      system (sm1);
      system (sm2);
    }
    
    delete sf1;
    delete sf2;
    delete sm1;
    delete sm2;
    delete sr1;
    delete sr2;
    // End Setting up Output Directories
    // That is a mess...
    // ---------------------------------
    
    
    // Initialize Physics Domain
    switch(parameters.phys_type)
    {
        // Initialize Neutronics Problem Only
        case(Constants::Physics_Neutronics_Only):
        {
          Triangulation<dim> mesh_ (upload_new_mesh ("kin_mesh_0.msh"));
          neutronics.initialize_problem (mesh_);
          neutronics.set_prob_type (Constants::steady);
          declare_kinetics_parameters (kin_handler, parameters.folder_name);
          break;
        }
        
        // Initialize TH Problem Only
        case(Constants::Physics_TH_Only):
        {
          Triangulation<dim> mesh_ (upload_new_mesh ("TH_mesh_0.msh"));
          TH.initialize_problem (mesh_);
          TH.set_prob_type (Constants::steady);
          declare_TH_parameters (TH_handler, parameters.folder_name);
          break;
        }
        
        // Initialize Both Neutronics and TH Problems
        case(Constants::Physics_Both):
        {
          Triangulation<dim> mesh_ (upload_new_mesh ("kin_mesh_0.msh"));
          neutronics.initialize_problem (mesh_);
          neutronics.set_prob_type (Constants::steady);
          neutronics.set_TH_pointer (&TH);
          declare_kinetics_parameters (kin_handler, parameters.folder_name);
        }
        {
          Triangulation<dim> mesh_ (upload_new_mesh ("TH_mesh_0.msh"));
          TH.initialize_problem (mesh_);
          TH.set_prob_type (Constants::steady);
          TH.set_heat_src_pointers (&neutronics.heat_src, &neutronics.dof_handler);
          declare_TH_parameters (TH_handler, parameters.folder_name);
          break;
        }
    }
    
    // Collect All Transient Information
    declare_transient_parameters (trans_handler, parameters.folder_name);
    
    // Allocate Shared Physics Vectors
    synch_times.resize (parameters.number_synch_points, 0.0);
    phys_coupling.resize (parameters.number_synch_points, 0); // Only matters if both physics are turned on
    dt_master_inp.resize (parameters.number_synch_points, 0.0);
    master_step_selection.resize (parameters.number_synch_points, 0);
    master_time_step_selection.resize (parameters.number_synch_points, 0);
    
    // Allocate Physics Specific Vectors
    if (parameters.phys_type == Constants::Physics_Neutronics_Only || parameters.phys_type == Constants::Physics_Both)
    {
      dt_master_kin_inp.resize (parameters.number_synch_points, 0.0);
      dt_kin_inp.resize (parameters.number_synch_points, vector<double> (kin_handler.get_integer ("n_erg"), 0.0));
      kin_coupling.resize (parameters.number_synch_points, 0);
      kin_step_selection.resize (parameters.number_synch_points, 0);
      kin_time_step_selection.resize (parameters.number_synch_points, 0);
      kin_time_step_estimator.resize (parameters.number_synch_points, 0);
      change_kin_mesh.resize (parameters.number_synch_points, false);
      kin_advance.resize (neutronics.get_n_groups ());
      
      kin_err_tol.resize (parameters.number_synch_points, 0.0);
      kin_min_time_step.resize (parameters.number_synch_points, 0.0);
      kin_max_time_step.resize (parameters.number_synch_points, 0.0);
    }
    if (parameters.phys_type == Constants::Physics_TH_Only || parameters.phys_type == Constants::Physics_Both)
    {
      dt_TH_inp.resize (parameters.number_synch_points, 0.0);
      TH_step_selection.resize (parameters.number_synch_points, 0);
      change_TH_mesh.resize (parameters.number_synch_points, false);
      TH_time_step_selection.resize (parameters.number_synch_points, 0);
      TH_time_step_estimator.resize (parameters.number_synch_points, 0);
      
      TH_err_tol.resize (parameters.number_synch_points, 0.0);
      TH_min_time_step.resize (parameters.number_synch_points, 0.0);
      TH_max_time_step.resize (parameters.number_synch_points, 0.0);
    }
    
    // Final Time
    final_time = trans_handler.get_double ("Final Time");
    
    // Initialize All Transient Values
    for (unsigned int i=0; i<parameters.number_synch_points; ++i)
    {
        const unsigned int j = i + 1;
        const std::string section_name = "Synch Point " + Utilities::int_to_string(j,Utilities::needed_digits(j));
        trans_handler.enter_subsection(section_name);
        
        synch_times.at(i) = trans_handler.get_double ("Time");
        dt_master_inp.at(i) = trans_handler.get_double ("Master Time Step");
        phys_coupling.at(i) = trans_handler.get_integer ("Physics Coupling");
        master_time_step_selection.at(i) = trans_handler.get_integer ("Master Time Step Selection");
        
        switch(parameters.phys_type)
        {
            // Initialize Neutronics Problem Only
            case(Constants::Physics_Neutronics_Only):
              dt_master_kin_inp.at(i) = trans_handler.get_double("Master Kinetics Time Step");
              change_kin_mesh.at(i) = trans_handler.get_bool("Change Kinetics Mesh");
              kin_coupling.at(i) = trans_handler.get_integer("Kinetics Coupling");
              kin_time_step_selection.at(i) = trans_handler.get_integer ("Kinetics Time Step Selection");
              kin_time_step_estimator.at(i) = trans_handler.get_integer ("Kinetics Time Step Estimator");
              for (unsigned int group=0;group<neutronics.get_n_groups ();++group)
                dt_kin_inp.at(i).at(group) = trans_handler.get_double("Kinetics Time Step Group "+Utilities::int_to_string(group+1,Utilities::needed_digits(group+1)));
              
              kin_err_tol.at(i) = trans_handler.get_double ("Kinetics Adaptive Error Tolerance");
              kin_min_time_step.at(i) = trans_handler.get_double ("Kinetics Min Time Step");
              kin_max_time_step.at(i) = trans_handler.get_double ("Kinetics Max Time Step");
              break;
            
            // Initialize TH Problem Only
            case(Constants::Physics_TH_Only):
              TH_time_step_selection.at(i) = trans_handler.get_integer ("TH Time Step Selection");
              TH_time_step_estimator.at(i) = trans_handler.get_integer ("KinetTHics Time Step Estimator");
              dt_TH_inp.at(i) = trans_handler.get_double("Master TH Time Step");
              change_TH_mesh.at(i) = trans_handler.get_bool("Change TH Mesh");
              
              TH_err_tol.at(i) = trans_handler.get_double ("TH Adaptive Error Tolerance");
              TH_min_time_step.at(i) = trans_handler.get_double ("TH Min Time Step");
              TH_max_time_step.at(i) = trans_handler.get_double ("TH Max Time Step");
              break;
            
            // Initialize Both Neutronics and TH Problems
            case(Constants::Physics_Both):
              dt_master_kin_inp.at(i) = trans_handler.get_double("Master Kinetics Time Step");
              kin_coupling.at(i) = trans_handler.get_integer("Kinetics Coupling");
              kin_time_step_selection.at(i) = trans_handler.get_integer ("Kinetics Time Step Selection");
              kin_time_step_estimator.at(i) = trans_handler.get_integer ("Kinetics Time Step Estimator");
              change_TH_mesh.at(i) = trans_handler.get_bool("Change Kinetics Mesh");
              for (unsigned int group=0;group<neutronics.get_n_groups ();++group)
                dt_kin_inp.at(i).at(group) = trans_handler.get_double("Kinetics Time Step Group "+Utilities::int_to_string(group+1,Utilities::needed_digits(group+1)));
              
              kin_err_tol.at(i) = trans_handler.get_double ("Kinetics Adaptive Error Tolerance");
              kin_min_time_step.at(i) = trans_handler.get_double ("Kinetics Min Time Step");
              kin_max_time_step.at(i) = trans_handler.get_double ("Kinetics Max Time Step");
              
              TH_time_step_selection.at(i) = trans_handler.get_integer ("TH Time Step Selection");
              TH_time_step_estimator.at(i) = trans_handler.get_integer ("KinetTHics Time Step Estimator");
              dt_TH_inp.at(i) = trans_handler.get_double("Master TH Time Step");
              change_TH_mesh.at(i) = trans_handler.get_bool("Change TH Mesh");
              
              TH_err_tol.at(i) = trans_handler.get_double ("TH Adaptive Error Tolerance");
              TH_min_time_step.at(i) = trans_handler.get_double ("TH Min Time Step");
              TH_max_time_step.at(i) = trans_handler.get_double ("TH Max Time Step");
              break;
        }
        
        trans_handler.leave_subsection ();
    }
    
    return;
}

//======================================================================================================================
// Execute Problem
template <int dim>
void Driver<dim>::execute_problem ()
{
  
  // Local Variables
  // ---------------
  bool transient_running (true);
  
  // Run Steady-State Solve if called for
  // Currently, this is always required
  // ------------------------------------
  steady_state_solve ();
  
  // Terminate Problem if no transient input is specified
  if (parameters.number_synch_points == 0)
  {
    terminate_problem ();
    return;
  }
  else if (fabs(synch_times[0]) < Constants::SMALL)
    update_next_synch_point ();
    
  // Perform Transient Steps
  // Exit when all physics reach final time
  // --------------------------------------
  do
  {
    dtstep ();
    run_time_step ();
    
    transient_running = check_transient_status ();
    if (check_next_synch_point ()) update_next_synch_point ();
  }
  while(transient_running);
  
  
  terminate_problem ();
  return;
}

//======================================================================================================================
// Steady-State Calculation Subroutine
template <int dim>
void Driver<dim>::steady_state_solve ()
{
    std::cout << std::setprecision (12) << std::fixed;
    
    if (parameters.phys_type == Constants::Physics_Neutronics_Only)
    {
      neutronics.eigenvalue_calc ();
      neutronics.set_initial_neutron_population (neutronics.get_total_neutron_population(0,0));
      neutronics.set_prob_type(Constants::transient);
      return;
    }
    else if (parameters.phys_type == Constants::Physics_TH_Only)
    {
      TH.steady_state_solve ();
      TH.set_initial_temp_integral (TH.get_total_energy_save ());
      TH.set_prob_type(Constants::transient);
      return;
    }
    
    double err_kin, err_TH;
    double kin_int, kin_int0;
    double TH_int, TH_int0;
    unsigned int iteration = 1;
    
    err_kin = 0.0; err_TH = 0.0;
    kin_int = 0.0; kin_int0 = 0.0;
    TH_int = 0.0; TH_int0 = 0.0;
    
    // Set Initial Heat Source to Allow for Program Execution
    neutronics.update_heat_src (0.0, 0);
    
    do
    {
        // TH Solves
        TH.steady_state_solve ();
        TH_int = fabs(TH.get_total_energy_save ());
        err_TH = fabs(TH_int - TH_int0) / fabs(TH_int);
        TH_int0 = TH_int;
        
        // Neutronics Solves
        neutronics.eigenvalue_calc ();
        neutronics.update_heat_src (0.0, 0);
        kin_int = fabs(neutronics.get_total_fission_source ());
        err_kin = fabs(kin_int - kin_int0)/fabs(kin_int);
        kin_int0 = kin_int;
        
        ++iteration;
    }
    while((err_kin > parameters.tolerance) && (err_TH > parameters.tolerance) && (iteration < Constants::MAX_ITERATION));
    
    // Clean up each physics
    neutronics.set_initial_neutron_population (neutronics.get_total_neutron_population(0,0));
    TH.set_initial_temp_integral (TH.get_total_energy_save ());
    neutronics.set_prob_type(Constants::transient);
    TH.set_prob_type(Constants::transient);
    
    return;
}

//======================================================================================================================
// Run a time step
template <int dim>
void Driver<dim>::run_time_step ()
{
  // Local Variables
  bool  any_advance (false);
  bool  kin_advance_local (false);
  bool  TH_advance_local  (false);
  
  // Check if any advance flags are set
  if (parameters.phys_type == Constants::Physics_TH_Only || parameters.phys_type == Constants::Physics_Both)
    if (TH_advance == 1)
    {
      any_advance = true;
      TH_advance_local = true;
    }
  if (parameters.phys_type == Constants::Physics_Neutronics_Only || parameters.phys_type == Constants::Physics_Both)
    for (unsigned int g=0; g<neutronics.get_n_groups (); ++g)
      if (kin_advance.at(g) == 1)
      {
        any_advance = true;
        kin_advance_local = true;
      }
  
  // Return if no advance flags are set - not sure here...
  if (!any_advance) return;
  
  // Go through logic to perform time step
  if (parameters.phys_couple == Constants::JFNK) // Not implemented yet...
  {
    run_JFNK_step ();
    return;
  }
  else if (parameters.phys_type == Constants::Physics_Both)
  {
    if (parameters.phys_couple == Constants::TwoStep)
    {
      run_two_step ();
      return;
    }
  }
  else if (parameters.phys_type == Constants::Physics_Neutronics_Only)
  {
    neutronics.set_all_advance (kin_advance);
    neutronics_time_step ();
    neutronics.update_solution ();
    neutronics.set_all_advance (0);
    return;
  }
  else if (parameters.phys_type == Constants::Physics_TH_Only)
  {
    TH.set_advance (1);
    TH_time_step ();
    TH.update_solution ();
    TH.set_advance (0);
    return;
  }
  
  // Perform Asynchronous Step
  if (parameters.phys_couple == Constants::asynchronous || parameters.phys_couple == Constants::full_asynchronous)
  {
    if (TH_advance_local)
    {
      neutronics.set_all_advance (0);
      neutronics.update_heat_src (TH.get_next_time (), 1);
      TH_time_step ();
      
      TH.update_solution ();
      TH.set_advance (0);
    }
    else if (kin_advance_local)
    {
      TH.set_advance (0);
      neutronics_time_step ();
      
      neutronics.update_solution ();
      neutronics.set_all_advance (0);
      
    }
  }
  else if (parameters.phys_couple == Constants::synchronous)
  {
    neutronics.set_all_advance (0);
    neutronics.update_heat_src (TH.get_next_time (), 1);
    
    TH.set_advance (1);
    TH_time_step ();
    TH.set_advance (0);
    
    neutronics.set_all_advance (1);
    neutronics_time_step ();
    
    neutronics.update_solution ();
    TH.update_solution ();
    
    neutronics.set_all_advance (0);
    TH.set_advance (0);
    
  }
  else if (parameters.phys_couple == Constants::full_synchronous)
    run_full_synchronous ();
  
  
  return;
}

//======================================================================================================================
// Run the JFNK Step
template <int dim>
void Driver<dim>::run_JFNK_step ()
{

}

//======================================================================================================================
// Run the fully-synchronous problem
template <int dim>
void Driver<dim>::run_full_synchronous ()
{
  std::cout << std::setprecision (12) << std::fixed;
  
  // Local Variables
  double dt = master_dt;
  
  double err_kin, err_TH;
  double kin_int, kin_int0;
  double TH_int, TH_int0;
  unsigned int iteration = 1;
  
  err_kin = 0.0; err_TH = 0.0;
  kin_int = 0.0; kin_int0 = 0.0;
  TH_int = 0.0; TH_int0 = 0.0;
    
  // Set Initial Heat Source to Allow for Program Execution
  neutronics.set_all_advance (0);
  neutronics.update_heat_src (TH.get_next_time (), 1);
  TH.set_advance (1);
  
  do
  {
    TH_time_step ();
    TH_int = fabs(TH.get_total_energy_save ());
    err_TH = fabs(TH_int - TH_int0) / fabs(TH_int);
    TH_int0 = TH_int;
    TH.set_advance (1);
    
    neutronics.set_all_advance (1);
    neutronics_time_step ();
    kin_int = fabs(neutronics.get_total_fission_source ());
    err_kin = fabs(kin_int - kin_int0)/fabs(kin_int);
    kin_int0 = kin_int;
    neutronics.update_heat_src (TH.get_next_time (), 0);
    
    ++iteration;
  }
  while((err_kin > parameters.tolerance) && (err_TH > parameters.tolerance) && (iteration < Constants::MAX_ITERATION));
  
  // Clean up Physics Classes
  neutronics.set_all_advance (1);
  TH.set_advance (1);
  neutronics.update_solution ();
  TH.update_solution ();
    
  return;
}

//======================================================================================================================
// Run the two-step time step
template <int dim>
void Driver<dim>::run_two_step ()
{
  // Local Variables
  double dt = master_dt;
  //if () 
  
  neutronics.set_all_advance (0);
  neutronics.update_heat_src (TH.get_next_time (), 1);
  
  // Run first TH step
  TH.set_advance (1);
  TH_time_step ();
  
  // Run first neutronics step
  neutronics.set_all_advance (1);
  neutronics_time_step ();
  neutronics.update_heat_src (TH.get_next_time (), 0);
  
  // Run second TH step
  TH.set_advance (1);
  TH_time_step ();
  
  // Run second neutronics step
  neutronics.set_all_advance (1);
  neutronics_time_step ();
  
  // Clean up Physics Classes
  neutronics.update_solution ();
  TH.update_solution ();
    
  return;
}

//======================================================================================================================
// Terminate Problem Subroutine
template <int dim>
void Driver<dim>::terminate_problem ()
{
  if (parameters.phys_type == Constants::Physics_Neutronics_Only || parameters.phys_type == Constants::Physics_Both)
    neutronics.terminate_problem ();
  if (parameters.phys_type == Constants::Physics_TH_Only || parameters.phys_type == Constants::Physics_Both)
    TH.terminate_problem ();
    
  return;
}

//======================================================================================================================
// Neutronics Time-Step Calculation Subroutine
template <int dim>
void Driver<dim>::neutronics_time_step ()
{
  neutronics.run_time_step ();
  return;
}

//======================================================================================================================
// Thermal Hydraulics Time-Step Calculation Subroutine
template <int dim>
void Driver<dim>::TH_time_step ()
{
  TH.run_time_step ();
  return;
}

//======================================================================================================================
// Sort Physics Time Steps
template <int dim>
void Driver<dim>::dtstep ()
{
    // Local Variables
    double          TH_time (0.0);
    double          TH_next_time (0.0);
    double          TH_dt_next (0.0);
    double          TH_dtest (0.0);
    double          kin_dtest (0.0);
    
    vector<double>  kin_dt_next;
    vector<double>  kin_time;
    vector<double>  kin_next_time;
    
    double          min_time (0.0);
    double          kin_min_time (1.0e32);
    unsigned int    min_group_counter;
    
    if (parameters.phys_type == Constants::Physics_Neutronics_Only)
    {
      kin_dt_next.resize (neutronics.get_n_groups (), 0.0);
      kin_next_time.resize (neutronics.get_n_groups (), 0.0);
      neutronics.set_all_advance (0);
      TH_advance = 0;
      kin_time = neutronics.get_all_time ();
      if (parameters.kin_couple != Constants::asynchronous && parameters.kin_couple != Constants::full_asynchronous) // Synchronous
      {
        kin_next_time = neutronics.get_all_time ();
        if (kin_time_step_selection.at(synch_counter) != Constants::manual)
        {
          double dt_temp (1.0e12);
          kin_dt_next = neutronics.get_all_dtest ();
          for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
            if (fabs(kin_dt_next.at(g)) < dt_temp)
              dt_temp = kin_dt_next.at(g);
          
          if (fabs(dt_temp) < kin_min_time_step.at(synch_counter))
            dt_temp = kin_min_time_step.at(synch_counter);
          else if (fabs(dt_temp) > kin_max_time_step.at(synch_counter))
            dt_temp = kin_max_time_step.at(synch_counter);
          
          for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
            kin_dt_next.at(g) = dt_temp;
        }
        else
        {
          for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          {
            if (fabs(master_dt) > Constants::SMALL)
              kin_dt_next.at(g) = master_dt;
            else
              kin_dt_next.at(g) = master_kin_dt;
              
            kin_next_time.at(g) += kin_dt_next.at(g);
          }
        }
        
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          kin_advance.at(g) = 1;
         
      }
      else  // Asynchronous Coupling
      {
        kin_dt_next.resize (neutronics.get_n_groups (), 0.0);
        kin_next_time.resize (neutronics.get_n_groups (), 0.0);
        kin_next_time = neutronics.get_all_time ();
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          kin_advance.at(g) = 0;
        
        // Determine dtnext and next_time quantities for all energy groups
        if (kin_time_step_selection.at(synch_counter) != Constants::manual)
        {
          kin_dt_next = neutronics.get_all_dtest ();
          for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
            if (fabs(kin_dt_next.at(g)) < kin_min_time_step.at(synch_counter))
              kin_dt_next.at(g) = kin_min_time_step.at(synch_counter);
            else if (fabs(kin_dt_next.at(g)) > kin_max_time_step.at(synch_counter))
              kin_dt_next.at(g) = kin_max_time_step.at(synch_counter);
        }
        else
        {
          for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
            kin_dt_next.at(g) = dt_kin_inp.at(synch_counter).at(g);
        }
        
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          kin_next_time.at(g) += kin_dt_next.at(g);
        
        // Find Minimum Next Time between Energy Groups
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
        {
          if (fabs(kin_next_time.at(g)) > Constants::SMALL && fabs(kin_next_time.at(g)) < fabs(kin_min_time))
          {
            kin_min_time = fabs(kin_next_time.at(g));
            min_group_counter = g;
          }
        }
        
        // Determine if Multiple Energy Groups Have same Minimum Time
        // Set their advance flags and dtnext values also
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          if ( fabs(kin_min_time - kin_next_time.at(g)) < Constants::SMALL)
            if (parameters.kin_couple == Constants::asynchronous)
              kin_advance.at(g) = 1;
            else if (parameters.kin_couple == Constants::full_asynchronous)
            {
              kin_advance.at(g) = 1;
              break;
            }
        
      }
      
      if (kin_min_time > synch_times.at(synch_counter))
        kin_min_time = synch_times.at(synch_counter);
      else if (kin_min_time > final_time)
        kin_min_time = final_time;
      
      for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
      {
        if (fabs(kin_time.at(g) - kin_min_time) < Constants::SMALL)
          kin_advance.at(g) = 0;
        else if (kin_time.at(g) > kin_min_time)
          kin_advance.at(g) = 0;
        else if (kin_time.at(g) < kin_min_time && kin_next_time.at(g) >= kin_min_time)
        {
          if (fabs(kin_next_time.at(g) - kin_min_time) < Constants::SMALL)
          {
            kin_advance.at(g) = 1;
            kin_dt_next.at(g) = kin_min_time - kin_time.at(g);
          }
          else if (fabs(kin_next_time.at(g)) >= synch_times.at(synch_counter) || fabs(kin_next_time.at(g)) >= final_time)
          {
            kin_advance.at(g) = 1;
            kin_dt_next.at(g) = kin_min_time - kin_time.at(g);
          }
        }
      }
      
      neutronics.set_all_dtnext (kin_dt_next);
      
    }
    else if (parameters.phys_type == Constants::Physics_TH_Only)
    {
      TH_time = TH.get_time ();
      TH_advance = 1;
      if (TH_time_step_selection.at(synch_counter) != Constants::manual)
        TH_dt_next = TH.get_dtest ();
      else
        if (fabs(master_dt) > Constants::SMALL)
          TH_dt_next = master_dt;
        else
          TH_dt_next = master_TH_dt;
      
      TH_next_time = TH_time + TH_dt_next;
      
      if (TH_next_time > synch_times.at(synch_counter))
        TH_dt_next = synch_times.at(synch_counter) - TH_time;
      else if (TH_next_time > final_time)
        TH_dt_next = final_time - TH_time;
      
      TH.set_dtnext (TH_dt_next);
    }
    else
    {
      double time_temp (1.0e14);
      kin_dt_next.resize (neutronics.get_n_groups (), 0.0);
      kin_next_time.resize (neutronics.get_n_groups (), 0.0);
      TH_time = TH.get_time ();
      kin_time = neutronics.get_all_time ();
      TH_advance = 0;
      kin_advance.resize (neutronics.get_n_groups (), 0);
      
      if (parameters.phys_couple != Constants::asynchronous && parameters.phys_couple != Constants::full_asynchronous)
      {
        Assert (fabs(master_dt) < Constants::SMALL, ExcZero ())
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          kin_dt_next.at(g) = master_dt;
        TH_dt_next = master_dt;
        
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          kin_next_time.at(g) = kin_time.at(g) + kin_dt_next.at(g);
        TH_next_time = TH_time + TH_dt_next;
        
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          kin_advance.at(g) = 1;
        TH_advance = 1;
        
      }
      else
      {
        // TH Times
        if (TH_time_step_selection.at(synch_counter) != Constants::manual)
          TH_dt_next = TH.get_dtest ();
        else
          TH_dt_next = master_TH_dt;
        TH_next_time = TH_time + TH_dt_next;
        
        // Neutronics Times
        if (parameters.kin_couple != Constants::asynchronous && parameters.kin_couple != Constants::full_asynchronous)
          kin_dt_next.resize (neutronics.get_n_groups (), master_kin_dt);
        else if (kin_time_step_selection.at(synch_counter) != Constants::manual)
        {
          kin_dt_next = neutronics.get_all_dtest ();
          for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
            if (fabs(kin_dt_next.at(g)) < kin_min_time_step.at(synch_counter))
              kin_dt_next.at(g) = kin_min_time_step.at(synch_counter);
            else if (fabs(kin_dt_next.at(g)) > kin_max_time_step.at(synch_counter))
              kin_dt_next.at(g) = kin_max_time_step.at(synch_counter);
        }
        else
          for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
            kin_dt_next.at(g) = dt_kin_inp.at(synch_counter).at(g);
        
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          kin_next_time.at(g) = kin_time.at(g) + kin_dt_next.at(g);
        
      }
      
      for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
        if (kin_next_time.at(g) < time_temp)
          time_temp = kin_next_time.at(g);
      if (TH_next_time < time_temp)
        time_temp = TH_next_time;
        
      min_time = time_temp;
      if (min_time > synch_times.at(synch_counter))
        min_time = synch_times.at(synch_counter);
      else if (min_time > final_time)
        min_time = final_time;
        
      if (parameters.phys_couple != Constants::asynchronous && parameters.phys_couple != Constants::full_asynchronous)
      {
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          kin_dt_next.at(g) = min_time - kin_time.at(g);
        TH_dt_next = min_time - TH_time;
        
        neutronics.set_all_dtnext (kin_dt_next);
        TH.set_dtnext (TH_dt_next);
        neutronics.set_all_advance (1);
        TH.set_advance (1);
      }
      else
      {
        if (fabs(TH_next_time - min_time) < Constants::SMALL)
          TH_advance = 1;
        else if (TH_next_time > min_time && (TH_next_time > synch_times.at(synch_counter) || TH_next_time > final_time))
        {
          TH_advance = 1;
          TH_dt_next = min_time - TH_time;
        }
        else
          TH_advance = 0;
        
        for (unsigned int g = 0; g < neutronics.get_n_groups (); ++g)
          if (fabs(TH_next_time - min_time) < Constants::SMALL)
            kin_advance.at(g) = 1;
          else if (TH_next_time > min_time && (TH_next_time > synch_times.at(synch_counter) || TH_next_time > final_time))
          {
            kin_advance.at(g) = 1;
            kin_dt_next.at(g) = min_time - kin_time.at(g);
          }
          else
            kin_advance.at(g) = 0;
        
      }
      
      neutronics.set_all_dtnext (kin_dt_next);
      TH.set_dtnext (TH_dt_next);
      neutronics.set_all_advance (kin_advance);
      TH.set_advance (TH_advance);
      
    }
      
    return;
}

//======================================================================================================================
// Calculate Kinetics Time Step
template <int dim>
void Driver<dim>::dtkin()
{
    
    // Local Variables
    std::vector<double> tnext (neutronics.get_n_groups (), 0.0);
    
    // Get Next times
    for (unsigned int group = 0; group < neutronics.get_n_groups (); ++group)
      tnext[group] = neutronics.energy_groups[group]->get_next_time ();
    
    return;
}

//======================================================================================================================
// Calculate Thermal Hydraulics Time Step
template <int dim>
void Driver<dim>::dtTH()
{
    
  // Local Variables
  double tnext (TH.get_next_time());

  return;
}

//======================================================================================================================
// Check if all physics have reach final time
template <int dim>
bool Driver<dim>::check_transient_status ()
{
  
  // Local Variables
  bool status (true);
  
  // Check Times for All Physics
  if (parameters.phys_type == Constants::Physics_Neutronics_Only || parameters.phys_type == Constants::Physics_Both)
    if (fabs(neutronics.get_min_time () - final_time) > Constants::SMALL && fabs(neutronics.get_min_time ()) < final_time)
      status = false;
  if (parameters.phys_type == Constants::Physics_TH_Only || parameters.phys_type == Constants::Physics_Both)
    if (fabs(TH.get_time () - final_time) > Constants::SMALL && fabs(TH.get_time ()) < final_time)
      status = false;
  
  // Swap to return true if transient is not finished
  if (status == false)
    status = true;
  else 
    status = false;
    
  return status;
}

//======================================================================================================================
// Check if all physics have reached next synch point
template <int dim>
bool Driver<dim>::check_next_synch_point ()
{
  
  // Local Variables
  bool status (true);
  
  // Check Times for All Physics
  if (parameters.phys_type == Constants::Physics_Neutronics_Only || parameters.phys_type == Constants::Physics_Both)
    if (fabs(neutronics.get_min_time () - synch_times[synch_counter]) > Constants::SMALL)
      status = false;
  if (parameters.phys_type == Constants::Physics_TH_Only || parameters.phys_type == Constants::Physics_Both)
    if (fabs(TH.get_time () - synch_times[synch_counter]) > Constants::SMALL)
      status = false;
  
  
  return status;
}

//======================================================================================================================
// Update problem at next synch point
template <int dim>
void Driver<dim>::update_next_synch_point ()
{
  
  // Return if Current Synch Time Ends at Final Time
  if (fabs(final_time - synch_times.at (synch_counter)) < Constants::SMALL) return;
  
  // Print header
  print_synch_point_header (synch_counter+1, synch_times.at(synch_counter));
  
  // Update Meshes First if Required
  if (parameters.phys_type == Constants::Physics_Neutronics_Only || parameters.phys_type == Constants::Physics_Both)
    if (change_kin_mesh[synch_counter])
    {
      std::string filename = "kin_mesh_" + Utilities::int_to_string(synch_counter+1, Utilities::needed_digits(synch_counter+1)) + ".msh";
      Triangulation<dim> mesh_ (upload_new_mesh (filename));
      neutronics.update_mesh (mesh_);
    }
  if (parameters.phys_type == Constants::Physics_TH_Only || parameters.phys_type == Constants::Physics_Both)
    if (change_TH_mesh[synch_counter])
    {
      std::string filename = "TH_mesh_" + Utilities::int_to_string(synch_counter+1, Utilities::needed_digits(synch_counter+1)) + ".msh";
      Triangulation<dim> mesh_ (upload_new_mesh (filename));
      TH.update_mesh (mesh_);
    }
  
  // Update Counter
  if (parameters.number_synch_points != synch_counter + 1) 
    ++synch_counter;
  else
    return;
  
  // Change Parameter Variables
  parameters.phys_couple = phys_coupling.at(synch_counter);
  
  // Change Other Common Variables
  master_dt = dt_master_inp.at(synch_counter);
  
  // Update Timing Controls with New Synch Point Information
  if (parameters.phys_type == Constants::Physics_Neutronics_Only || parameters.phys_type == Constants::Physics_Both)
  {
    master_kin_dt = dt_master_kin_inp.at(synch_counter);
    if (parameters.phys_couple != Constants::asynchronous && parameters.phys_couple != Constants::full_asynchronous)
      parameters.kin_couple = parameters.phys_couple;
    else
      parameters.kin_couple = kin_coupling.at(synch_counter);
    
    neutronics.set_time_step_selection (kin_time_step_selection.at(synch_counter));
    neutronics.set_dt_estimator (kin_time_step_estimator.at(synch_counter));
    neutronics.set_adaptive_err_tol (kin_err_tol.at(synch_counter));
  }
  if (parameters.phys_type == Constants::Physics_TH_Only || parameters.phys_type == Constants::Physics_Both)
  {
    master_TH_dt = dt_TH_inp.at(synch_counter);
    TH.set_time_step_selection (TH_time_step_selection.at(synch_counter));
    TH.set_dt_estimator (TH_time_step_estimator.at(synch_counter));
    TH.set_adaptive_err_tol (TH_err_tol.at(synch_counter));
  }
  
  return;
}

//======================================================================================================================


#endif

