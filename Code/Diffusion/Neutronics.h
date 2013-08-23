#ifndef neutronics_h_
#define neutronics_h_

//Global Headers
#include "../Global/GlobalHeaders.h"
#include "../Global/Deal2Headers.h"
#include "../Global/Constants.h"

#include "../InputOutput/PrintStatements.h"
#include "../InputOutput/ParameterHandlerFunctions.h"

#include "MaterialData.h"
#include "EnergyGroup.h"
//#include "QS_Tools.h"

#include "../TH/HeatConduction.h"

//======================================================================================================================
  template <int dim>
  class Neutronics
  {
    public:
      class Parameters
      {
        public:
          Parameters ();

          void get_parameters (ParameterHandler &prm);
          
          unsigned int n_mat;
          unsigned int n_groups;
          unsigned int n_prec;
          unsigned int fe_degree;
          
          unsigned int coupling;
          unsigned int prob_type;
          unsigned int phys_type;
          unsigned int extrap;
          unsigned int time_step_selection;
          unsigned int dt_estimator;
          unsigned int PRKE_amplitude_count;
          
          std::string folder_name;
          std::string project_name;
          std::string problem_name;

          double tolerance;
          double initial_power;
          double adaptive_err_tol;
          
          bool  temp_dep;
          bool  use_adjoint_weighting;
          bool  print_output;
      };
    
    public:
      // Constructor
      Neutronics (const std::string inp_path, const unsigned int fe_deg);

      // Destructor
      ~Neutronics ();
      
      // Initialization Routines
      void set_parameters ();
      void initialize_problem(Triangulation<dim> &mesh_);
      void update_mesh (Triangulation<dim> &mesh_);

      // Calculation Subroutines
      void eigenvalue_calc ();
      void adjoint_calc ();
      void run_time_step ();
      void run_PCQS_time_step ();
      void run_IQS_time_step ();
      void update_precursors ();
      void output_precursors ();
      void update_solution ();
      void update_powers (const double t, const unsigned int flag);
      void compile_kinetics_physics (const double t);
      void terminate_problem ();
      
      // QS Tools
      void update_PCQS_amplitude ();

      // Manipulators
      inline void set_prob_type (const unsigned int p_type);
      inline void set_phys_type (unsigned int p_type) {parameters.phys_type=p_type;}
      inline void set_time_step_selection (unsigned int p);
      void set_dt_estimator (unsigned int p);
      inline void set_initial_neutron_population (double p) {initial_neutron_population = p;}
      void set_TH_pointer (HeatConduction<dim> *th) {TH = th;}
      void set_all_dtnext (const double dt_in);
      void set_all_dtnext (const vector<double> &dt_in);
      void set_all_advance (const unsigned int in);
      void set_all_advance (const vector<unsigned int> &in);
      void set_adaptive_err_tol (const double in_);
            
      // Accessors
      inline double     get_n_groups () {return parameters.n_groups;};
      inline double     get_n_prec () {return parameters.n_prec;};
      double            get_max_time ();
      double            get_min_time ();
      double            get_max_next_time ();
      double            get_min_next_time ();
      double            get_max_dt ();
      double            get_min_dt ();
      double            get_global_frequency () {return global_frequency;};
      vector<double>    get_all_time ();
      vector<double>    get_all_next_time ();
      vector<double>    get_all_dtest ();
      vector<double>    get_energy_group_frequencies () {return energy_group_frequencies;};
      vector<double>    get_min_energy_group_frequencies () {return min_energy_group_frequencies;};
      
      // Heat Source Manipulations
      Vector<double>    get_heat_src () {return heat_src;}
      void              update_heat_src (const double t, const unsigned int flag);
      void              linearize_heat_src ();
      double            get_total_neutron_population (const double t, const unsigned int flag);
      double            get_total_fission_source () const;
    
    public:
      DoFHandler<dim>   dof_handler;
      Vector<double>    heat_src;
    
    private:
      
      double get_total_energy (const double t, const unsigned int flag);
      double get_total_energy_source (double t);
      
      Parameters      parameters;
      MaterialData    material_data;
      FE_Q<dim>       fe;

      double    k_eff;
      double    initial_neutron_population;

      vector<EnergyGroup<dim>*> energy_groups;
    
    private: // Physics Variables
      double            power_time;
      double            power_time0;
      double            power;
      double            power0;
      double            global_frequency;
      vector<double>    energy_group_frequencies;
      vector<double>    min_energy_group_frequencies;
      
      unsigned int              temp_iter_counter;
      vector<double>            all_update_times;
      vector<double>            all_update_powers;
      vector<unsigned int>      all_iteration_counter;
      vector<vector<bool>>      all_update_egroups;
      vector<double>            all_update_global_frequency;
      vector<vector<double>>    all_update_group_frequencies;
      vector<vector<double>>    all_update_min_group_frequencies;
      
    private: // Precursor Parameters and Kinetics Mesh
      unsigned int          prec_count;
      double                prec_time;
      double                prec_time_old;
      Triangulation<dim>    mesh;
      
      ConstraintMatrix          hanging_node_constraints;
      vector<Vector<double>>    precursors;
      vector<Vector<double>>    precursors_old;
      map<unsigned int,double>  boundary_values;
      
      
    private: // TH Properties
      HeatConduction<dim>*   TH;
      
  };

//======================================================================================================================
  // Parameters Constructor
template <int dim>
Neutronics<dim>::Parameters::Parameters ()  {}

//======================================================================================================================
// get_parameters
template <int dim>
void Neutronics<dim>::Parameters::get_parameters (ParameterHandler &prm)
{
    n_groups                = prm.get_integer ("n_erg");
    n_prec                  = prm.get_integer ("n_prec");
    coupling                = prm.get_integer ("Kinetics Coupling");
    extrap                  = prm.get_integer ("Extrapolation Method");
    phys_type               = prm.get_integer ("Physics Type");
    tolerance               = prm.get_double ("tolerance");
    initial_power           = prm.get_double ("Initial Power");
    use_adjoint_weighting   = prm.get_bool ("Adjoint Weighting");
    
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
}

//======================================================================================================================
// Neutronics Constructor
template <int dim>
Neutronics<dim>::Neutronics (const std::string input_path, const unsigned int FE_deg)        
    :
    material_data (Constants::input_str + input_path + "/"),
    fe (FE_deg)
{
    parameters.fe_degree = FE_deg;
    parameters.project_name = input_path;
    parameters.folder_name = Constants::input_str + input_path + "/";
    parameters.prob_type = Constants::steady;
    set_parameters ();
}

//======================================================================================================================
// Neutronics Destructor
template <int dim>
Neutronics<dim>::~Neutronics ()
{
    dof_handler.clear ();
    mesh.clear ();
    TH = NULL;
    
    for (unsigned int group=0; group<energy_groups.size(); ++group)
      delete energy_groups[group];
    
    energy_groups.resize (0);
}

//======================================================================================================================
// Set all kinetics parameters
template <int dim>
void Neutronics<dim>::set_parameters ()
{
    // Local Variables
    ParameterHandler    proj_handler;
    ParameterHandler    kin_handler;
    unsigned int        phys_type;
    
    declare_project_parameters (proj_handler, parameters.folder_name);
    
    phys_type = proj_handler.get_integer ("Physics Type");
    
    // Assign Physics Type
    switch (phys_type) 
    {
      case 0:
        phys_type = Constants::Physics_Neutronics_Only;
        break;
      case 1:
        phys_type = Constants::Physics_TH_Only;
        break;
      case 2:
        phys_type = Constants::Physics_Both;
        break;
    }
    
    if (phys_type == Constants::Physics_TH_Only) return;
    
    declare_kinetics_parameters (kin_handler, parameters.folder_name);
    
    parameters.get_parameters (kin_handler);
    parameters.problem_name = proj_handler.get ("Problem Name");
    parameters.print_output = proj_handler.get_bool ("Print Output");
    
    return;
}

//======================================================================================================================
// initialize_problem
template <int dim>
void Neutronics<dim>::initialize_problem(Triangulation<dim> &mesh_)
{
    // Reset Precursor
    precursors.resize (parameters.n_prec);
    precursors_old.resize (parameters.n_prec);
    prec_time = 0.0;
    prec_count = 0;
    
    // Copy Mesh
    dof_handler.clear ();
    mesh.clear ();
    mesh.copy_triangulation (mesh_);
    
    // Initialize DOF for Precursors
    dof_handler.initialize (mesh,fe);
    dof_handler.distribute_dofs (fe);
    
    // Remake hanging node constraints
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
    hanging_node_constraints.close ();
    
    // Initialize Heat Source Vector
    heat_src.reinit (dof_handler.n_dofs ());
    
    // Allocate Memory Space for Precursors
    for (unsigned int i = 0; i < parameters.n_prec; ++i) {
      precursors[i].reinit (dof_handler.n_dofs ());
      precursors_old[i].reinit (dof_handler.n_dofs ());
    }
    
    // Allocate Energy Groups
    energy_groups.resize (parameters.n_groups);
    for (unsigned int group=0; group<parameters.n_groups; ++group)
    {
      energy_groups[group] = new EnergyGroup<dim> (group, material_data, mesh, fe);
      energy_groups[group]->setup_linear_system ();
      energy_groups[group]->set_extrap (parameters.extrap);
      energy_groups[group]->set_prob_type (parameters.prob_type);
      energy_groups[group]->set_phys_type (parameters.phys_type);
      energy_groups[group]->set_dtnext (0.0);
      energy_groups[group]->set_folder_name (parameters.folder_name);
      energy_groups[group]->set_project_name (parameters.project_name);
      energy_groups[group]->set_problem_name (parameters.problem_name);
      energy_groups[group]->set_use_adjoint_weighting (parameters.use_adjoint_weighting);
      energy_groups[group]->set_print_output (parameters.print_output);
    }
    
    // Allocate Physics Variables
    power_time  = 0.0;
    power_time0 = 0.0;
    power0      = 0.0;
    power       = parameters.initial_power;
    energy_group_frequencies.resize (parameters.n_groups, 0.0);
    min_energy_group_frequencies.resize (parameters.n_groups, 0.0);
    all_update_times.push_back(0.0);
    all_update_powers.push_back(parameters.initial_power);
    all_update_egroups.resize (1, vector<bool>(parameters.n_groups,true));
    all_update_global_frequency.push_back (0.0);
    all_update_group_frequencies.resize (1, vector<double>(parameters.n_groups,0.0));
    all_update_min_group_frequencies.resize (1, vector<double>(parameters.n_groups,0.0));
    
    // Reset Boundary Values
    // Presently, only reflective and zero flux boundaries are working
    // ---------------------------------------------------------------
    boundary_values.clear();
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(), boundary_values);
    
    // Reset Output Folder Information
    // WARNING - THIS WILL DELETE INFORMATION THERE
    // --------------------------------------------
    std::string out_folder1 = Constants::output_str + parameters.project_name + "/";
    std::string out_folder2 = Constants::output_str + parameters.project_name + "/" + parameters.problem_name + "/";
    
    
    
    return;
}

//======================================================================================================================
// update_mesh
template <int dim>
void Neutronics<dim>::update_mesh (Triangulation<dim> &mesh_)
{
  // Local Variables
  Threads::ThreadGroup<>    threads;
  
  // Update Remaining Mesh Information
  dof_handler.clear ();
  mesh.clear ();
  mesh.copy_triangulation(mesh_);
  dof_handler.initialize (mesh,fe);
  dof_handler.distribute_dofs (fe);
  boundary_values.clear();
  VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(), boundary_values);
    
  // Update Energy Group Flux Solutions
  for (unsigned int group=0; group<parameters.n_groups; ++group)
    threads += Threads::new_thread (&EnergyGroup<dim>::update_mesh, *energy_groups[group], mesh, fe, dof_handler);
  
  threads.join_all ();
  
  return;
}


//======================================================================================================================
// get_total_fission_source
template <int dim>
double Neutronics<dim>::get_total_fission_source () const
{
    std::vector<Threads::Thread<double> > threads;
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      threads.push_back (Threads::new_thread (&EnergyGroup<dim>::get_fission_source, *energy_groups[group]));

    double fission_source = 0;
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      fission_source += threads[group].return_value ();

    return fission_source;
}
  
//======================================================================================================================
// get_total_energy
template <int dim>
double Neutronics<dim>::get_total_energy (const double t, const unsigned int flag) 
{
    double energy_source = 0;
    
    //Threads::ThreadGroup<double> threads;
    //for (unsigned int group=0; group<parameters.n_groups; ++group)
      //threads += Threads::new_thread (&EnergyGroup<dim>::get_energy_source, *energy_groups[group], t, flag);
      
    //for (unsigned int group=0; group<parameters.n_groups; ++group)
      //energy_source += threads[group].return_value (); 
      
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      energy_source += energy_groups[group]->get_energy_source (t,flag);

    return energy_source;
}

//======================================================================================================================
// Eigenvalue Calculation Subroutine
template <int dim>
void Neutronics<dim>::eigenvalue_calc () 
{
    std::cout << std::setprecision (12) << std::fixed;
    
    Timer timer;
    timer.start ();
    
    double error;
    double k_eff_old = k_eff;
    double power_integral, norm_const;
    unsigned int cycle = 1;
    unsigned int iteration = 1;

    print_eigenvalue_header ();
    
    Threads::ThreadGroup<> threads;
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      threads += Threads::new_thread (&EnergyGroup<dim>::assemble_system_matrix, *energy_groups[group]);
    threads.join_all ();
    
    // Begin Iteration Loop
    do
      {
        for (unsigned int group=0; group<parameters.n_groups; ++group)
          {
            energy_groups[group]->assemble_ingroup_rhs (ZeroFunction<dim>());

            for (unsigned int bgroup=0; bgroup<parameters.n_groups; ++bgroup)
              energy_groups[group]->assemble_cross_group_rhs (*energy_groups[bgroup]);
              
            energy_groups[group]->solve ();
          }

        k_eff = get_total_fission_source();
        error = fabs(k_eff-k_eff_old)/fabs(k_eff);
        std::cout << "   Iteration " << iteration << ":  k_eff = " << k_eff << std::endl;
        k_eff_old=k_eff;

        for (unsigned int group=0; group<parameters.n_groups; ++group)
          {
            energy_groups[group]->sol = energy_groups[group]->solsave;
            energy_groups[group]->sol /= k_eff;
          }

        ++iteration;
      }
    while((error > parameters.tolerance) && (iteration < Constants::MAX_ITERATION));
    
    all_iteration_counter.resize (1,iteration - 1);
    power_integral = get_total_energy(0.0, 0);
    norm_const = parameters.initial_power/power_integral;
    
    for (unsigned int group=0; group<parameters.n_groups; ++group)
    {
        energy_groups[group]->mult_solsave(norm_const);
        energy_groups[group]->set_sol(energy_groups[group]->get_solsave ());
        energy_groups[group]->set_keff (k_eff);
        energy_groups[group]->output_results ();
    }
    
    // Compute Steady-State Precursor Concentrations
    update_precursors ();
    output_precursors ();
    
    //Compute Steady-State Heat Source
    update_heat_src (0.0, 0);
    
    std::cout << std::endl;
    std::cout << "   Finished Kinetics Eigenvalue Calculation,  k_eff=" << k_eff << ", time=" << timer() << std::endl;
    std::cout << std::endl << std::endl;
    
    // Adjoint Calculation if Required
    if (parameters.use_adjoint_weighting)
      adjoint_calc ();
    
}

//======================================================================================================================
// Adjoint Calculation Subroutine
template<int dim>
void Neutronics<dim>::adjoint_calc ()
{
    std::cout << std::setprecision (12) << std::fixed;
    
    // Timers
    Timer timer;
    timer.start ();
    
    // Print Header
    print_neutronics_adjoint_header ();
    
    // Assemble LHS Identically to Forward Calculation
    Threads::ThreadGroup<> threads;
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      threads += Threads::new_thread (&EnergyGroup<dim>::assemble_adjoint_system_matrix, *energy_groups[group]);
    threads.join_all ();
    
    // Begin Iteration Loop
    double error, err, err0;
    unsigned int iteration = 1;
    do
    {
        
        for (unsigned int group=0; group<parameters.n_groups; ++group)
        {
          energy_groups[group]->assemble_adjoint_ingroup_rhs (ZeroFunction<dim>());

          for (unsigned int bgroup=0; bgroup<parameters.n_groups; ++bgroup)
            energy_groups[group]->assemble_adjoint_cross_group_rhs (*energy_groups[bgroup]);
              
          energy_groups[group]->solve ();
        }
        
        err = get_total_fission_source();
        error = fabs(err - err0)/fabs(err);
        err0 = err;
        
        ++iteration;
    }
    while((error > parameters.tolerance) && (iteration < Constants::MAX_ITERATION));
    
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      energy_groups[group]->set_adjoint (energy_groups[group]->solsave);
    
    std::cout << std::endl;
    std::cout << "   Finished Kinetics Adjoint Calculation,  # Iterations = " << iteration-1 << ", time = " << timer() << std::endl;
    std::cout << std::endl << std::endl;
    
    
    
    return;
}

//======================================================================================================================
// Time-Step Calculation Subroutine
template<int dim>
void Neutronics<dim>::run_time_step () 
{
    std::cout << std::setprecision (12) << std::fixed;
    
    Timer timer;
    timer.start ();
    
    bool            synch (true);
    double          error, err, err0;
    unsigned int    ngrp (0), iteration (1);
    vector<double>  next_times (parameters.n_groups, 0.0);
    vector<double>  times (parameters.n_groups, 0.0);
    
    Threads::ThreadGroup<> threads;
    
    std::string advance_e_groups = "   Advancing Energy Groups:         ";
    std::string advance_ncount   = "   Advancing EGroup Step Count:     ";
    std::string advance_dt       = "   Advancing EGroup Dts:            ";
    std::string current_times    = "   Advancing EGroup Current Times:  ";
    std::string advance_times    = "   Advancing EGroup Times:          ";
    std::string iter_line;
    
    // Determine Advancing Energy Groups
    for (unsigned int i=0; i< parameters.n_groups; ++i)
    {
        if (energy_groups[i]->get_advance () == 1)
        {
            ++ngrp;
            advance_e_groups += Utilities::int_to_string (i+1, Utilities::needed_digits(i+1)) + "   ";
            advance_ncount += Utilities::int_to_string (energy_groups[i]->get_ncount ()+1, Utilities::needed_digits(energy_groups[i]->get_ncount ()+1)) + "   ";
            
            std::ostringstream strs1, strs2, strs3;
            times.at(i) = energy_groups[i]->get_time ();
            next_times.at(i) = energy_groups[i]->get_next_time ();
            strs1 << next_times.at(i) << "   ";
            strs2 << energy_groups[i]->get_dtnext () << "   ";
            strs3 << times.at(i) << "   ";
            current_times += strs3.str ();
            advance_times += strs1.str ();
            advance_dt += strs2.str ();
            
            energy_groups[i]->set_advance (0);
            energy_groups[i]->set_solsave (energy_groups[i]->get_flux (next_times.at(i)));
            energy_groups[i]->set_advance (1);
            
            threads += Threads::new_thread (&EnergyGroup<dim>::assemble_system_matrix, *energy_groups[i]);
        }
        else
        {
            synch = false;
        }
        energy_groups[i]->set_prec_time (prec_time);
        energy_groups[i]->set_prob_type (Constants::transient);
    }
    
    threads.join_all ();
    print_neutronics_time_step_header ();
    cout << std::endl << advance_e_groups << std::endl;
    cout << advance_ncount << std::endl;
    cout << current_times << std::endl;
    cout << advance_dt << std::endl;
    cout << advance_times << std::endl << std::endl;
    
    // Perform Solution Iterations
    do 
    {
        err = 0;
        
        for (unsigned int i = 0; i < parameters.n_groups; ++i)
        {
            if (energy_groups[i]->get_advance () == 1)
            {
                energy_groups[i]->assemble_ingroup_rhs (ZeroFunction<dim>());
                energy_groups[i]->assemble_precursor_rhs (precursors, dof_handler);
                for (unsigned int bgroup=0; bgroup<parameters.n_groups; ++bgroup)
                  energy_groups[i]->assemble_cross_group_rhs (*energy_groups[bgroup]);
                energy_groups[i]->solve ();
                err += energy_groups[i]->get_fission_source ();
            }
        }
        
        error = fabs(err - err0)/fabs(err);
        err0 = err;
        
        std::cout << "   Iteration " << iteration << ":  error = " << error << std::endl;
        
        ++iteration;
    }
    while((error > parameters.tolerance) && (iteration < Constants::MAX_ITERATION)); 
    
    temp_iter_counter = iteration - 1;
    
    std::cout << std::endl;
    std::cout << "   Finished Neutronics Time Step,  # Iterations = " << iteration-1 << ", time = " << timer() << std::endl;
    std::cout << std::endl << std::endl;
    
    return;
}

//======================================================================================================================
// Time-Step Calculation Subroutine
template<int dim>
void Neutronics<dim>::run_PCQS_time_step () 
{
    std::cout << std::setprecision (12) << std::fixed;
    
    // Timers
    Timer timer;
    timer.start ();
    
    // Local Variables
    double          error, err, err0;
    unsigned int    iteration (1);
    
    // Assemble System Matrices
    Threads::ThreadGroup<> threads1;
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      threads1 += Threads::new_thread (&EnergyGroup<dim>::assemble_system_matrix, *energy_groups[group]);
    threads1.join_all ();
    
    // Perform Iterations
    do
    {
      
      for (unsigned int g1=0; g1<parameters.n_groups; ++g1)
      {
        energy_groups[g1]->assemble_ingroup_rhs ( ZeroFunction<dim>() );
        energy_groups[g1]->assemble_precursor_rhs( precursors,dof_handler );
        for (unsigned int g2=0; g2<parameters.n_groups; ++g2)
          energy_groups[g1]->assemble_cross_group_rhs (*energy_groups[g2]);
        
        energy_groups[g1]->solve ();
      }
            
      err = get_total_fission_source();
      error = fabs(err - err0)/fabs(err);
      err0 = err;
      std::cout << "   Iteration " << iteration << ":  error = " << error << std::endl;
      
      ++iteration;
      
    }
    while((error > parameters.tolerance) && (iteration < Constants::MAX_ITERATION)); 
    
    temp_iter_counter = iteration - 1;
    
    // update amplitude
    update_PCQS_amplitude ();
    
    std::cout << std::endl;
    std::cout << "   Finished PCQS Neutronics Time Step,  # Iterations = " << iteration-1 << ", time = " << timer() << std::endl;
    std::cout << std::endl << std::endl;
    
    return;
}

//======================================================================================================================
// update PCQS amplitudes
template <int dim>
void Neutronics<dim>::update_PCQS_amplitude ()
{
  // Local Variables
  unsigned int      n_prec (parameters.n_prec + 1);
  double            current_pop, normalization;
  
  // PRKE Variables
  Table<2, double>          reactivity;
  vector<Table<2,double>>   beta;
  vector<double>            mean_gen_time;
  vector<double>            lambda;
  vector<vector<double>>    precurs_vals;
  vector<double>            amp;
  
  // Allocate PRKE Variables
  if (parameters.PRKE_amplitude_count == Constants::amplitude_count_single)
  {
    reactivity.reinit (1,1);
    beta.resize (n_prec);
    mean_gen_time.resize (1, 0.0);
    lambda.resize (n_prec);
    precurs_vals.resize (n_prec, vector<double> (1, 0.0));
  }
  else
  {
    reactivity.reinit (parameters.n_groups,parameters.n_groups);
    beta.resize (n_prec);
    mean_gen_time.resize (parameters.n_groups);
    lambda.resize (n_prec);
    precurs_vals.resize (parameters.n_groups, vector<double> (n_prec, 0.0));
  }
  
  // Get Energy Group Values
  for (unsigned int i = 0; i < n_prec; ++i)
  {
    lambda.at(i) = material_data.get_lambda(i);
  }
  
  // Normalize Solution
  current_pop = get_total_neutron_population (get_max_next_time (), 0);
  normalization = current_pop / initial_neutron_population;
  for (unsigned int g = 0; g < parameters.n_groups; ++g)
    energy_groups[g]->mult_solsave (normalization);
  
  // get PRKE parameters
  if (parameters.PRKE_amplitude_count == Constants::amplitude_count_single)
  {
    mean_gen_time.at(0) = get_total_neutron_population (get_max_next_time (), 0);
    for (unsigned int g = 0; g < parameters.n_groups; ++g)
    {
      
    }
  }
  else
  {
    for (unsigned int g = 0; g < parameters.n_groups; ++g)
    {
      mean_gen_time.at(g) = energy_groups[g]->get_neutron_population(get_max_next_time (), 0);
    }
  }
  
  return;
}

//======================================================================================================================
// Update Precursor Solution Subroutine
template<int dim>
void Neutronics<dim>::update_precursors () 
{
    
    // Local Variables
    const QGauss<dim>     quadrature_formula(fe.degree + 1);
    FEValues<dim>         fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values);
    SparsityPattern       sparsity_pattern;
    SparseMatrix<double>  matrix;
    Vector<double>        rhs;
    
    const unsigned int    n_dofs = dof_handler.n_dofs ();
    const unsigned int    nprec = parameters.n_prec;
    const unsigned int    dofs_per_cell = fe.dofs_per_cell;
    const unsigned int    n_q_points    = quadrature_formula.size();
    vector<double>        prec_frac (nprec, 0.0);
    vector<double>        prec_values (n_q_points, 0.0);
    vector<double>        temp_values (n_q_points, 1.0);
        
    vector<vector<double>> fiss_XS (parameters.n_groups, vector<double> (n_q_points, 0.0));
    vector<vector<double>> sol_values (parameters.n_groups, vector<double> (n_q_points, 0.0));
    
    FullMatrix<double>    cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>        cell_rhs (dofs_per_cell);
    vector<unsigned int>  local_dof_indices (dofs_per_cell);
    Vector<double>        temperature;
    
    // Local Transient Variables
    double t_min = get_min_time ();
    double dt_prec = fabs(t_min - prec_time);
    if (dt_prec < Constants::SMALL && parameters.prob_type == Constants::transient) return;
    
    if (parameters.phys_type == Constants::Physics_Both) 
      temperature.reinit (TH->get_temperature (t_min));
    
    sparsity_pattern.reinit (n_dofs, n_dofs, dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    hanging_node_constraints.condense (sparsity_pattern);
    sparsity_pattern.compress ();
    
    // Precursor Fraction for the Fission Term
    for (unsigned int i = 0; i < nprec; ++i)
      if (parameters.prob_type == Constants::transient)
        prec_frac.at(i) = material_data.get_beta(i) * dt_prec;
      else
        prec_frac.at(i) = material_data.get_beta(i) / material_data.get_lambda(i);
    
    // Update Old Solutions for Transient Steps
    if (parameters.prob_type == Constants::transient)
      for (unsigned int q_prec = 0; q_prec < parameters.n_prec; ++q_prec)
        precursors_old[q_prec] = precursors[q_prec];
    
    // Loop through Precursor Groups
    for (unsigned int p = 0; p < nprec; ++p)
    {
        matrix.reinit (sparsity_pattern);
        rhs.reinit (n_dofs);
        
        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for (; cell!=endc; ++cell)
        {
          cell_matrix = 0;
          cell_rhs = 0;
          fe_values.reinit (cell);
          cell->get_dof_indices (local_dof_indices);
          
          if (parameters.phys_type == Constants::Physics_Both) 
            fe_values.get_function_values (temperature, temp_values);
          
          if (parameters.prob_type == Constants::transient)
            fe_values.get_function_values (precursors_old[p], prec_values);
          for (unsigned int j=0; j<parameters.n_groups; ++j)
          {
            fe_values.get_function_values (energy_groups[j]->get_flux(t_min), sol_values[j]);
            for (unsigned int k=0; k<n_q_points; ++k)
              if (parameters.phys_type == Constants::Physics_Both) 
                fiss_XS[j][k] = material_data.get_fission_XS (j,cell->material_id(), temp_values.at(k));
              else
                fiss_XS[j][k] = material_data.get_fission_XS (j,cell->material_id(), temp_values.at(k));
          }
          
          // Assemble local matrix and rhs
          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i) 
            {
              if (parameters.prob_type == Constants::transient)
                cell_rhs(i) += prec_values.at(q_point) * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
              for (unsigned int group = 0; group<parameters.n_groups; ++group)
                cell_rhs(i) += prec_frac.at(p) * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point) * fiss_XS[group][q_point] * sol_values[group][q_point];
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                if (parameters.prob_type == Constants::transient)
                  cell_matrix(i,j) += (1 + dt_prec * material_data.get_lambda(i)) * fe_values.shape_value(i,q_point) * fe_values.shape_value(j,q_point) * fe_values.JxW(q_point);
                else
                  cell_matrix(i,j) += fe_values.shape_value(i,q_point) * fe_values.shape_value(j,q_point) * fe_values.JxW(q_point);
            }
          
          for (unsigned int i=0; i<dofs_per_cell; ++i) {
            rhs(local_dof_indices[i]) += cell_rhs(i);
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i,j));
          }
          
        }
        
        
        // Solve Precursor Equation
        hanging_node_constraints.condense (rhs);
        MatrixTools::apply_boundary_values (boundary_values, matrix, precursors[p], rhs);
        
        SolverControl   solver_control (matrix.m(), 1e-12*rhs.l2_norm());
        SolverCG<>      cg (solver_control);
        
        PreconditionSSOR<> preconditioner;
        preconditioner.initialize(matrix, 1.2);
        
        cg.solve (matrix, precursors[p], rhs, preconditioner);
        hanging_node_constraints.distribute (precursors[p]);
        
    }
    
    if (parameters.prob_type == Constants::steady)
      for (unsigned int p = 0; p < nprec; ++p)
        precursors_old[p] = precursors[p];
    if (parameters.prob_type == Constants::transient)
    {
        prec_time_old = prec_time;
        prec_time += dt_prec;
        ++prec_count;
        output_precursors ();
    }
    
    return;
}

//======================================================================================================================
// Update Heat Source Subroutine
template<int dim>
void Neutronics<dim>::update_heat_src (const double t, const unsigned int flag) 
{
  
  // Local Variables
  const QGauss<dim>     quadrature_formula(fe.degree + 1);
  FEValues<dim>         fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values);
  SparsityPattern       sparsity_pattern;
  SparseMatrix<double>  matrix;
  Vector<double>        rhs;
  Vector<double>        temperature;
  
  const unsigned int    n_dofs        =   dof_handler.n_dofs ();
  const unsigned int    dofs_per_cell =   fe.dofs_per_cell;
  const unsigned int    n_q_points    =   quadrature_formula.size();
  const double          energy_conversion (Constants::ENERGY_FROM_FISS*Constants::MEV_TO_J);
  
  vector<vector<double>> fiss_XS (parameters.n_groups, vector<double> (n_q_points, 0.0));
  vector<vector<double>> sol_values (parameters.n_groups, vector<double> (n_q_points, 0.0));
  std::vector<double>    temp_vals (n_q_points, 1.0);
  
  FullMatrix<double>    cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>        cell_rhs (dofs_per_cell);
  vector<unsigned int>  local_dof_indices (dofs_per_cell);
  
  sparsity_pattern.reinit (n_dofs, n_dofs, dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  hanging_node_constraints.condense (sparsity_pattern);
  sparsity_pattern.compress ();
  
  if (parameters.phys_type == Constants::Physics_Both)
    temperature.reinit (TH->get_temperature (get_max_time ()));
  
  matrix.reinit (sparsity_pattern);
  rhs.reinit (n_dofs);
        
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    
    cell_matrix = 0;
    cell_rhs = 0;
    fe_values.reinit (cell);
    cell->get_dof_indices (local_dof_indices);
    
    if (parameters.phys_type == Constants::Physics_Both) 
      fe_values.get_function_values (temperature, temp_vals);
    
    for (unsigned int group = 0; group<parameters.n_groups; ++group)
    {
      //cout << energy_groups[group]->solsave.size() << std::endl;
      if (parameters.prob_type == Constants::steady) 
        fe_values.get_function_values (energy_groups[group]->solsave, sol_values[group]);
      else if (parameters.prob_type == Constants::transient && flag == 0)
        fe_values.get_function_values (energy_groups[group]->solsave, sol_values[group]);
      else if (parameters.prob_type == Constants::transient && flag != 0)
        fe_values.get_function_values (energy_groups[group]->get_flux(t), sol_values[group]);
      
      for (unsigned int q_node = 0; q_node < n_q_points; ++q_node)
        if (parameters.phys_type == Constants::Physics_Both)
          fiss_XS[group][q_node] = material_data.get_fission_XS (group,cell->material_id(), temp_vals.at(q_node)) / material_data.get_nu (group,cell->material_id());
        else
          fiss_XS[group][q_node] = material_data.get_fission_XS (group,cell->material_id(), temp_vals.at(q_node)) / material_data.get_nu (group,cell->material_id());
    }
    
    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int group = 0; group<parameters.n_groups; ++group)
          cell_rhs(i) +=  sol_values[group][q_point] * fiss_XS[group][q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point)* energy_conversion;
    
    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          cell_matrix(i,j) += fe_values.shape_value(i,q_point) * fe_values.shape_value(j,q_point) * fe_values.JxW(q_point);
    
    for (unsigned int i=0; i<dofs_per_cell; ++i)
    {
      rhs(local_dof_indices[i]) += cell_rhs(i);
      for (unsigned int j=0; j<dofs_per_cell; ++j)
        matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i,j));
    }
    
  }
  
  // Solve Heat Source Equation
  hanging_node_constraints.condense (rhs);
  MatrixTools::apply_boundary_values (boundary_values, matrix, heat_src, rhs);
  
  SolverControl   solver_control (matrix.m(), 1e-12*rhs.l2_norm());
  SolverCG<>      cg (solver_control);
  
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(matrix, 1.2);
  
  cg.solve (matrix, heat_src, rhs, preconditioner);
  hanging_node_constraints.distribute (heat_src);
  
  return;
}

//======================================================================================================================
// Update Heat Source Subroutine
template<int dim>
void Neutronics<dim>::linearize_heat_src () 
{
  
  // Local Variables
  
  
  
  return;
}

//======================================================================================================================
// Update Solution Subroutine
template<int dim>
void Neutronics<dim>::update_solution () 
{
  // Local Variables
  double energy (0.0);
  
  if (parameters.prob_type == Constants::steady)
    for (unsigned int group=0; group<parameters.n_groups;++group)
      energy_groups[group]->update_solution ();
  else
  {
    for (unsigned int group=0; group<parameters.n_groups;++group)
      if (energy_groups[group]->get_advance () == 1)
        energy_groups[group]->update_solution ();
    
    unsigned int n = all_update_egroups.size ();
    all_iteration_counter.push_back (temp_iter_counter);
    all_update_egroups.resize(n + 1);
    all_update_egroups.back().resize (parameters.n_groups, false);
    for (unsigned int g = 0; g < parameters.n_groups; ++g)
      all_update_egroups.back().at(g) = energy_groups[g]->get_advance ();
    
    // Update Power Information
    power0 = power;
    power_time0 = power_time;
    power_time = get_max_time ();
    set_all_advance (0);
    update_precursors ();
    power = get_total_energy (power_time, 1);
    compile_kinetics_physics (power_time);
    
    all_update_times.push_back (power_time);
    all_update_powers.push_back (power);
    
    all_update_global_frequency.push_back (global_frequency);
    all_update_group_frequencies.resize (n + 1);
    all_update_min_group_frequencies.resize (n + 1);
    all_update_group_frequencies.back().resize (parameters.n_groups, 0.0);
    all_update_min_group_frequencies.back().resize (parameters.n_groups, 0.0);
    
    for (unsigned int g = 0; g < parameters.n_groups; ++g)
    {
      all_update_group_frequencies.back().at(g) = energy_group_frequencies.at(g);
      all_update_min_group_frequencies.back().at(g) = min_energy_group_frequencies.at(g);
    }
  }
  
  return;
}

//======================================================================================================================
// Update Powers Subroutine
template<int dim>
void Neutronics<dim>::update_powers (const double t, const unsigned int flag) 
{
  // Not sure anymore why this is here...
  
  return;
}

//======================================================================================================================
// Update Heat Source Subroutine
template<int dim>
double Neutronics<dim>::get_total_neutron_population (const double t, const unsigned int flag) 
{
  
  // Local Variables
  double    neutron_pop;
  
  for (unsigned int group=0; group<parameters.n_groups; ++group)
    neutron_pop += energy_groups[group]->get_neutron_population (t, flag);
  
  return neutron_pop;
}

//======================================================================================================================
// Compile Physics Quantities Subroutine
// This is performed after a time step
// and solution update.
template<int dim>
void Neutronics<dim>::compile_kinetics_physics (const double t)
{
  // Local Variables
  double    dt (0.0);
  
  global_frequency = 0.0;
  energy_group_frequencies.resize (parameters.n_groups, 0.0);
  min_energy_group_frequencies.resize (parameters.n_groups, 1.0e32);
  
  for (unsigned int g = 0; g < parameters.n_groups; ++g)
  {
    
    // Work with Local Frequency Values for each Energy Group
    vector<double>  periods (energy_groups[g]->get_local_frequencies (t));
    for (unsigned int i = 0; i < periods.size (); ++i)
    {
      if (fabs(periods.at(i)) > Constants::SMALL && fabs(periods.at(i)) < fabs(min_energy_group_frequencies.at(g)))
        min_energy_group_frequencies.at(g) = periods.at(i);
    }
    
    // Compile Energy Group Frequencies
    if (t > energy_groups[g]->get_time0 ())
    {
      double  tp  = energy_groups[g]->get_total_power ();
      double  tp0 = energy_groups[g]->get_total_power0 ();
      double  dtt = energy_groups[g]->get_dt ();
      if (fabs(dtt) > Constants::SMALL && fabs(tp0) > Constants::SMALL)
        energy_group_frequencies[g] = log(tp/tp0)/dtt;
      else
        energy_group_frequencies[g] = 0.0;
    }
    else if (t >= energy_groups[g]->get_time1 () && t <= energy_groups[g]->get_time0 ())
    {
      double  tp  = energy_groups[g]->get_total_power0 ();
      double  tp0 = energy_groups[g]->get_total_power1 ();
      double  dtt = energy_groups[g]->get_dt0 ();
      if (fabs(dtt) > Constants::SMALL && fabs(tp0) > Constants::SMALL)
        energy_group_frequencies[g] = log(tp/tp0)/dtt;
      else
        energy_group_frequencies[g] = 0.0;
    }
  }
  
  // Get Global Frequency Based on Total Reactor Power
  dt = power_time - power_time0;
  if (t > power_time0)
    if (fabs(dt) > Constants::SMALL && fabs(power0) > Constants::SMALL)
      global_frequency = log(power/power0)/dt;
    else
      global_frequency = 0.0;
      
  return;
}

//======================================================================================================================
// Output Precursor Solution Subroutine
template<int dim>
void Neutronics<dim>::output_precursors () 
{
    if (parameters.prob_type == Constants::transient && !parameters.print_output) return;
    
    // Local Variables
    std::string fname1, fname2;
    
    if (parameters.prob_type == Constants::steady)
    {
        fname1 = Constants::output_str + parameters.project_name + "/" + parameters.problem_name + "/" + "precursor_mesh.eps";
        fname2 = Constants::output_str + parameters.project_name + "/" + parameters.problem_name + "/" + "precursor_sol_0_";
        std::ofstream mesh_output (fname1);
        
        GridOut grid_out;
        grid_out.write_eps (mesh, mesh_output);
    }
    else
    {
        fname2 = Constants::output_str + parameters.project_name + "/" + parameters.problem_name + "/" + 
                "precursor_sol_" + Utilities::int_to_string(prec_count,Utilities::needed_digits(prec_count)) + "_";
    }
    
    for (unsigned int prec=0; prec < parameters.n_prec; ++prec)
    {
        std::string filename = fname2 + Utilities::int_to_string(prec+1,Utilities::needed_digits(prec+1)) + ".gmv";
        std::ofstream sol_output (filename);
        
        DataOut<dim> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (precursors[prec], "precursor " + Utilities::int_to_string(prec+1,Utilities::needed_digits(prec+1)));
        data_out.build_patches ();
        data_out.write_gnuplot (sol_output);
    }
    
    return;
}

//======================================================================================================================
// Maximum Time Calculation Subroutine
template<int dim>
double Neutronics<dim>::get_max_time () 
{
  
  double tempt;
  double ttime (0.0);
  unsigned int nerg = parameters.n_groups;
  
  for (unsigned int i=0;i<nerg;++i)
  {
    tempt = energy_groups[i]->get_time ();
    if (tempt > ttime) {ttime = tempt;}
  }
  
  return ttime;
}

//======================================================================================================================
// Minimum Time Calculation Subroutine
template<int dim>
double Neutronics<dim>::get_min_time () 
{
  
  double tempt;
  double ttime = Constants::LARGE;
  unsigned int nerg = parameters.n_groups;
  
  for (unsigned int i=0;i<nerg;++i)
  {
    tempt = energy_groups[i]->get_time ();
    if (tempt < ttime) {ttime = tempt;}
  }
  
  return ttime;
}

//======================================================================================================================
// Maximum Next Time Calculation Subroutine
template<int dim>
double Neutronics<dim>::get_max_next_time () 
{
  
  double tempt;
  double ttime (0.0);
  unsigned int nerg = parameters.n_groups;
  
  for (unsigned int i=0;i<nerg;++i)
  {
    tempt = energy_groups[i]->get_next_time ();
    if (tempt > ttime) {ttime = tempt;}
  }
  
  return ttime;
}

//======================================================================================================================
// Minimum Next Time Calculation Subroutine
template<int dim>
double Neutronics<dim>::get_min_next_time () 
{
  
  double tempt;
  double ttime = Constants::LARGE;
  unsigned int nerg = parameters.n_groups;
  
  for (unsigned int i=0;i<nerg;++i)
  {
    tempt = energy_groups[i]->get_next_time ();
    if (tempt < ttime) {ttime = tempt;}
  }
  
  return ttime;
}

//======================================================================================================================
// Maximum Dt Calculation Subroutine
template<int dim>
double Neutronics<dim>::get_max_dt () 
{
  
  double tempt;
  double dtt (0.0);
  unsigned int nerg = parameters.n_groups;
  
  for (unsigned int i=0;i<nerg;++i)
  {
    tempt = energy_groups[i]->dtnext;
    if (tempt > dtt) {dtt = tempt;}
  }
  
  
  return dtt;
}

//======================================================================================================================
// Minimum Dt Calculation Subroutine
template<int dim>
double Neutronics<dim>::get_min_dt () 
{
  
  double tempt;
  double dtt = Constants::LARGE;
  unsigned int nerg = parameters.n_groups;
  
  for (unsigned int i=0;i<nerg;++i)
  {
    tempt = energy_groups[i]->dtnext;
    if (tempt < dtt) {dtt = tempt;}
  }
  
  
  return dtt;
}

//======================================================================================================================
// Set All Energy Group Advance Flags Subroutine
template<int dim>
void Neutronics<dim>::set_all_advance (const unsigned int int_in) 
{
  
  for (unsigned int group=0; group < parameters.n_groups; ++group)
    energy_groups[group]->set_advance (int_in);
  
  return;
}

//======================================================================================================================
// Set All Energy Group dtnext Variables Subroutine
template<int dim>
void Neutronics<dim>::set_all_advance (const vector<unsigned int> &in) 
{
  for(unsigned int group=0; group<parameters.n_groups; ++group)
    energy_groups[group]->set_advance(in[group]);
  
  return;
}

//======================================================================================================================
// Set All Energy Group dtnext Variables Subroutine
template<int dim>
void Neutronics<dim>::set_all_dtnext (const double dt_in) 
{
  
  for (unsigned int group=0; group < parameters.n_groups; ++group)
    energy_groups[group]->set_dtnext (dt_in);
  
  return;
}

//======================================================================================================================
// Set All Energy Group dtnext Variables Subroutine
template<int dim>
void Neutronics<dim>::set_all_dtnext (const vector<double> &dt_in) 
{
  
  for (unsigned int group=0; group < parameters.n_groups; ++group)
    energy_groups[group]->set_dtnext (dt_in.at(group));
  
  return;
}

//======================================================================================================================
// Get All Energy Group time Variables Subroutine
template<int dim>
vector<double> Neutronics<dim>::get_all_time () 
{
  // Local Variables
  vector<double>    times (parameters.n_groups, 0.0);
  
  for (unsigned int group = 0; group < parameters.n_groups; ++group)
    times.at(group) = energy_groups[group]->get_time ();
  
  return times;
}

//======================================================================================================================
// Get All Energy Group dtnext Variables Subroutine
template<int dim>
vector<double> Neutronics<dim>::get_all_next_time () 
{
  // Local Variables
  vector<double>    times (parameters.n_groups, 0.0);
  
  for (unsigned int group = 0; group < parameters.n_groups; ++group)
    times.at(group) = energy_groups[group]->get_next_time ();
  
  return times;
}

//======================================================================================================================
// Get All Energy Group dtest Variables Subroutine
template<int dim>
vector<double> Neutronics<dim>::get_all_dtest () 
{
  // Local Variables
  vector<double>    times (parameters.n_groups, 0.0);
  
  for (unsigned int group = 0; group < parameters.n_groups; ++group)
    times.at(group) = energy_groups[group]->get_dtest ();
  
  return times;
}

//======================================================================================================================
// Set All Energy Group dtnext Variables Subroutine
template<int dim>
void Neutronics<dim>::set_prob_type (const unsigned int in_) 
{
  parameters.prob_type = in_;
  for (unsigned int group=0; group < parameters.n_groups; ++group)
    energy_groups[group]->set_prob_type (in_);
  
  return;
}

//======================================================================================================================
// 
template<int dim>
void Neutronics<dim>::set_time_step_selection (const unsigned int in_) 
{
  parameters.time_step_selection = in_;
  for (unsigned int group=0; group < parameters.n_groups; ++group)
    energy_groups[group]->set_time_step_selection (in_);
  
  return;
}

//======================================================================================================================
// 
template<int dim>
void Neutronics<dim>::set_dt_estimator (const unsigned int in_) 
{
  parameters.dt_estimator = in_;
  for (unsigned int group=0; group < parameters.n_groups; ++group)
    energy_groups[group]->set_dt_estimator (in_);
  
  return;
}

//======================================================================================================================
// 
template<int dim>
void Neutronics<dim>::set_adaptive_err_tol (const double in_) 
{
  parameters.adaptive_err_tol = in_;
  for (unsigned int group=0; group < parameters.n_groups; ++group)
    energy_groups[group]->set_adaptive_err_tol (in_);
  
  return;
}

//======================================================================================================================
// Terminate Problem Subroutine
template<int dim>
void Neutronics<dim>::terminate_problem ()
{
  // Local Variables
  unsigned int err;
  std::string filename1 = Constants::output_str + parameters.project_name + "/" + parameters.problem_name + "/kinetics_power_data";
  std::string filename2 = Constants::output_str + parameters.project_name + "/" + parameters.problem_name + "/kinetics_groups_data";
  std::string filename3 = Constants::output_str + parameters.project_name + "/" + parameters.problem_name + "/kinetics_physics_data";
  char * strs1 = new char[filename1.length() + 1];
  char * strs2 = new char[filename2.length() + 1];
  char * strs3 = new char[filename3.length() + 1];
  std::strcpy (strs1,filename1.c_str());
  std::strcpy (strs2,filename2.c_str());
  std::strcpy (strs3,filename3.c_str());
  std::ifstream infile1 (filename1);
  std::ifstream infile2 (filename2);
  std::ifstream infile3 (filename3);
  if (infile1.good ())
    err = remove (strs1);
  if (infile2.good ())
    err = remove (strs2);
  if (infile3.good ())
    err = remove (strs3);
    
  std::ofstream myfile1;
  std::ofstream myfile2;
  std::ofstream myfile3;
  myfile1.open (filename1, ios::out | ios::ate);
  myfile2.open (filename2, ios::out | ios::ate);
  myfile3.open (filename3, ios::out | ios::ate);
      
  for (unsigned int i = 0; i < all_update_times.size (); ++i)
  {
    std::ostringstream ss1;
    std::ostringstream ss2;
    std::ostringstream ss3;
    ss1 << "Step: ";
    ss1 << ss1.width(8) << std::right << i;
    ss2 << "Step: ";
    ss2 << ss1.width(8) << std::right << i;
    ss3 << "Step: ";
    ss3 << ss1.width(8) << std::right << i;
    ss1 << ",    Advance Groups: ";
    ss2 << ",    Advance Groups: ";
    ss3 << ",    Advance Groups: ";
    
    for (unsigned int g = 0; g < parameters.n_groups; ++g)
      if (all_update_egroups.at(i).at(g))
        ss1 << std::right << 1 << ";";
      else
        ss1 << std::right << 0 << ";";
        
    for (unsigned int g = 0; g < parameters.n_groups; ++g)
      if (all_update_egroups.at(i).at(g))
        ss2 << std::right << 1 << ";";
      else
        ss2 << std::right << 0 << ";";
    
    for (unsigned int g = 0; g < parameters.n_groups; ++g)
      if (all_update_egroups.at(i).at(g))
        ss3 << std::right << 1 << ";";
      else
        ss3 << std::right << 0 << ";";
        
    ss1 << "    Time: ";
    ss1 << std::setprecision (8) << std::fixed << all_update_times.at(i);
    ss3 << "    Time: ";
    ss3 << std::setprecision (8) << std::fixed << all_update_times.at(i);
    ss1 << "    #Iterations: ";
    ss1 << ss1.width(6) << std::right << all_iteration_counter.at(i);
    ss1 << ",    Power: ";
    ss1 << std::setprecision (8) << std::scientific << all_update_powers.at(i);
    ss1 << std::endl;
    ss2 << std::endl;
    
    ss3 << ",    Glob Frequency: ";
    ss3 << std::setprecision (8) << std::fixed << all_update_global_frequency.at(i);
    ss3 << ",    Group Frequencies: ";
    for (unsigned int g = 0; g < parameters.n_groups; ++g)
      ss3 << all_update_group_frequencies.at(i).at(g) << ";";
      
    ss3 << "   Min Frequencies: ";
    for (unsigned int g = 0; g < parameters.n_groups; ++g)
      ss3 << all_update_min_group_frequencies.at(i).at(g) << ";";
    
    ss3 << std::endl;
    
    myfile1 << ss1.str();
    myfile2 << ss2.str();
    myfile3 << ss3.str();
  }
  
  myfile1.close ();
  myfile2.close ();
  myfile3.close ();
  
  delete strs1;
  delete strs2;
  delete strs3;
  return;
}

//======================================================================================================================
#endif


