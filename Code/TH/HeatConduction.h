//======================================================================================================================
/*!
 * \class HeatConduction
 * \brief
 *
 * Routine to print out heading and various outputs.
 */
//======================================================================================================================

#ifndef HEATCONDUCTION_H
#define HEATCONDUCTION_H

#include "../Global/GlobalHeaders.h"
#include "../Global/Deal2Headers.h"
#include "../Global/Constants.h"

#include "THMaterialData.h"
#include "../Diffusion/Neutronics.h"


template <int dim>
class HeatConduction 
{

  public:

  // Constructor
  HeatConduction (const std::string inp_path, const unsigned int fe_deg);
  void set_parameters ();
  
  // Destructor
  ~HeatConduction ();

  // Mesh Routines
  unsigned int n_active_cells () const;
  unsigned int n_dofs () const;
  void update_mesh (Triangulation<dim> &mesh_);
  
  // Matrix and Vector Operations
  void initialize_problem (Triangulation<dim> &mesh_);
  void assemble_rhs ();
  void assemble_system_matrix ();
  
  // Solution Routines
  void steady_state_solve ();
  void run_time_step ();
  void solve ();
  void update_solution ();
  void output_results () const;
  void terminate_problem ();
  
  //Manipulators
  void set_heat_src_pointers (Vector<double> *vec_in, DoFHandler<dim> *dof_in) {heat_src=vec_in; dof_kin=dof_in;};
  void set_initial_temp_integral (double p) {initial_temp_integral = p;};
  void set_prob_type (unsigned int p) {prob_type = p;};
  void set_phys_type (unsigned int p) {phys_type = p;};
  void set_extrap (unsigned int p) {extrap = p;};
  void set_advance (unsigned int p) {advance = p;};
  void set_time_step_selection (unsigned int p) {time_step_selection = p;};
  void set_dt_estimator (unsigned int p) {dt_estimator = p;};
  
  void set_solsave (Vector<double> val) {solsave=val;};
  void set_sol (Vector<double> val) {sol=val;};
  void set_sol0 (Vector<double> val) {sol0=val;};
  void set_sol1 (Vector<double> val) {sol1=val;};
  void set_dtnext (double t) {dtnext = t;};
  void set_adaptive_err_tol (double t) {adaptive_err_tol = t;};
  
  // Accessors
  Vector<double> get_temperature (double t);
  double get_total_energy (double t);
  double get_total_energy ();
  double get_total_energy_save ();
  
  // Accessors - time variables
  double get_time () {return time;};
  double get_time0 () {return time0;};
  double get_time1 () {return time1;};
  double get_dt () {return dt;};
  double get_dt0 () {return dt0;};
  double get_dtnext () {return dtnext;};
  double get_dtest () {return dtest;};
  double get_next_time () {return time+dtnext;};
  
  // Accessors - solution variables
  Vector<double> get_solsave () {return solsave;};
  Vector<double> get_sol () {return sol;};
  Vector<double> get_sol0 () {return sol0;};
  Vector<double> get_sol1 () {return sol1;};
  
  double tot_energy, tot_energy0, tot_energy1;
  
  private:
  
  // solution vectors
  Vector<double> sol;
  Vector<double> sol0;
  Vector<double> sol1;
  Vector<double> solsave;
  Vector<double> adjoint;
  
  // time variables
  double time;
  double time0;
  double time1;
  double dt;
  double dt0;
  double dtnext;
  double dtest;
  
  // Variables
  unsigned int  prob_type;
  unsigned int  phys_type;
  unsigned int  extrap;
  unsigned int  TH_Type;
  unsigned int  advance;
  unsigned int  ncount;
  unsigned int  time_step_selection;
  unsigned int  dt_estimator;
  
  double	    sim_time;
  double        tolerance;
  double        initial_temp_integral;
  double        adaptive_err_tol;
  
  bool          temp_dep;
  bool          print_output;
  
  double        bound_temp;
  double        bulk_temp;
  double        heat_trans;
  
  std::string   folder_name;
  std::string   project_name;
  std::string   problem_name;
  
  THMaterialData        material_data;
  Triangulation<dim>    mesh;
  DoFHandler<dim>       dof_handler;
  FE_Q<dim>             fe;
  
  // Sparsity Pattern
  SparsityPattern	    sparsity_pattern;
  SparseMatrix<double>	system_matrix;

  // RHS Vectors
  Vector<double>    system_rhs;
  Vector<double>    rhs_src;

  std::map<unsigned int,double> boundary_values;
  ConstraintMatrix              hanging_node_constraints;
  
  private:
  Vector<double>*   heat_src;
  DoFHandler<dim>*  dof_kin;
  
  void      calculate_dtest ();
  double    calculate_maximum_absolute_secder_stencil ();
  double    calculate_maximum_relative_secder_stencil ();
  double    calculate_minimum_local_period ();
  
};

//======================================================================================================================
// HeatConduction Constructor
template <int dim>
HeatConduction<dim>::HeatConduction (const std::string proj_name, const unsigned int FE_deg)
                  :
                  material_data (Constants::input_str + proj_name + "/"), 
                  fe (FE_deg)
{
    project_name = proj_name;
    folder_name = Constants::input_str + proj_name + "/";
    ncount = 0;
    set_parameters ();
}

//======================================================================================================================
// HeatConduction Destructor
template <int dim>
HeatConduction<dim>::~HeatConduction ()
{
dof_handler.clear ();
mesh.clear ();

heat_src = NULL;
dof_kin = NULL;
}

//======================================================================================================================
// n_active_cells
template <int dim>
unsigned int HeatConduction<dim>::n_active_cells () const {
    return mesh.n_active_cells ();
}

//======================================================================================================================
// n_dofs
template <int dim>
unsigned int HeatConduction<dim>::n_dofs () const {
    return dof_handler.n_dofs ();
}

//======================================================================================================================
// set_parameters
template <int dim>
void HeatConduction<dim>::set_parameters () 
{
    // Local Variables
    ParameterHandler  proj_handler;
    ParameterHandler  TH_handler;

    declare_project_parameters (proj_handler, folder_name);

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
    
    if (phys_type == Constants::Physics_Neutronics_Only) return;
    
    declare_TH_parameters (TH_handler, folder_name);

    problem_name = proj_handler.get ("Problem Name");
    print_output = proj_handler.get_bool ("Print Output");
    bound_temp = TH_handler.get_double ("Dirichlet Temperature");
    bulk_temp = TH_handler.get_double ("Bulk Coolant Temperature");
    heat_trans = TH_handler.get_double ("Heat Transfer Coefficient");
    extrap = TH_handler.get_integer ("Extrapolation Method");
    
    return;
}
  
//======================================================================================================================
// Set Mesh Subroutine
template<int dim>
void HeatConduction<dim>::update_mesh (Triangulation<dim> &mesh_)
{

  // Update Remaining Mesh Information
  dof_handler.clear ();
  mesh.clear ();
  mesh.copy_triangulation(mesh_);
  dof_handler.initialize (mesh,fe);
  dof_handler.distribute_dofs (fe);
  boundary_values.clear();
  VectorTools::interpolate_boundary_values (dof_handler, 1, ConstantFunction<dim>(bound_temp), boundary_values);
    
  return;
}

//======================================================================================================================
// Initialize Problem Subroutine
template<int dim>
void HeatConduction<dim>::initialize_problem(Triangulation<dim> &mesh_)
{
    if (phys_type == Constants::Physics_Neutronics_Only) return;
    
    // Read in Mesh
    dof_handler.clear ();
    mesh.clear ();
    mesh.copy_triangulation(mesh_);
    
    // Initialize DoF Handler
    dof_handler.initialize (mesh,fe);
    dof_handler.distribute_dofs (fe);
    
    // Setup All Arrays and Finite Element Variables
    const unsigned int n_dofs = dof_handler.n_dofs();
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
    hanging_node_constraints.close ();
    
    // Reset Sparsity Pattern Variable
    sparsity_pattern.reinit (n_dofs, n_dofs, dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    hanging_node_constraints.condense (sparsity_pattern);
    sparsity_pattern.compress ();
    
    // Clear Arrays
    system_matrix.clear ();
    system_matrix.reinit (sparsity_pattern);
    system_rhs.reinit (n_dofs);
    
    // Resize and Reset Solution Vectors
    if (sol.size() == 0)
    {
        sol.reinit (n_dofs);
        sol0.reinit(n_dofs);
        sol1.reinit(n_dofs);
        solsave.reinit(n_dofs);
        sol1 = bound_temp;
        sol0 = bound_temp;
        sol = bound_temp;
        solsave = bound_temp;
    }
    
    time = 0.0;
    time0 = 0.0;
    time1 = 0.0;
    dt = 0.0;
    dt0 = 0.0;
    dtnext = 0.0;
    dtest = 0.0;
      
    // Reset Boundary Values
    boundary_values.clear();
    VectorTools::interpolate_boundary_values (dof_handler, 1, ConstantFunction<dim>(bound_temp), boundary_values);
    
    
    return;
}

//======================================================================================================================
// Assemble Right Hand Side Subroutine
template<int dim>
void HeatConduction<dim>::assemble_rhs()
{
    // Setup Local Finite Element Variables
    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    const unsigned int n_dofs = dof_handler.n_dofs ();
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values);
                             
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    
    // Setup Local Arrays
    Vector<double>          zero_src;
    Vector<double>          cell_rhs (dofs_per_cell);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    vector<double>          pCp (n_q_points);
    vector<double>          temp (n_q_points);
    vector<double>          old_sol (n_q_points);
    vector<double>          src (n_q_points);
    
    if (phys_type == Constants::Physics_TH_Only)
    {
      zero_src.reinit (n_dofs);
      zero_src = 0.01;
    }
    
    system_rhs.reinit (n_dofs);
    system_rhs = 0;
    
    // Loop through all cells
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
      
        cell_rhs = 0;
        cell->get_dof_indices (local_dof_indices);
        
        temp.resize(n_q_points, 0.0);
        pCp.resize(n_q_points, 0.0);
        
        fe_values.reinit (cell);
        fe_values.get_function_values(solsave,temp);
        fe_values.get_function_values(sol,old_sol);
        
        if (phys_type != Constants::Physics_TH_Only)
          fe_values.get_function_values(*heat_src,src);
        else
          fe_values.get_function_values(zero_src,src);

        for (unsigned int i=0; i<n_q_points; ++i)
          pCp.at(i) = material_data.get_pCp(cell->material_id(),temp.at(i));

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i) 
          {
            cell_rhs(i) += src.at(q_point)*fe_values.shape_value(i,q_point)*fe_values.JxW(q_point);
            if (prob_type == Constants::transient)
              cell_rhs(i) += pCp.at(q_point)*old_sol.at(i)*fe_values.shape_value(i,q_point)*fe_values.JxW(q_point)/dtnext;
          }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
      }
    
    return;
}


//======================================================================================================================
// Assemble System Matrix Subroutine
template<int dim>
void HeatConduction<dim>::assemble_system_matrix ()
{
    // Setup Local Finite Element Variables
    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values);
                             
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    
    // Setup Local Arrays
    FullMatrix<double>      cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>          cell_rhs (dofs_per_cell);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    vector<double>          k (n_q_points);
    vector<double>          pCp (n_q_points);
    vector<double>          temp (n_q_points);
    
    system_matrix = 0.0;
    
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        temp.resize(n_q_points, 0);
        k.resize(n_q_points, 0);
        pCp.resize(n_q_points, 0);
        fe_values.reinit (cell);
        fe_values.get_function_values(solsave,temp);

        for (unsigned int i=0; i<n_q_points; ++i)
        {
            k.at(i) = material_data.get_k(cell->material_id(),temp.at(i));
            pCp.at(i) = material_data.get_pCp(cell->material_id(),temp.at(i));
            
        }

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
          for (unsigned int i=0; i<dofs_per_cell; ++i) {
            for (unsigned int j=0; j<dofs_per_cell; ++j) {
              cell_matrix(i,j) += k.at(q_point) *fe_values.shape_grad(i,q_point) *fe_values.shape_grad(j,q_point)*fe_values.JxW(q_point);
              if (prob_type == Constants::transient)
                cell_matrix(i,j) += pCp.at(q_point)*fe_values.shape_value(i,q_point)*fe_values.shape_value(j,q_point)*fe_values.JxW(q_point)/dtnext;
            }
          }
        }

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i,j));
      }

    hanging_node_constraints.condense (system_matrix);
    
    return;
}

//================================================================
// Steady-State Solution Subroutine
template<int dim>
void HeatConduction<dim>::steady_state_solve()
{
    cout << std::setprecision (12) << std::scientific;
    print_TH_ss_header ();
    Timer timer;
    timer.start ();
    
    // Setup Initial Arrays
    assemble_system_matrix ();
    assemble_rhs ();
    
    // Begin Iteration Loop
    double error, int1, int0;
    unsigned int iteration = 1;
    int0 = 0.0;
    int1 = 0.0;
    
    do
    {
        if (temp_dep)
        {
          assemble_system_matrix ();
          assemble_rhs ();
        }
        
        solve ();
        
        int1 = get_total_energy_save();
        error = fabs(int1 - int0)/fabs(int1);
        std::cout << "   Iteration " << iteration << ":  energy = " << int1 << std::endl;
        int0 = int1;
        
        sol = solsave;
        ++iteration;
    }
    while((error > tolerance) && (iteration < Constants::MAX_ITERATION));
    
    output_results ();
    
    std::cout << std::endl;
    std::cout << "   Finished TH Steady-State Calculation,  # Iterations = " << iteration-1 << ", time = " << timer() << std::endl;
    std::cout << std::endl << std::endl;
    
    return;
}

//======================================================================================================================
// Time Step Solution Subroutine
template<int dim>
void HeatConduction<dim>::run_time_step()
{
    Timer timer;
    timer.start ();
    
    // Begin Iteration Loop
    double error, int1, int0;
    unsigned int iteration = 1;
    int0 = 0.0;
    int1 = 0.0;
    
    std::string advance_ncount   = "   Step Count:      " + Utilities::int_to_string(ncount + 1, Utilities::needed_digits(ncount+1));
    std::string current_time     = "   Current Time:    ";
    std::string advance_dt       = "   Dt:              ";
    std::string advance_time     = "   Advancing Time:  ";
    
    std::ostringstream strs1, strs2, strs3;
    strs1 << time << "   ";
    strs2 << dtnext << "   ";
    strs3 << time + dtnext << "   ";
    
    current_time += strs1.str ();
    advance_dt += strs2.str ();
    advance_time += strs3.str ();
    
    print_TH_time_step_header ();
    cout << std::endl;
    cout << advance_ncount << std::endl;
    cout << current_time << std::endl;
    cout << advance_dt << std::endl;
    cout << advance_time << std::endl << std::endl;
    
    Threads::ThreadGroup<> threads;
    
    //threads += Threads::new_thread (&HeatConduction<dim>::assemble_system_matrix, *energy_groups[group], mesh, fe, dof_handler);
    assemble_system_matrix ();
    assemble_rhs ();
    
    do 
    {
        if (temp_dep)
        {
          assemble_system_matrix ();
          assemble_rhs ();
        }
        
        int1 = get_total_energy_save();
        error = fabs(int1 - int0)/fabs(int1);
        std::cout << "   Iteration " << iteration << ":  energy = " << int1 << std::endl;
        int0 = int1;
        
        ++iteration;
    }
    while((error > tolerance) && (iteration < Constants::MAX_ITERATION));
    
    std::cout << std::endl;
    std::cout << "   Finished TH Time Step,  # Iterations = " << iteration-1 << ", time = " << timer() << std::endl;
    std::cout << std::endl << std::endl;
        
    
    return;
}

//======================================================================================================================
// Solve Subroutine
template<int dim>
void HeatConduction<dim>::solve()
{
    hanging_node_constraints.condense (system_rhs);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, solsave, system_rhs);

    SolverControl   solver_control (system_matrix.m(), 1e-12*system_rhs.l2_norm());
    SolverCG<>      cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve (system_matrix, solsave, system_rhs, preconditioner);

    hanging_node_constraints.distribute (solsave);
    
    return;
}

//======================================================================================================================
// Update Solutions and Timing Variables Subroutine
template<int dim>
void HeatConduction<dim>::update_solution()
{
    // Solutions
    sol1 = sol0;
    sol0 = sol;
    sol = solsave;
    
    // Time Variables
    time1 = time0;
    time0 = time;
    time = time + dtnext;
    dt0 = dt;
    dt = dtnext;
    
    if (prob_type == Constants::transient)
    {
      ++ncount;
      calculate_dtest ();
      output_results ();
    }
    
    return;
}

//======================================================================================================================
// output_results
template <int dim>
void HeatConduction<dim>::output_results () const 
{
    if (prob_type == Constants::transient && !print_output) return;
    
    // Local Strings
    std::string count_str = Utilities::int_to_string(ncount,Utilities::needed_digits(ncount));
    std::string fname1;
    std::string fname2;
    
    if (prob_type == Constants::steady)
    {
        fname1 = Constants::output_str + project_name + "/" + problem_name + "/TH_mesh.eps";
        fname2 = Constants::output_str + project_name + "/" + problem_name + "/Temperature_steady.gmv";
        std::ofstream mesh_output (fname1);
        
        // Ouput Grid
        GridOut grid_out;
        grid_out.write_eps (mesh, mesh_output);
    }
    else
    {
        fname2 = Constants::output_str + project_name + "/" + problem_name + "/Temperature_step-" + count_str + ".gmv";
    }
    
    // Set out data ofstream
    std::ofstream sol_output (fname2);
    
    // Output Solution Vector
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (sol, "sol");
    data_out.build_patches ();
    data_out.write_gnuplot (sol_output);
    
}

//======================================================================================================================
// Calculate Time_Dependent Temperature Subroutine
template<int dim>
Vector<double> HeatConduction<dim>::get_temperature(double t)
{
    // Local Variables
    const unsigned int n = HeatConduction<dim>::n_dofs ();
    Vector<double> temp;
    double w, mult;
    double dtt1, dtt2;
    double val1, val2;
    
    // See if input time is a solution time
    if (fabs(t - time) < Constants::SMALL)       return sol;
    else if (fabs(t - time0) < Constants::SMALL) return sol0;
    else if (fabs(t - time1) < Constants::SMALL) return sol1;
    else if (fabs(time + dtnext - t) < Constants::SMALL && advance == 1 && prob_type == Constants::transient) return solsave;
    
    // Now, do some extrapolation/interpolation instead
    if (t > time)
    {
        if (extrap == Constants::no_extrap) return sol;
        if (extrap == Constants::global) w = tot_energy/tot_energy0;
        dtt1 = dt;
        dtt2 = t - time;
    }
    else if (t > time0 && t < time)
    {
        dtt1 = dt;
        dtt2 = t - time0;
        if (extrap == Constants::no_extrap) return sol0;
        if (extrap == Constants::global)
          if (fabs(tot_energy0) > Constants::SMALL && fabs(dtt1) > Constants::SMALL)
            w = log(tot_energy/tot_energy0)/dtt1;
    }
    else if (t > time1 && t < time0)
    {
        dtt1 = dt;
        dtt2 = t - time1;
        if (extrap == Constants::no_extrap) return sol1;
        if (extrap == Constants::global) 
          if (fabs(tot_energy1) > Constants::SMALL && fabs(dtt1) > Constants::SMALL)
            w = log(tot_energy0/tot_energy1)/dtt1;
        
    }
    
    // Loop through all degrees of freedom and interpolate/extrapolate
    temp.reinit(n);
    for (unsigned int i = 0; i < n; ++i)
    {
        if ( t >= time0 )
        {
            val1 = sol0[i];
            val2 = sol[i];
            if (extrap == Constants::linear)
            {
                w = (val2 - val1)/dtt1;
                temp[i] = val2 + w*dtt2;
            }
            else if (extrap == Constants::pointwise)
            {
                if (fabs(val1) > Constants::SMALL && fabs(dtt1) > Constants::SMALL)
                  w = log(val2/val1)/dtt1;
                else
                  w = 0.0;
                temp[i] = val2*exp(w*dtt2);
            }
            else if (extrap == Constants::global) temp[i] = val2*exp(w*dtt2);
        }
        else if (t >= time1 && t < time0)
        {
            val1 = sol1[i];
            val2 = sol0[i];
            if (extrap == Constants::linear)
            {
                w = (val2 - val1)/dtt1;
                temp[i] = val2 + w*dtt2;
            }
            else if (extrap == Constants::pointwise)
            {
                if (fabs(val1) > Constants::SMALL && fabs(dtt1) > Constants::SMALL)
                  w = log(val2/val1)/dtt1;
                else
                  w = 0.0;
                temp[i] = val2*exp(w*dtt2);
            }
            else if (extrap == Constants::global) temp[i] = val2*exp(w*dtt2);
        }
    }
    
    return temp;
}

//======================================================================================================================
// Calculate Total Energy Subroutine
template<int dim>
double HeatConduction<dim>::get_total_energy(double t)
{
    // Local Variables
    double total_energy (0.0);
    const QGauss<dim>  quadrature_formula (fe.degree + 1);
    const unsigned int n_q_points    = quadrature_formula.size();
    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_JxW_values);
    
    Vector<double>  temp_vec;
    vector<double>  pCp (n_q_points);
    vector<double>  temp (n_q_points);
    
    temp_vec = get_temperature(t);
    
    // Begin Execution
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        temp.resize(n_q_points, 0);
        pCp.resize(n_q_points, 0);
        fe_values.reinit (cell);
        fe_values.get_function_values(temp_vec,temp);

        for (unsigned int i=0; i<n_q_points; ++i)
            pCp.at(i) = material_data.get_pCp(cell->material_id(),temp.at(i));

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
                total_energy += pCp.at(q_point)*temp.at(q_point)*fe_values.shape_value(i,q_point)*fe_values.JxW(q_point);
    }

    return total_energy;
}

//======================================================================================================================
// Calculate Total Energy Subroutine
template<int dim>
double HeatConduction<dim>::get_total_energy_save()
{
    // Local Variables
    double total_energy (0.0);
    const QGauss<dim>  quadrature_formula (fe.degree + 1);
    const unsigned int n_q_points    = quadrature_formula.size();
    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_JxW_values);
    
    vector<double>  pCp (n_q_points);
    vector<double>  temp (n_q_points);
    
    // Begin Execution
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        temp.resize(n_q_points, 0);
        pCp.resize(n_q_points, 0);
        fe_values.reinit (cell);
        fe_values.get_function_values(solsave,temp);

        for (unsigned int i=0; i<n_q_points; ++i)
            pCp.at(i) = material_data.get_pCp(cell->material_id(),temp.at(i));

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
                total_energy += pCp.at(q_point)*temp.at(q_point)*fe_values.shape_value(i,q_point)*fe_values.JxW(q_point);
    }

    return total_energy;
}

//======================================================================================================================
// calculate dtest 
template <int dim>
void HeatConduction<dim>::calculate_dtest () 
{
  // Local Variables
  double dt_temp, sderest;
  
  switch(dt_estimator) {
    case Constants::absolute_trunc_err:
      sderest = calculate_maximum_absolute_secder_stencil ();
      if (fabs(sderest) > Constants::SMALL)
        dt_temp = 2*adaptive_err_tol/fabs(sderest);
      else
        dt_temp = 1.0e10;
      break;
    case Constants::relative_trunc_err:
      sderest = calculate_maximum_relative_secder_stencil ();
      if (fabs(sderest) > Constants::SMALL)
        dt_temp = sqrt(2*adaptive_err_tol/fabs(sderest));
      else
        dt_temp = 1.0e10;
      break;
    case Constants::max_local_period:
      dt_temp = calculate_minimum_local_period ();
      break;
  }
  
  dtest = dt_temp;
  
  return;
}

//======================================================================================================================
// calculate maximum second derivative 
template <int dim>
double HeatConduction<dim>::calculate_maximum_absolute_secder_stencil () 
{
    if (dt < Constants::SMALL || dt0 < Constants::SMALL) return 0.0;
    
    double sdermax, sderest;
    sdermax = 0.0; sderest = 0.0;
    
    for (unsigned int i = 0; i < sol.size(); ++i)
    {
        sderest = (dt0*sol[i] - (dt+dt0)*sol0[i] + dt*sol1[i])/(dt*dt*dt0);
        if (fabs(sderest) > sdermax)
          sdermax = fabs(sderest);
    }
    
    return sderest;
}

//======================================================================================================================
// calculate maximum relative second derivative 
template <int dim>
double HeatConduction<dim>::calculate_maximum_relative_secder_stencil () 
{
    if (dt < Constants::SMALL || dt0 < Constants::SMALL) return 0.0;
    
    double sdermax (0.0);
    double sderest (0.0);
    
    for (unsigned int i = 0; i < sol.size(); ++i)
    {
        sderest = (dt0*sol[i] - (dt+dt0)*sol0[i] + dt*sol1[i])/(dt*dt*dt0);
        if (fabs(sderest)/fabs(sol[i]) > sdermax && sol[i] != 0.0)
          sdermax = fabs(sderest/sol[i]);
    }
    
    return sderest;
}

//======================================================================================================================
// calculate minimum local period
template <int dim>
double HeatConduction<dim>::calculate_minimum_local_period () 
{
    if (dt < Constants::SMALL) return 0.0;
    
    double dt_est (0.0);
    double max_temp (0.0);
    
    for (unsigned int i = 0; i < sol.size(); ++i)
    {
        if (fabs(sol[i] + sol0[i]) < Constants::SMALL || fabs(dt) < Constants::SMALL)
          continue;
        else
        {
          max_temp = (sol[i] - sol0[i])/(sol[i] + sol0[i])*(2/dt);
          if (fabs(max_temp) > dt_est)
            dt_est = fabs(max_temp);
        }
    }
    
    if (fabs(dt_est) > Constants::SMALL)
      dt_est = 1.0 / dt_est;
    else
      dt_est = 1.0e10;
    
    return dt_est;
}

//======================================================================================================================
// Terminate Problem Subroutine
template<int dim>
void HeatConduction<dim>::terminate_problem ()
{
  // Local Variables
  
  return;
}

//======================================================================================================================


#endif
