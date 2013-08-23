#ifndef energy_group_h_
#define energy_group_h_

//Global Headers
#include "../Global/GlobalHeaders.h"
#include "../Global/Deal2Headers.h"
#include "../Global/Constants.h"

#include "MaterialData.h"

// Thermal Hydraulics Header
#include "../TH/HeatConduction.h"

//======================================================================================================================
template <int dim>
class EnergyGroup {

    public:

      // Constructor
      EnergyGroup (const unsigned int group,
                   const MaterialData &material_data,
                   Triangulation<dim> &mesh,
                   FiniteElement<dim> &fe);
      
      // Destructor
      ~EnergyGroup ();
      
      // Mesh Routines
      void setup_linear_system ();
      void update_mesh (const Triangulation<dim> &mesh_, const FiniteElement<dim> &fe_, const DoFHandler<dim> &dof_);
      unsigned int n_active_cells () const;
      unsigned int n_dofs () const;
      
      // Matrix and Vector Assembly Routines
      void assemble_system_matrix ();
      void assemble_adjoint_system_matrix ();
      void assemble_ingroup_rhs (const Function<dim> &extraneous_source);
      void assemble_adjoint_ingroup_rhs (const Function<dim> &extraneous_source);
      void assemble_cross_group_rhs (EnergyGroup<dim> &g_prime);
      void assemble_adjoint_cross_group_rhs (EnergyGroup<dim> &g_prime);
      void assemble_precursor_rhs (const vector<Vector<double>> &prec_in, const DoFHandler<dim> &dof_in);
      void assemble_previous_flux_rhs ();

      double get_fission_source ();
      double get_energy_source (const double t, const unsigned int flag);
      double get_neutron_population (const double t, const unsigned int flag);
      Vector<double> get_flux (double t) const;
      vector<double> get_local_frequencies (double t) const;
      
      // Solution Routines
      void solve ();
      void update_solution ();
      void update_powers ();
      void output_results () const;
      
      // Manipulators
      inline void set_folder_name (const std::string &inp_str) {folder_name = inp_str;};
      inline void set_project_name (const std::string &inp_str) {project_name = inp_str;};
      inline void set_problem_name (const std::string &inp_str) {problem_name = inp_str;};
      void set_TH_pointer (HeatConduction<dim> *th) {TH = th;}
      
      inline void set_keff (double k) {keff=k;};
      inline void set_dtnext (double t) {dtnext=t;};
      inline void set_prec_time (double t) {prec_time=t;};
      inline void set_adaptive_err_tol (double t) {adaptive_err_tol=t;};
      
      inline void set_prob_type (unsigned int p) {prob_type = p;};
      inline void set_phys_type (unsigned int p) {phys_type = p;};
      inline void set_extrap (unsigned int p) {extrap = p;};
      inline void set_advance (unsigned int p) {advance = p;};
      inline void set_time_step_selection (unsigned int p) {time_step_selection=p;}
      inline void set_dt_estimator (unsigned int p) {dt_estimator = p;};
      inline void set_use_adjoint_weighting (bool p) {use_adjoint_weighting = p;};
      inline void set_print_output (bool p) {print_output = p;};
      
      inline void mult_solsave (double val) {solsave *= val;};
      inline void mult_sol (double val) {sol *= val;};
      inline void mult_sol0 (double val) {sol0 *= val;};
      inline void mult_sol1 (double val) {sol1 *= val;};
      
      void set_adjoint (Vector<double> val) {adjoint=val;};
      void set_solsave (Vector<double> val) {solsave=val;};
      void set_sol (Vector<double> val) {sol=val;};
      void set_sol0 (Vector<double> val) {sol0=val;};
      void set_sol1 (Vector<double> val) {sol1=val;};
            
      // Accessors
      Vector<double> get_solsave () {return solsave;};
      Vector<double> get_sol () {return sol;};
      Vector<double> get_sol0 () {return sol0;};
      Vector<double> get_sol1 () {return sol1;};
      Vector<double> get_adjoint1 () {return adjoint;};
      
      inline double get_time () {return time;};
      inline double get_time0 () {return time0;};
      inline double get_time1 () {return time1;};
      inline double get_dt () {return dt;};
      inline double get_dt0 () {return dt0;};
      inline double get_dtnext () {return dtnext;};
      inline double get_dtest () {return dtest;};
      inline double get_next_time () {return time+dtnext;};
      
      inline double get_total_power () {return total_power;};
      inline double get_total_power0 () {return total_power0;};
      inline double get_total_power1 () {return total_power1;};
      
      inline unsigned int get_prob_type () {return prob_type;};
      inline unsigned int get_phys_type () {return phys_type;};
      inline unsigned int get_extrap () {return extrap;};
      inline unsigned int get_advance () {return advance;};
      inline unsigned int get_ncount () {return ncount;};

    public:

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
      double prec_time;
      
    private:
      
      // control variables
      unsigned int  prob_type;
      unsigned int  phys_type;
      unsigned int  extrap;
      unsigned int  advance;
      unsigned int  time_step_selection;
      unsigned int  dt_estimator;
      unsigned int  ncount;
      
      bool          use_adjoint_weighting;
      bool          print_output;
      
      const unsigned int    group;
      double                keff;
      
      double total_power;
      double total_power0;
      double total_power1;
      double adaptive_err_tol;
      
      std::string           folder_name;
      std::string           project_name;
      std::string           problem_name;
      
      const MaterialData    &material_data;

      Triangulation<dim>    mesh;
      FiniteElement<dim>    &fe;
      DoFHandler<dim>       dof_handler;

      SparsityPattern       sparsity_pattern;
      SparseMatrix<double>  system_matrix;

      Vector<double>        system_rhs;

      std::map<unsigned int,double> boundary_values;
      ConstraintMatrix              hanging_node_constraints;
      
    private: // functions
      void      calculate_dtest ();
      double    calculate_maximum_absolute_secder_stencil ();
      double    calculate_maximum_relative_secder_stencil ();
      double    calculate_minimum_local_period ();
      double    get_fiss_frac () const;
      Vector<double> flux_extrap (double t) const;
    
    private: // TH Classes
      HeatConduction<dim>*      TH;
};

//======================================================================================================================
// EnergyGroup Constructor
template <int dim>
EnergyGroup<dim>::EnergyGroup (const unsigned int group,
                               const MaterialData &material_data,
                               Triangulation<dim> &mesh_,
                               FiniteElement<dim> &fe_)
                  :
                  group (group),
                  material_data (material_data),
                  fe (fe_),
                  dof_handler (mesh_)
                  {
                    mesh.copy_triangulation (mesh_);
                    dof_handler.distribute_dofs (fe_);
                    ncount = 0;
                  }

//======================================================================================================================
// EnergyGroup Destructor
template <int dim>
EnergyGroup<dim>::~EnergyGroup () 
{
dof_handler.clear ();
mesh.clear ();
TH = NULL;
}

//======================================================================================================================
// n_active_cells
template <int dim>
unsigned int EnergyGroup<dim>::n_active_cells () const
{
    return mesh.n_active_cells ();
}

//======================================================================================================================
// n_dofs
template <int dim>
unsigned int EnergyGroup<dim>::n_dofs () const
{
    return dof_handler.n_dofs ();
}


//======================================================================================================================
// setup_linear_system
template <int dim>
void EnergyGroup<dim>::setup_linear_system ()
{
    const unsigned int n_dofs = EnergyGroup<dim>::n_dofs ();

    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
    hanging_node_constraints.close ();

    system_matrix.clear ();

    sparsity_pattern.reinit (n_dofs, n_dofs, dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    hanging_node_constraints.condense (sparsity_pattern);
    sparsity_pattern.compress ();

    system_matrix.reinit (sparsity_pattern);

    system_rhs.reinit (n_dofs);

    if (sol.size() == 0)
      {
        sol.reinit (n_dofs);
        sol0.reinit(n_dofs);
        sol1.reinit(n_dofs);
        solsave.reinit(n_dofs);
        sol1 = 0.0;
        sol0 = 0.0;
        sol = 1.0;
        solsave = 1.0;
      }
    
    time = 0.0;
    time0 = 0.0;
    time1 = 0.0;
    dt = 0.0;
    dt0 = 0.0;
    dtnext = 0.0;
    dtest = 0.0;

    boundary_values.clear();
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(), boundary_values);
}

//======================================================================================================================
// update_mesh
template <int dim>
void EnergyGroup<dim>::update_mesh (const Triangulation<dim> &mesh_, const FiniteElement<dim> &fe_, const DoFHandler<dim> &dof_)
{
    // Local Variables
    Vector<double>      temp_sol (dof_.n_dofs ());
    
    // Interpolate Solution Vectors
    //temp_sol.swap (sol1);
    //VectorTools::interpolate_to_different_mesh (dof_handler, sol1, dof_, temp_sol);
    //temp_sol.swap (sol0);
    //VectorTools::interpolate_to_different_mesh (dof_handler, sol0, dof_, temp_sol);
    //temp_sol.swap (sol);
    //VectorTools::interpolate_to_different_mesh (dof_handler, sol, dof_, temp_sol);
    //temp_sol.swap (solsave);
    //VectorTools::interpolate_to_different_mesh (dof_handler, solsave, dof_, temp_sol);
    
    // Reset Mesh Data
    dof_handler.clear ();
    mesh.clear ();
    mesh.copy_triangulation (mesh_);
    dof_handler.initialize (mesh,fe);
    dof_handler.distribute_dofs (fe);
    
    // Reset Boundary Values
    boundary_values.clear();
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(), boundary_values);
    
    return;
}

//======================================================================================================================
// assemble_system_matrix
template <int dim>
void EnergyGroup<dim>::assemble_system_matrix () 
{
    
    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values);
    
    const unsigned int  dofs_per_cell = fe.dofs_per_cell;
    const unsigned int  n_q_points    = quadrature_formula.size();
    const double        vel = material_data.get_velocity(group);
    
    FullMatrix<double>      cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>          cell_rhs (dofs_per_cell);
    Vector<double>          temperature;
    vector<double>          diff_XS (n_q_points, 0.0);
    vector<double>          rem_XS (n_q_points, 0.0);
    vector<double>          temp_values (n_q_points, 1.0);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    
    if (phys_type == Constants::Physics_Both) 
      temperature.reinit (TH->get_temperature (get_next_time ()));
    
    system_matrix = 0.0;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        fe_values.reinit (cell);
        
        if (phys_type == Constants::Physics_Both)
          fe_values.get_function_values (temperature, temp_values);
        
        for (unsigned int q_point=0; q_point < n_q_points; ++q_point)
        {
          if (phys_type == Constants::Physics_Both)
          {
            diff_XS.at(q_point) = material_data.get_diffusion_coefficient (group, cell->material_id(), temp_values.at(q_point));
            rem_XS.at(q_point) = material_data.get_removal_XS (group, cell->material_id(), temp_values.at(q_point));
          }
          else
          {
            diff_XS.at(q_point) = material_data.get_diffusion_coefficient (group, cell->material_id(), temp_values.at(q_point));
            rem_XS.at(q_point) = material_data.get_removal_XS (group, cell->material_id(), temp_values.at(q_point));
          }
        }
        
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
          for (unsigned int i=0; i<dofs_per_cell; ++i) {
            for (unsigned int j=0; j<dofs_per_cell; ++j) {
              if (prob_type == Constants::transient)
                cell_matrix(i,j) += fe_values.shape_value(i,q_point) * fe_values.shape_value(j,q_point) * fe_values.JxW(q_point) / (vel*dtnext);
              cell_matrix(i,j) += diff_XS.at(q_point) * fe_values.shape_grad(i,q_point) * fe_values.shape_grad(j,q_point) * fe_values.JxW(q_point);
              cell_matrix(i,j) += rem_XS.at(q_point) * fe_values.shape_value(i,q_point) * fe_values.shape_value(j,q_point) * fe_values.JxW(q_point);
            }
          }
        }

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i,j));
      }

    hanging_node_constraints.condense (system_matrix);
}

//======================================================================================================================
// assemble_ingroup_rhs
template <int dim>
void EnergyGroup<dim>::assemble_ingroup_rhs (const Function<dim> &extraneous_source) 
{

    system_rhs.reinit (dof_handler.n_dofs());
    system_rhs = 0.0;
    if (prob_type == Constants::transient) assemble_previous_flux_rhs ();

    const QGauss<dim>  quadrature_formula (fe.degree + 1);

    // Local Constant Variables
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const double vel = material_data.get_velocity(group);
    const double fiss_frac = get_fiss_frac ();
    
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);
    
    // Local Cell-wise Variables
    Vector<double>          cell_rhs (dofs_per_cell);
    Vector<double>          temperature;
    vector<double>          extraneous_source_values (n_q_points, 0.0);
    vector<double>          sol_values (n_q_points, 0.0);
    vector<double>          sol_old_values (n_q_points, 0.0);
    vector<double>          temp_values (n_q_points, 1.0);
    vector<double>          fiss_XS (n_q_points, 0.0);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    
    if (phys_type == Constants::Physics_Both)
      temperature.reinit (TH->get_temperature (get_next_time ()));
    
    // Loop Through all Cells and Build RHS
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        cell_rhs = 0;
        fe_values.reinit (cell);
        
        if (phys_type == Constants::Physics_Both)
          fe_values.get_function_values (temperature, temp_values);
        
        for (unsigned int q_point=0; q_point < n_q_points; ++q_point)
          if (phys_type == Constants::Physics_Both)
            fiss_XS.at(q_point) = material_data.get_fission_dist_XS (group, group, cell->material_id(), temp_values.at(q_point));
          else
            fiss_XS.at(q_point) = material_data.get_fission_dist_XS (group, group, cell->material_id(), temp_values.at(q_point));

        cell->get_dof_indices (local_dof_indices);
        extraneous_source.value_list (fe_values.get_quadrature_points(), extraneous_source_values);
        if (prob_type == Constants::steady)
          fe_values.get_function_values (sol, sol_values);
        else
          fe_values.get_function_values (solsave, sol_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            // In-Group Fission Terms
            if (prob_type == Constants::transient)
              cell_rhs(i) += fiss_frac * fiss_XS.at(q_point)/keff * sol_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
            else
              cell_rhs(i) += fiss_XS.at(q_point) * sol_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
            // Extraneous Source
            cell_rhs(i) += extraneous_source_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
          }
        }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
}

//======================================================================================================================
// assemble_cross_group_rhs
template <int dim>
void EnergyGroup<dim>::assemble_cross_group_rhs (EnergyGroup<dim> &g_prime)
{
    if (group == g_prime.group) return;
      
    const QGauss<dim>  quadrature_formula (fe.degree + 1);

    // Local Constant Variables
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int gp = g_prime.group;
    const double vel = material_data.get_velocity(group);
    const double fiss_frac = get_fiss_frac ();
    double t = time + dtnext;
    
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);
    
    // Local Cell-wise Variables
    Vector<double>          cell_rhs (dofs_per_cell);
    Vector<double>          flux;
    Vector<double>          temperature;
    vector<double>          extraneous_source_values (n_q_points, 0.0);
    vector<double>          gp_values (n_q_points, 0.0);
    vector<double>          gp_old_values (n_q_points, 0.0);
    vector<double>          temp_values (n_q_points, 1.0);
    vector<double>          fiss_XS (n_q_points, 0.0);
    vector<double>          scatt_XS (n_q_points, 0.0);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    
    if (phys_type == Constants::Physics_Both) 
      temperature.reinit (TH->get_temperature (get_next_time ()));
    
    if (prob_type == Constants::transient && g_prime.get_advance() == 0)
      flux.reinit (flux = g_prime.get_flux (t));
    
    // Loop Through all Cells and Build RHS
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        cell_rhs = 0;
        fe_values.reinit (cell);
        
        if (phys_type == Constants::Physics_Both)
          fe_values.get_function_values (temperature, temp_values);
        
        for (unsigned int q_point=0; q_point < n_q_points; ++q_point)
        {
          fiss_XS.at(q_point) = material_data.get_fission_dist_XS (group, gp, cell->material_id(), temp_values.at(q_point));
          scatt_XS.at(q_point) = material_data.get_scattering_XS (gp, group, cell->material_id(), temp_values.at(q_point));
        }

        cell->get_dof_indices (local_dof_indices);
        if (prob_type == Constants::steady)
          if (gp < group)
            fe_values.get_function_values (g_prime.solsave, gp_values);
          else
            fe_values.get_function_values (g_prime.sol, gp_values);
        else
        {
          if (g_prime.get_advance() == 1)
            fe_values.get_function_values (g_prime.solsave, gp_values);
          else
            fe_values.get_function_values (flux, gp_values);
        }

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            if (prob_type == Constants::transient)
            {
              cell_rhs(i) += fiss_frac * fiss_XS.at(q_point)/keff * gp_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
              cell_rhs(i) += scatt_XS.at(q_point) * gp_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
            }
            else
            {
              cell_rhs(i) += fiss_frac * fiss_XS.at(q_point) * gp_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
              cell_rhs(i) += scatt_XS.at(q_point) * gp_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
            }
          }
        }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
    
    return;
}

//======================================================================================================================
// assemble_precursor_rhs
template <int dim>
void EnergyGroup<dim>::assemble_precursor_rhs (const vector<Vector<double>> &prec_in, const DoFHandler<dim> &dof_in) 
{
    if (prob_type == Constants::steady) return;
    const QGauss<dim>  quadrature_formula (fe.degree + 1);

    // Local Constant Variables
    const unsigned int  dofs_per_cell = fe.dofs_per_cell;
    const unsigned int  n_q_points = quadrature_formula.size();
    const unsigned int  n_prec = material_data.get_n_prec ();
    const double        fiss_frac = get_fiss_frac ();
    const double        dt_prec = get_next_time () - prec_time;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);
    
    // Local Cell-wise Variables
    Vector<double>          cell_rhs (dofs_per_cell);
    vector<double>          extraneous_source_values (n_q_points);
    vector<double>          sol_values (n_q_points);
    vector<double>          prec_frac (n_prec);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    vector< vector<double>> prec_values (n_prec, vector<double>(n_q_points,0.0));
    
    // Loop Through all Cells and Build RHS
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        fe_values.reinit (cell);
        cell_rhs = 0;
        
        // Compute Precursor Fractions if in Transient Mode
        if (prob_type == Constants::transient)
          for (unsigned int i = 0; i < n_prec; ++i)
            prec_frac.at(i) = material_data.get_chi_d(group,i,cell->material_id())*material_data.get_lambda(i)/(1 + material_data.get_lambda(i)*dt_prec);
                
        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i = 0; i < n_prec; ++i)
          fe_values.get_function_values (prec_in[i], prec_values[i]);
        
        // Loop through precursors and for local rhs
        for (unsigned int q_prec=0; q_prec<n_prec; ++q_prec)
          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              cell_rhs(i) += prec_frac[q_prec] * prec_values[q_prec][q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
        
        // Add local rhs vector to global rhs vector
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
}

//======================================================================================================================
// assemble_previous_flux_rhs
template <int dim>
void EnergyGroup<dim>::assemble_previous_flux_rhs () 
{
    if (prob_type == Constants::steady) return;
    const QGauss<dim>  quadrature_formula (fe.degree + 1);

    // Local Constant Variables
    const unsigned int  n_dofs = dof_handler.n_dofs ();
    const unsigned int  dofs_per_cell = fe.dofs_per_cell;
    const unsigned int  n_q_points = quadrature_formula.size();
    const double        vel = material_data.get_velocity(group);

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);
    
    // Local Cell-wise Variables
    Vector<double>          cell_rhs (dofs_per_cell);
    vector<double>          sol_values (n_q_points);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    
    // Loop Through all Cells and Build RHS
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        cell_rhs = 0.0;
        cell->get_dof_indices (local_dof_indices);
        fe_values.reinit (cell);
        fe_values.get_function_values (sol, sol_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            cell_rhs(i) += sol_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point) / (vel*dtnext);
        
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
    
    
    return;
}

//======================================================================================================================
// assemble_adjoint_system_matrix
template <int dim>
void EnergyGroup<dim>::assemble_adjoint_system_matrix () 
{
    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

    const unsigned int  dofs_per_cell = fe.dofs_per_cell;
    const unsigned int  n_q_points    = quadrature_formula.size();

    FullMatrix<double>      cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>          cell_rhs (dofs_per_cell);
    Vector<double>          temperature;
    vector<double>          diff_XS (n_q_points, 0.0);
    vector<double>          rem_XS (n_q_points, 0.0);
    vector<double>          temp_values (n_q_points, 1.0);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    
    if (phys_type == Constants::Physics_Both)
      temperature.reinit (TH->get_temperature (get_next_time ()));
    
    system_matrix = 0.0;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        fe_values.reinit (cell);
        
        if (phys_type == Constants::Physics_Both)
          fe_values.get_function_values (temperature, temp_values);
        
        for (unsigned int q_point=0; q_point < n_q_points; ++q_point)
        {
          diff_XS.at(q_point) = material_data.get_diffusion_coefficient (group, cell->material_id(), temp_values.at(q_point));
          rem_XS.at(q_point) = material_data.get_removal_XS (group, cell->material_id(), temp_values.at(q_point));
        }
        
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
          for (unsigned int i=0; i<dofs_per_cell; ++i) {
            for (unsigned int j=0; j<dofs_per_cell; ++j) {
              cell_matrix(i,j) += diff_XS.at(q_point) * fe_values.shape_grad(i,q_point) * fe_values.shape_grad(j,q_point) * fe_values.JxW(q_point);
              cell_matrix(i,j) += rem_XS.at(q_point) * fe_values.shape_value(i,q_point) * fe_values.shape_value(j,q_point) * fe_values.JxW(q_point);
            }
          }
        }

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i,j));
      }

    hanging_node_constraints.condense (system_matrix);
}

//======================================================================================================================
// assemble_adjoint_ingroup_rhs
template <int dim>
void EnergyGroup<dim>::assemble_adjoint_ingroup_rhs (const Function<dim> &extraneous_source) 
{

    system_rhs.reinit (dof_handler.n_dofs());
    system_rhs = 0.0;

    const QGauss<dim>  quadrature_formula (fe.degree + 1);

    // Local Constant Variables
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);
    
    // Local Cell-wise Variables
    Vector<double>          cell_rhs (dofs_per_cell);
    Vector<double>          temperature;
    vector<double>          extraneous_source_values (n_q_points, 0.0);
    vector<double>          sol_values (n_q_points, 0.0);
    vector<double>          temp_values (n_q_points, 1.0);
    vector<double>          fiss_XS (n_q_points, 0.0);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    
    if (phys_type == Constants::Physics_Both)
      temperature.reinit (TH->get_temperature (get_next_time ()));
    
    // Loop Through all Cells and Build RHS
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        cell_rhs = 0;
        fe_values.reinit (cell);
        
        if (phys_type == Constants::Physics_Both)
          fe_values.get_function_values (temperature, temp_values);
        
        for (unsigned int q_point=0; q_point < n_q_points; ++q_point)
          fiss_XS.at(q_point) = material_data.get_fission_dist_XS (group, group, cell->material_id(), temp_values.at(q_point));

        cell->get_dof_indices (local_dof_indices);
        extraneous_source.value_list (fe_values.get_quadrature_points(), extraneous_source_values);
        fe_values.get_function_values (solsave, sol_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            cell_rhs(i) += fiss_XS.at(q_point) / keff * sol_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
            cell_rhs(i) += extraneous_source_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
          }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
}

//======================================================================================================================
// assemble_adjoint_cross_group_rhs
template <int dim>
void EnergyGroup<dim>::assemble_adjoint_cross_group_rhs (EnergyGroup<dim> &g_prime)
{
    if (group == g_prime.group) return;
      
    const QGauss<dim>  quadrature_formula (fe.degree + 1);

    // Local Constant Variables
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int gp = g_prime.group;
    const double vel = material_data.get_velocity(group);
    double t = time + dtnext;
    
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);
    
    // Local Cell-wise Variables
    Vector<double>          cell_rhs (dofs_per_cell);
    Vector<double>          temperature;
    vector<double>          extraneous_source_values (n_q_points, 0.0);
    vector<double>          gp_values (n_q_points, 0.0);
    vector<double>          temp_values (n_q_points, 1.0);
    vector<double>          fiss_XS (n_q_points, 0.0);
    vector<double>          scatt_XS (n_q_points, 0.0);
    vector<unsigned int>    local_dof_indices (dofs_per_cell);
    
    if (phys_type == Constants::Physics_Both)
      temperature.reinit (TH->get_temperature (get_next_time ()));
    
    // Loop Through all Cells and Build RHS
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        cell_rhs = 0;
        fe_values.reinit (cell);
        
        if (phys_type == Constants::Physics_Both)
          fe_values.get_function_values (temperature, temp_values);
        
        for (unsigned int q_point=0; q_point < n_q_points; ++q_point)
        {
          fiss_XS.at(q_point) = material_data.get_fission_XS (group, cell->material_id(), temp_values.at(q_point)) * material_data.get_chi (gp, cell->material_id());
          scatt_XS.at(q_point) = material_data.get_scattering_XS (group, gp, cell->material_id(), temp_values.at(q_point));
        }

        cell->get_dof_indices (local_dof_indices);
        fe_values.get_function_values (g_prime.solsave, gp_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            cell_rhs(i) += fiss_XS.at(q_point) / keff * gp_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
            cell_rhs(i) += scatt_XS.at(q_point) * gp_values[q_point] * fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);
          }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
    
    return;
}

//======================================================================================================================
// get_fission_source
template <int dim>
double EnergyGroup<dim>::get_fission_source () 
{
    const QGauss<dim>  quadrature_formula (fe.degree + 1);
    const unsigned int n_q_points    = quadrature_formula.size();

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_JxW_values);

    vector<double>  sol_values (n_q_points, 0.0);
    vector<double>  temp_values (n_q_points, 1.0);
    vector<double>  fiss_XS (n_q_points, 0.0);
    Vector<double>  temperature;

    double fission_source = 0;
    
    if (phys_type == Constants::Physics_Both)
      temperature.reinit (TH->get_temperature (get_next_time ()));

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        
        if (phys_type == Constants::Physics_Both)
          fe_values.get_function_values (temperature, temp_values);

        for (unsigned int q_point=0; q_point < n_q_points; ++q_point)
          fiss_XS.at(q_point) = material_data.get_fission_dist_XS (group, group, cell->material_id(), temp_values.at(q_point));

        fe_values.get_function_values (solsave, sol_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          fission_source += (fiss_XS.at(q_point) * sol_values[q_point] * fe_values.JxW(q_point));
      }

    return fission_source;
}

//======================================================================================================================
// get_energy_source
template <int dim>
double EnergyGroup<dim>::get_energy_source (const double t, const unsigned int flag) 
{
    const QGauss<dim>  quadrature_formula (fe.degree + 1);
    const unsigned int n_q_points    = quadrature_formula.size();
    const double energy_conversion (Constants::ENERGY_FROM_FISS*Constants::MEV_TO_J);

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_JxW_values);
    
    double          energy_source = 0;
    vector<double>  sol_values (n_q_points, 0.0);
    vector<double>  fiss_XS (n_q_points, 0.0);
    vector<double>  temp_values (n_q_points, 1.0);
    Vector<double>  temperature;
    Vector<double>  flux;
    
    if (prob_type == Constants::transient && flag != 0 && advance == 0)
      flux = get_flux(t);
    
    if (phys_type == Constants::Physics_Both)
      temperature.reinit (TH->get_temperature (get_next_time ()));
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);

        for (unsigned int q_point=0; q_point < n_q_points; ++q_point)
        {
          fiss_XS.at(q_point) = material_data.get_fission_dist_XS (group, group, cell->material_id(), temp_values.at(q_point));
          fiss_XS.at(q_point) /= material_data.get_nu(group, cell->material_id());
        }

        if (prob_type == Constants::transient && flag != 0 && advance == 0)
          fe_values.get_function_values (flux, sol_values);
        else
          fe_values.get_function_values (solsave, sol_values);
        
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          energy_source += (fiss_XS.at(q_point) * energy_conversion * sol_values[q_point] * fe_values.JxW(q_point));
      }

    return energy_source;
}

//======================================================================================================================
// get_neutron_population
template <int dim>
double EnergyGroup<dim>::get_neutron_population (const double t, const unsigned int flag) 
{
    const QGauss<dim>  quadrature_formula (fe.degree + 1);
    const unsigned int n_q_points    = quadrature_formula.size();
    const double       velocity = material_data.get_velocity (group);

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_JxW_values);

    vector<double>  sol_values (n_q_points);
    vector<double>  wgt_values (n_q_points);
    Vector<double>  flux;
    
    if (prob_type == Constants::transient && flag != 0 && advance == 0)
      flux = get_flux(t);

    double neutron_pop = 0;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);


        if (prob_type == Constants::transient && flag != 0 && advance == 0)
          fe_values.get_function_values (flux, sol_values);
        else
          fe_values.get_function_values (solsave, sol_values);
        
        if (use_adjoint_weighting)
          fe_values.get_function_values (adjoint, wgt_values);
        else
          wgt_values.resize(n_q_points, 1.0);
        
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          neutron_pop += (sol_values[q_point] * wgt_values[q_point] * fe_values.JxW(q_point) / velocity);
      }

    return neutron_pop;
}

//======================================================================================================================
// solve
template <int dim>
void EnergyGroup<dim>::solve () 
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
void EnergyGroup<dim>::update_solution()
{
    
    if (prob_type == Constants::transient)
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
    
      ++ncount;
      update_powers ();
      total_power = get_energy_source (time, 1);
      calculate_dtest ();
      output_results ();
    }
    else
    {
      sol = solsave;
      total_power = get_energy_source (0.0, 0);
    }
    
    return;
}

//======================================================================================================================
// Update Reactor Powers Subroutine
template<int dim>
void EnergyGroup<dim>::update_powers()
{
  total_power1 = total_power0;
  total_power0 = total_power;
  return;
}

//======================================================================================================================
// output_results
template <int dim>
void EnergyGroup<dim>::output_results () const 
{
    if (prob_type == Constants::transient && !print_output) return;
    
    // Local Strings
    std::string group_str = Utilities::int_to_string(group+1,Utilities::needed_digits(group+1));
    std::string fname1;
    std::string fname2;
    
    if (prob_type == Constants::steady)
    {
        fname1 = Constants::output_str + project_name + "/" + problem_name + "/" + std::string("kinetics_mesh_group-") + group_str + ".eps";
        fname2 = Constants::output_str + project_name + "/" + problem_name + "/" + std::string("kinetics_group-") + group_str + "_steady.gmv";
        std::ofstream mesh_output (fname1);
        
        // Ouput Grid
        GridOut grid_out;
        grid_out.write_eps (mesh, mesh_output);
    }
    else
    {
        fname2 = Constants::output_str + project_name + "/" + problem_name + "/" + std::string("kinetics_group-") 
                + group_str + "_step-" + Utilities::int_to_string(ncount,Utilities::needed_digits(ncount)) + ".gmv";
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
// get_flux
template <int dim>
Vector<double> EnergyGroup<dim>::get_flux (double t) const 
{
  
  if (prob_type == Constants::steady) {return solsave;}
  else if (advance == 1) {return solsave;}
  
  Vector<double> flux;
  
  if (fabs(t-time) < Constants::SMALL) {return sol;}
  else if (fabs(t-time0) < Constants::SMALL) {return sol0;}
  else if (fabs(t-time1) < Constants::SMALL) {return sol1;}
  else
  {
    if (t > time && fabs(dt) < Constants::SMALL) {return sol;}
    else if (t > time0 && t < time && fabs(dt) < Constants::SMALL) {return sol;}
    else if (t > time1 && t < time0 && fabs(dt) < Constants::SMALL) {return sol;}
    else if (t > time1 && t < time0 && fabs(dt0) < Constants::SMALL) {return sol0;}
    else
    {
      flux = EnergyGroup<dim>::flux_extrap (t);
      return flux;
    }
  }

}

//======================================================================================================================
// flux_extrap
template <int dim>
Vector<double> EnergyGroup<dim>::flux_extrap (double t) const 
{

  double dt, dtt, w;
  Vector<double> flux(sol.size());
  Vector<double> flux0(sol.size());
  
  if (t >= time ) 
  {
    dt = t - time;
    dtt = dt;
    flux = sol;
    flux0 = sol0;
  }
  else if (t > time0 && t < time)
  {
    dt = t - time0;
    dtt = dt;
    flux = sol;
    flux0 = sol0;
  }
  else if (t > time1 && t < time0)
  {
    dt = t - time1;
    dtt = dt0;
    flux = sol0;
    flux0 = sol1;
  }
  
  if (extrap == Constants::global)
  {
    if (t >= time0 )
    {
      if (fabs(total_power0) > Constants::SMALL && fabs(dt) > Constants::SMALL)
        w = log(total_power/total_power0)/dt;
      else
        w = 0.0;
    }
    else
    {
      if (fabs(total_power1) > Constants::SMALL && fabs(dt0) > Constants::SMALL)
        w = log(total_power0/total_power1)/dt0;
      else
        w = 0.0;
    }
  }
  
  for (unsigned int i=0; i < flux.size(); ++i) {
    if (extrap == Constants::linear)
    {
      w = (flux[i] - flux0[i])/dtt;
      flux[i] = flux[i] + w*dt;
    }
    else if (extrap == Constants::pointwise) 
    {
      if (fabs(flux0[i]) > Constants::SMALL && fabs(dtt) > Constants::SMALL)
        w = log(flux[i]/flux0[i])/dtt;
      else
        w = 0.0;
      flux[i] = flux[i]*exp(w*dt);
    }
    else if (extrap == Constants::global) 
    {
      flux[i] = flux[i]*exp(w*dt);
    }
    
  }
  
  return flux;

}

//======================================================================================================================
// flux_extrap
template <int dim>
vector<double> EnergyGroup<dim>::get_local_frequencies (double t) const 
{
  // Local Variables
  const unsigned int    n = n_dofs ();
  vector<double>        periods (n, 0.0);
  
  if (t > time0)
  {
    if (fabs(dt) < Constants::SMALL)
      return periods;
    else
    {
      for (unsigned int i = 0; i < n; ++i)
        if (fabs(sol0[i]) > Constants::SMALL)
          periods.at(i) = log(sol[i]/sol0[i])/dt;
    }
  }
  else if (t >= time1 && t < time0)
  {
    if (fabs(dt0) < Constants::SMALL)
      return periods;
    else
    {
      for (unsigned int i = 0; i < n; ++i)
        if (fabs(sol1[i]) > Constants::SMALL)
          periods.at(i) = log(sol0[i]/sol1[i])/dt0;
    }
  }
  
  return periods;
}

//======================================================================================================================
// get_fiss_frac
template <int dim>
double EnergyGroup<dim>::get_fiss_frac () const
{

  double frac (1.0);
  double bet;
  double lam;
  
  double dt_prec = (time + dtnext) - prec_time;
  unsigned int n_prec = material_data.get_n_prec ();

  switch (prob_type) {
    case Constants::steady:

        frac = 1.0;
        break;

    case Constants::transient:
        
        frac = 1.0;
        for (unsigned int i=0; i < n_prec; ++i) {
            bet = material_data.get_beta (i);
            lam = material_data.get_lambda (i);
            frac += lam*dt_prec/(1.0 + lam*dt_prec)*bet - bet;
        }
        break;
  }
  return frac;
}

//======================================================================================================================
// calculate dtest 
template <int dim>
void EnergyGroup<dim>::calculate_dtest () 
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
double EnergyGroup<dim>::calculate_maximum_absolute_secder_stencil () 
{
    if (dt < Constants::SMALL || dt0 < Constants::SMALL) return 0.0;
    
    double sdermax, sderest;
    const double vel = material_data.get_velocity (group);
    sdermax = 0.0; sderest = 0.0;
    
    for (unsigned int i = 0; i < sol.size(); ++i)
    {
        sderest = (dt0*sol[i] - (dt+dt0)*sol0[i] + dt*sol1[i])/(dt*dt*dt0);
        sderest /= vel;
        if (fabs(sderest) > sdermax)
          sdermax = fabs(sderest);
    }
    
    return sderest;
}

//======================================================================================================================
// calculate maximum relative second derivative 
template <int dim>
double EnergyGroup<dim>::calculate_maximum_relative_secder_stencil () 
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
double EnergyGroup<dim>::calculate_minimum_local_period () 
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

#endif
