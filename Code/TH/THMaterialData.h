//======================================================================================================================
/*!
 * \class Thermal Hydraulics Material Data
 * \brief
 *
 * Routine to print out heading and various outputs.
 */
//======================================================================================================================

#ifndef THmaterial_data_h_
#define THmaterial_data_h_

//Global Headers
#include "../Global/GlobalHeaders.h"
#include "../Global/Deal2Headers.h"
#include "../Global/Constants.h"



class THMaterialData
{
    public:
    
    // Constructor
    THMaterialData (const std::string path);
    
    //Destructor
    ~THMaterialData () {};
    
    // Accessors
    double get_k (const unsigned int mat, const double temp) const;
    double get_pCp (const unsigned int mat, const double temp) const;
    double get_Cp (const unsigned int mat, const double temp) const;
    double get_dens (const unsigned int mat, const double temp) const;
    
    double get_k_der (const unsigned int mat, const double temp) const;
    double get_pCp_der (const unsigned int mat, const double temp) const;
    double get_Cp_der (const unsigned int mat, const double temp) const;
    double get_dens_der (const unsigned int mat, const double temp) const;
    
    //Manipulators
    void set_material_properties ();
    void set_properties ();
    void set_temp_dep_properties ();
    
    private:
    ParameterHandler    TH_mat_prop;
    std::string         folder_path;
    unsigned int        n_materials;
    unsigned int        phys_type;
    bool                temp_dep;
    
    vector<double>  k_temp_count;
    vector<double>  Cp_temp_count;
    vector<double>  dens_temp_count;
    
    // Material Data Tables
    vector< vector<double>> k_val;
    vector< vector<double>> k_temp;
    vector< vector<double>> pCp_val;
    vector< vector<double>> pCp_temp;
    vector< Table<2,double>> thermal_conductivity;
    vector< Table<2,double>> heat_capacity;
    vector< Table<2,double>> density;
};

//======================================================================================================================
// THMaterialData Constructor
THMaterialData::THMaterialData (const std::string name) 
: folder_path (name)  
{
  set_material_properties ();
}

//======================================================================================================================
// Get Energy Density Coefficient
double THMaterialData::get_k (const unsigned int mat, const double temp) const
{
    // Local Variables
    double val (0.0);
    unsigned int n = k_temp_count.at(mat);
    
    if (!temp_dep || n == 1)
    {
        val = thermal_conductivity[mat][0][1];
        return val;
    }
    
    // Linearly interpolate if input value is less than lowest in table
    // or greater than highest in table
    if (temp <= thermal_conductivity[mat][0][0])
    {
      double temp_left = thermal_conductivity[mat][0][0];
      double temp_right = thermal_conductivity[mat][1][0];
      double val_left = thermal_conductivity[mat][0][1];
      double val_right = thermal_conductivity[mat][1][1];
      
      val = val_left - (val_right - val_left)/(temp_right - temp_left) * (temp_left - temp);
      return val;
    }
    else if (temp >= thermal_conductivity[mat][n-1][0])
    {
      double temp_left = thermal_conductivity[mat][n-2][0];
      double temp_right = thermal_conductivity[mat][n-1][0];
      double val_left = thermal_conductivity[mat][n-2][1];
      double val_right = thermal_conductivity[mat][n-1][1];
      
      val = val_right + (val_right - val_left)/(temp_right - temp_left) * (temp - temp_right);
      return val;
    }
    
    for (unsigned int i=1; i<n; ++i) {
        if (temp <= thermal_conductivity[mat][i][0]) {
            double temp_left = thermal_conductivity[mat][i-1][0];
            double temp_right = thermal_conductivity[mat][i][0];
            double val_left = thermal_conductivity[mat][i-1][1];
            double val_right = thermal_conductivity[mat][i][1];
            
            val = val_left + (val_right - val_left)/(temp_right - temp_left) * (temp - temp_left);
            return val;
        }
    }
    
    return val;
}

//======================================================================================================================
// Get pCp Coefficient
double THMaterialData::get_pCp (const unsigned int mat, const double temp) const
{
    // Local Variables
    double val (0.0);
    double Cp, dens;
    
    Cp = get_Cp (mat, temp);
    dens = get_dens (mat, temp);
    
    val = Cp * dens;
    
    return val;
}

//======================================================================================================================
// Get Heat Capacity Coefficient
double THMaterialData::get_Cp (const unsigned int mat, const double temp) const
{
    // Local Variables
    double val (0.0);
    unsigned int n = Cp_temp_count.at(mat);
    
    if (!temp_dep || n == 1)
    {
        val = heat_capacity[mat][0][1];
        return val;
    }
    
    // Linearly interpolate if input value is less than lowest in table
    // or greater than highest in table
    if (temp <= heat_capacity[mat][0][0])
    {
      double temp_left = heat_capacity[mat][0][0];
      double temp_right = heat_capacity[mat][1][0];
      double val_left = heat_capacity[mat][0][1];
      double val_right = heat_capacity[mat][1][1];
      
      val = val_left - (val_right - val_left)/(temp_right - temp_left) * (temp_left - temp);
      return val;
    }
    else if (temp >= heat_capacity[mat][n-1][0])
    {
      double temp_left = heat_capacity[mat][n-2][0];
      double temp_right = heat_capacity[mat][n-1][0];
      double val_left = heat_capacity[mat][n-2][1];
      double val_right = heat_capacity[mat][n-1][1];
      
      val = val_right + (val_right - val_left)/(temp_right - temp_left) * (temp - temp_right);
      return val;
    }
    
    for (unsigned int i=1; i<n; ++i) {
        if (temp <= heat_capacity[mat][i][0]) {
            double temp_left = heat_capacity[mat][i-1][0];
            double temp_right = heat_capacity[mat][i][0];
            double val_left = heat_capacity[mat][i-1][1];
            double val_right = heat_capacity[mat][i][1];
            
            val = val_left + (val_right - val_left)/(temp_right - temp_left) * (temp - temp_left);
            return val;
        }
    }
    
    return val;
}

//======================================================================================================================
// Get Density Coefficient
double THMaterialData::get_dens (const unsigned int mat, const double temp) const
{
    // Local Variables
    double val (0.0);
    unsigned int n = dens_temp_count.at(mat);
    
    if (!temp_dep || n == 1)
    {
        val = density[mat][0][1];
        return val;
    }
    
    // Linearly interpolate if input value is less than lowest in table
    // or greater than highest in table
    if (temp <= density[mat][0][0])
    {
      double temp_left = density[mat][0][0];
      double temp_right = density[mat][1][0];
      double val_left = density[mat][0][1];
      double val_right = density[mat][1][1];
      
      val = val_left - (val_right - val_left)/(temp_right - temp_left) * (temp_left - temp);
      return val;
    }
    else if (temp >= density[mat][n-1][0])
    {
      double temp_left = density[mat][n-2][0];
      double temp_right = density[mat][n-1][0];
      double val_left = density[mat][n-2][1];
      double val_right = density[mat][n-1][1];
      
      val = val_right + (val_right - val_left)/(temp_right - temp_left) * (temp - temp_right);
      return val;
    }
    
    for (unsigned int i=1; i<n; ++i) {
        if (temp <= density[mat][i][0]) {
            double temp_left = density[mat][i-1][0];
            double temp_right = density[mat][i][0];
            double val_left = density[mat][i-1][1];
            double val_right = density[mat][i][1];
            
            val = val_left + (val_right - val_left)/(temp_right - temp_left) * (temp - temp_left);
            
            return val;
        }
    }
    
    return val;
}

//======================================================================================================================
// Get Thermal Conductivity Coefficient Derivative
double THMaterialData::get_k_der (const unsigned int mat, const double temp) const
{
    // Local Variables
    double val (0.0);
    unsigned int n = k_temp_count.at(mat);
    
    if (!temp_dep || n == 1)
    {
        val = 0.0;
        return val;
    }
    
    // Linearly interpolate if input value is less than lowest in table
    // or greater than highest in table
    if (temp <= thermal_conductivity[mat][0][0])
    {
      double temp_left = thermal_conductivity[mat][0][0];
      double temp_right = thermal_conductivity[mat][1][0];
      double val_left = thermal_conductivity[mat][0][1];
      double val_right = thermal_conductivity[mat][1][1];
      
      val = (val_right - val_left)/(temp_right - temp_left);
      return val;
    }
    else if (temp >= thermal_conductivity[mat][n-1][0])
    {
      double temp_left = thermal_conductivity[mat][n-2][0];
      double temp_right = thermal_conductivity[mat][n-1][0];
      double val_left = thermal_conductivity[mat][n-2][1];
      double val_right = thermal_conductivity[mat][n-1][1];
      
      val = (val_right - val_left)/(temp_right - temp_left);
      return val;
    }
    
    for (unsigned int i=1; i<n; ++i) {
        if (temp <= thermal_conductivity[mat][i][0]) {
            double temp_left = thermal_conductivity[mat][i-1][0];
            double temp_right = thermal_conductivity[mat][i][0];
            double val_left = thermal_conductivity[mat][i-1][1];
            double val_right = thermal_conductivity[mat][i][1];
            
            val = (val_right - val_left)/(temp_right - temp_left);
            return val;
        }
    }
    
    return val;
}

//======================================================================================================================
// Get Thermal Conduction Coefficient Derivative
double THMaterialData::get_Cp_der (const unsigned int mat, const double temp) const
{
    // Local Variables
    double val (0.0);
    unsigned int n = Cp_temp_count.at(mat);
    
    if (!temp_dep || n == 1)
    {
        val = 0.0;
        return val;
    }
    
    // Linearly interpolate if input value is less than lowest in table
    // or greater than highest in table
    if (temp <= heat_capacity[mat][0][0])
    {
      double temp_left = heat_capacity[mat][0][0];
      double temp_right = heat_capacity[mat][1][0];
      double val_left = heat_capacity[mat][0][1];
      double val_right = heat_capacity[mat][1][1];
      
      val = (val_right - val_left)/(temp_right - temp_left);
      return val;
    }
    else if (temp >= heat_capacity[mat][n-1][0])
    {
      double temp_left = heat_capacity[mat][n-2][0];
      double temp_right = heat_capacity[mat][n-1][0];
      double val_left = heat_capacity[mat][n-2][1];
      double val_right = heat_capacity[mat][n-1][1];
      
      val = (val_right - val_left)/(temp_right - temp_left);
      return val;
    }
    
    for (unsigned int i=1; i<n; ++i) {
        if (temp <= heat_capacity[mat][i][0]) {
            double temp_left = heat_capacity[mat][i-1][0];
            double temp_right = heat_capacity[mat][i][0];
            double val_left = heat_capacity[mat][i-1][1];
            double val_right = heat_capacity[mat][i][1];
            
            val = (val_right - val_left)/(temp_right - temp_left);
            return val;
        }
    }
    
    return val;
}

//======================================================================================================================
// Get Density Coefficient Derivative
double THMaterialData::get_dens_der (const unsigned int mat, const double temp) const
{
    // Local Variables
    double val (0.0);
    unsigned int n = dens_temp_count.at(mat);
    
    if (!temp_dep || n == 1)
    {
        val = 0.0;
        return val;
    }
    
    // Linearly interpolate if input value is less than lowest in table
    // or greater than highest in table
    if (temp <= density[mat][0][0])
    {
      double temp_left = density[mat][0][0];
      double temp_right = density[mat][1][0];
      double val_left = density[mat][0][1];
      double val_right = density[mat][1][1];
      
      val = (val_right - val_left)/(temp_right - temp_left);
      return val;
    }
    else if (temp >= density[mat][n-1][0])
    {
      double temp_left = density[mat][n-2][0];
      double temp_right = density[mat][n-1][0];
      double val_left = density[mat][n-2][1];
      double val_right = density[mat][n-1][1];
      
      val = (val_right - val_left)/(temp_right - temp_left);
      return val;
    }
    
    for (unsigned int i=1; i<n; ++i) {
        if (temp <= density[mat][i][0]) {
            double temp_left = density[mat][i-1][0];
            double temp_right = density[mat][i][0];
            double val_left = density[mat][i-1][1];
            double val_right = density[mat][i][1];
            
            val = (val_right - val_left)/(temp_right - temp_left);
            return val;
        }
    }
    
    return val;
}

//======================================================================================================================
// Get Density Coefficient Derivative
double THMaterialData::get_pCp_der (const unsigned int mat, const double temp) const
{
    // Local Variables
    double val (0.0);
    
    val = get_Cp_der(mat,temp)*get_dens(mat,temp) + get_Cp(mat,temp)*get_dens_der(mat,temp);
    
    return val;
}

//======================================================================================================================
// Set Material Properties
void THMaterialData::set_material_properties ()
{
    
    // Use ParameterHandlers to get material information
    ParameterHandler proj_handler;
    ParameterHandler TH_handler;
    
    declare_project_parameters (proj_handler, folder_path);
    
    // Set material information
    n_materials = proj_handler.get_integer ("n_materials");
    temp_dep    = proj_handler.get_bool ("Temperature Dependence");
    phys_type   = proj_handler.get_integer ("Physics Type");
    
    if (phys_type == Constants::Physics_Neutronics_Only) return;
    
    declare_TH_parameters (TH_handler, folder_path);
    
    // Allocate Arrays
    thermal_conductivity.resize (n_materials);
    heat_capacity.resize (n_materials);
    density.resize (n_materials);
    k_temp_count.resize (n_materials);
    Cp_temp_count.resize (n_materials);
    dens_temp_count.resize (n_materials);
    
    for (unsigned int m = 0; m < n_materials; ++m)
    {
      if (temp_dep)
      {
        TH_handler.enter_subsection ("Number of Temperature Points");
          k_temp_count.at(m) = TH_handler.get_integer ("Material " + Utilities::int_to_string(m+1,Utilities::needed_digits(m+1)) + " Thermal Conductivity Point Count");
          Cp_temp_count.at(m) = TH_handler.get_integer ("Material " + Utilities::int_to_string(m+1,Utilities::needed_digits(m+1)) + " Heat Capacity Point Count");
          dens_temp_count.at(m) = TH_handler.get_integer ("Material " + Utilities::int_to_string(m+1,Utilities::needed_digits(m+1)) + " Density Point Count");
        TH_handler.leave_subsection ();
        
        thermal_conductivity[m].reinit (k_temp_count.at(m), 2);
        heat_capacity[m].reinit (Cp_temp_count.at(m), 2);
        density[m].reinit (dens_temp_count.at(m), 2);
      }
      else 
      {
        k_temp_count.at(m) = 1;
        Cp_temp_count.at(m) = 1;
        dens_temp_count.at(m) = 1;
        
        thermal_conductivity[m].reinit (1, 2);
        heat_capacity[m].reinit (1, 2);
        density[m].reinit (1, 2);
      }
    }
    
    // Declare TH Data
    if (temp_dep)
    {
      declare_temp_dep_TH_data (TH_mat_prop, folder_path);
      set_temp_dep_properties ();
    }
    else
    {
      declare_TH_data (TH_mat_prop, folder_path);
      set_properties ();
    }
    
    // Clean up and return
    TH_mat_prop.clear ();
    return;
}

//======================================================================================================================
// Set Non-Temperature Dependent Properties
void THMaterialData::set_properties ()
{
    // Local Variables
    
    for (unsigned int m = 0; m < n_materials; ++m)
    {
      TH_mat_prop.enter_subsection ("Material " + Utilities::int_to_string(m+1,Utilities::needed_digits(m+1)));
        thermal_conductivity[m][0][0] = 1.0;
        thermal_conductivity[m][0][1] = TH_mat_prop.get_double ("Thermal Conductivity");
        heat_capacity[m][0][0] = 1.0;
        heat_capacity[m][0][1] = TH_mat_prop.get_double ("Heat Capacity");
        density[m][0][0] = 1.0;
        density[m][0][1] = TH_mat_prop.get_double ("Density");
      TH_mat_prop.leave_subsection ();
    }
  
    return;
}

//======================================================================================================================
// Set Temperature Dependent Properties
void THMaterialData::set_temp_dep_properties ()
{
    // Local Variables
    
    for (unsigned int m = 0; m < n_materials; ++m)
    {
      std::string           m_num = Utilities::int_to_string(m+1, Utilities::needed_digits(m+1));
      const unsigned int    k_num = k_temp_count.at(m);
      const unsigned int    Cp_num = Cp_temp_count.at(m);
      const unsigned int    dens_num = dens_temp_count.at(m);
      
      // Thermal Conductivity
      TH_mat_prop.enter_subsection ("Material " + m_num + " Thermal Conductivity");
        
        for (unsigned int i = 0; i < k_num; ++i)
        {
          thermal_conductivity[m][i][0] = TH_mat_prop.get_double ("Temperature " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
          thermal_conductivity[m][i][1] = TH_mat_prop.get_double ("Value " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
        }
        
      TH_mat_prop.leave_subsection ();
      
      // Heat Capacity
      TH_mat_prop.enter_subsection ("Material " + m_num + " Heat Capacity");
        
        for (unsigned int i = 0; i < Cp_num; ++i)
        {
          heat_capacity[m][i][0] = TH_mat_prop.get_double ("Temperature " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
          heat_capacity[m][i][1] = TH_mat_prop.get_double ("Value " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
        }
      
      TH_mat_prop.leave_subsection ();
      
      // Density
      TH_mat_prop.enter_subsection ("Material " + m_num + " Density");
        
        for (unsigned int i = 0; i < dens_num; ++i)
        {
          density[m][i][0] = TH_mat_prop.get_double ("Temperature " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
          density[m][i][1] = TH_mat_prop.get_double ("Value " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
        }
        
      TH_mat_prop.leave_subsection ();
      
      
    }
  
    return;
}

//======================================================================================================================

#endif
