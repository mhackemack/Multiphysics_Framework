#ifndef material_data_h_
#define material_data_h_

//Global Headers
#include "../Global/GlobalHeaders.h"
#include "../Global/Deal2Headers.h"
#include "../Global/Constants.h"

class MaterialData
  {
    public:
      // Constructor
      MaterialData (std::string fname) : folder_path (fname) {set_data ();};
      
      // Load Information
      void set_data ();
      void set_kinetics_data ();
      void set_temp_dep_kinetics_data ();
      
      // Temperature Dependent Cross-Section Accessors
      double get_diffusion_coefficient (const unsigned int group, const unsigned int material_id) const;
      double get_removal_XS (const unsigned int group, const unsigned int material_id) const;
      double get_fission_XS (const unsigned int group, const unsigned int material_id) const;
      double get_fission_dist_XS (const unsigned int group_1, const unsigned int group_2, const unsigned int material_id) const;
      double get_scattering_XS (const unsigned int group_1, const unsigned int group_2, const unsigned int material_id) const;
      
      double get_diffusion_coefficient (const unsigned int group, const unsigned int material_id, const double temp) const;
      double get_removal_XS (const unsigned int group, const unsigned int material_id, const double temp) const;
      double get_fission_XS (const unsigned int group, const unsigned int material_id, const double temp) const;
      double get_fission_dist_XS (const unsigned int group_1, const unsigned int group_2, const unsigned int material_id, const double temp) const;
      double get_scattering_XS (const unsigned int group_1, const unsigned int group_2, const unsigned int material_id, const double temp) const;
      
      // Temperature Independent Cross-Sections
      double get_fission_spectrum (const unsigned int group, const unsigned int material_id) const;
      double get_chi (const unsigned int group, const unsigned int material_id) const;
      double get_chi_d (const unsigned int g, const unsigned int prec, const unsigned int material_id) const;
      double get_nu (const unsigned int group, const unsigned int material_id) const;

      // Accessors
      double          get_lambda (const unsigned int prec) const;
      double          get_beta (const unsigned int prec) const;
      double          get_velocity (const unsigned int g) const;
      
      unsigned int    get_n_groups () const;
      unsigned int    get_n_materials () const;
      unsigned int    get_n_prec () const;

    private:
      unsigned int  n_groups;
      unsigned int  n_materials;
      unsigned int  n_prec;
      unsigned int  phys_type;
      bool          temp_dep;
      std::string   folder_path;
      
      ParameterHandler XS_handler;
      
      // Temperature Dependent Values
      vector<vector<Table<2,double>>>         diffusion;
      vector<vector<Table<2,double>>>         sigma_r;
      vector<vector<Table<2,double>>>         nu_sigma_f;
      vector<vector<vector<Table<2,double>>>> sigma_s;
      
      // Temperature Independent Values
      Table<2,double>           chi;
      Table<2,double>           nu;
      Table<1,double>           lambda;
      Table<1,double>           beta;
      vector<Table<2,double>>   chi_d;
      Table<1,double>           velocity;
      
      vector<double>    material_temp_count;
      
      
    private:
      double get_value (const double &val_in, const Table<2,double> &tab_in) const;
      double get_derivative (const double &val_in, const Table<2,double> &tab_in) const;

  };


//======================================================================================================================
void MaterialData::set_data ()
{
    // Use ParameterHandlers to get material information
    ParameterHandler kin_handler;
    ParameterHandler proj_handler;
    
    declare_project_parameters (proj_handler, folder_path);
    
    // Set material information
    n_materials = proj_handler.get_integer ("n_materials");
    temp_dep    = proj_handler.get_bool ("Temperature Dependence");
    phys_type   = proj_handler.get_integer ("Physics Type");
    
    if (phys_type == Constants::Physics_TH_Only) return;
    if (phys_type == Constants::Physics_Neutronics_Only) temp_dep = false;
    
    declare_kinetics_parameters (kin_handler, folder_path);
    n_groups    = kin_handler.get_integer ("n_erg");
    n_prec      = kin_handler.get_integer ("n_prec");
    
    // Load in cross-section data
    if (temp_dep)
    {
      //eclare_temp_dep_kinetics_XS (XS_handler, folder_path);
      declare_kinetics_XS (XS_handler, folder_path);
      material_temp_count.resize (n_materials);
      kin_handler.enter_subsection ("Number of Temperature Points");
      for (unsigned int m = 0; m < n_materials; ++m)
        material_temp_count.at(m) = kin_handler.get_integer ("Material " + Utilities::int_to_string(m+1,Utilities::needed_digits(m+1)) + " Point Count");
      kin_handler.leave_subsection ();
    }
    else
    {
      material_temp_count.resize (n_materials, 1);
      declare_kinetics_XS (XS_handler, folder_path);
    }
        
    // Allocate Tables
    diffusion.resize (n_materials);
    sigma_r.resize (n_materials);
    nu_sigma_f.resize (n_materials);
    sigma_s.resize (n_materials);
    
    
    for (unsigned int i = 0; i < n_materials; ++i)
    {
      diffusion.at(i).resize(n_groups);
      sigma_r.at(i).resize(n_groups);
      nu_sigma_f.at(i).resize(n_groups);
      sigma_s.at(i).resize(n_groups);
      for (unsigned int g = 0; g < n_groups; ++g)
        sigma_s.at(i).at(g).resize(n_groups);
    }
    
    for (unsigned int i = 0; i < n_materials; ++i)
    {
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        if (temp_dep)
        {
          diffusion.at(i).at(g).reinit (material_temp_count.at(i),2);
          sigma_r.at(i).at(g).reinit (material_temp_count.at(i),2);
          nu_sigma_f.at(i).at(g).reinit (material_temp_count.at(i),2);
        }
        else
        {
          diffusion.at(i).at(g).reinit (1,2);
          sigma_r.at(i).at(g).reinit (1,2);
          nu_sigma_f.at(i).at(g).reinit (1,2);
        }
        for (unsigned int gg = 0; gg < n_groups; ++gg)
          if (temp_dep)
            sigma_s.at(i).at(g).at(gg).reinit (material_temp_count.at(i),2);
          else
            sigma_s.at(i).at(g).at(gg).reinit (1,2);
      }
    }
    
    // Temperature Independent Tables
    chi.reinit (n_materials, n_groups);
    nu.reinit (n_materials, n_groups);  
    lambda.reinit (n_prec);
    beta.reinit (n_prec);
    velocity.reinit (n_groups); 
    
    // Allocate Scattering XSection - WILL BE FIXED LATER
    chi_d.resize (n_materials);
    for (unsigned int i = 0; i < n_materials; ++i)
      chi_d[i].reinit (n_prec, n_groups);
    
    // Set Cross-Section Data
    if (temp_dep)
      //set_temp_dep_kinetics_data ();
      set_kinetics_data ();
    else
      set_kinetics_data ();
    
    XS_handler.clear ();
    
    return;
}

//======================================================================================================================
void MaterialData::set_kinetics_data ()
{
  
  // Get all Material Data
  for (unsigned int m = 0; m < n_materials; ++m)
  {
  
    XS_handler.enter_subsection ("Material " + Utilities::int_to_string(m+1,Utilities::needed_digits(m+1)));
    
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        const unsigned int dig = Utilities::needed_digits(g+1);
        diffusion[m][g][0][1] = XS_handler.get_double ("Diffusion "+ Utilities::int_to_string(g+1,dig));
        sigma_r[m][g][0][1] = XS_handler.get_double ("Removal "+ Utilities::int_to_string(g+1,dig));
        nu_sigma_f[m][g][0][1] = XS_handler.get_double ("Production "+ Utilities::int_to_string(g+1,dig));
        chi[m][g] = XS_handler.get_double ("Chi "+ Utilities::int_to_string(g+1,dig));
        nu[m][g] = XS_handler.get_double ("Nu "+ Utilities::int_to_string(g+1,dig));
        for (unsigned int p = 0; p < n_prec; ++p)
          chi_d[m][p][g] = XS_handler.get_double ("Delayed Chi from Precursor " + Utilities::int_to_string(p+1,Utilities::needed_digits(p+1)) + " to Group " + Utilities::int_to_string(g+1,dig));
        for (unsigned int gg = 0; gg < n_groups; ++gg)
        {
          if (g==gg) continue;
          sigma_s[m][g][gg][0][1] = XS_handler.get_double ("Scattering "+ Utilities::int_to_string(g+1,dig) + " to " + Utilities::int_to_string(gg+1,Utilities::needed_digits(gg+1)));
        }
      }
    
    XS_handler.leave_subsection ();
  
  }
  
  // Get Remaining Transient Data
  XS_handler.enter_subsection ("Transient Data");
    
    for (unsigned int g = 0; g < n_groups; ++g)
      velocity[g] = XS_handler.get_double ("Velocity "+ Utilities::int_to_string(g+1,Utilities::needed_digits(g+1)));
    
    for (unsigned int p = 0; p < n_prec; ++p)
    {
      lambda[p] = XS_handler.get_double ("Lambda "+ Utilities::int_to_string(p+1,Utilities::needed_digits(p+1)));
      beta[p] = XS_handler.get_double ("Beta "+ Utilities::int_to_string(p+1,Utilities::needed_digits(p+1)));
    }
    
  XS_handler.leave_subsection ();
  
  
  return;
}

//======================================================================================================================
void MaterialData::set_temp_dep_kinetics_data ()
{
  // Get all Material Data
  for (unsigned int m = 0; m < n_materials; ++m)
  {
    std::string   m_num = Utilities::int_to_string(m+1, Utilities::needed_digits(m+1));
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      std::string   g_num = Utilities::int_to_string(g+1, Utilities::needed_digits(g+1));
      
      // Diffusion
      XS_handler.enter_subsection ("Material " + m_num + " Diffusion " + g_num);
      for (unsigned int i = 0; i < material_temp_count.at(m); ++i)
      {
        diffusion[m][g][i][0] = XS_handler.get_double ("Temperature " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
        diffusion[m][g][i][1] = XS_handler.get_double ("Value " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
      }
      XS_handler.leave_subsection ();
      
      // Removal
      XS_handler.enter_subsection ("Material " + m_num + " Removal " + g_num);
      for (unsigned int i = 0; i < material_temp_count.at(m); ++i)
      {
        sigma_r[m][g][i][0] = XS_handler.get_double ("Temperature " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
        sigma_r[m][g][i][1] = XS_handler.get_double ("Value " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
      }
      XS_handler.leave_subsection ();
      
      // Production
      XS_handler.enter_subsection ("Material " + m_num + " Production " + g_num);
      for (unsigned int i = 0; i < material_temp_count.at(m); ++i)
      {
        nu_sigma_f[m][g][i][0] = XS_handler.get_double ("Temperature " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
        nu_sigma_f[m][g][i][1] = XS_handler.get_double ("Value " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
      }
      XS_handler.leave_subsection ();
      
      // Scattering
      for (unsigned int gg = 0; gg < n_groups; ++gg)
      {
        if (g == gg) continue;
        std::string   gg_num = Utilities::int_to_string(gg+1, Utilities::needed_digits(gg+1));
        XS_handler.enter_subsection ("Material " + m_num + " Scattering " + g_num + " to " + gg_num);
        for (unsigned int i = 0; i < material_temp_count.at(m); ++i)
        {
          sigma_s[m][g][g][i][0] = XS_handler.get_double ("Temperature " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
          sigma_s[m][g][g][i][1] = XS_handler.get_double ("Value " + Utilities::int_to_string(i+1,Utilities::needed_digits(i+1)));
        }
      XS_handler.leave_subsection ();
      }
      
      // Temperature Independent Data
      XS_handler.enter_subsection ("Material " + m_num + " Temperature Independent Data");
      
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        std::string g_num = Utilities::int_to_string(g+1, Utilities::needed_digits(g+1));
        chi[m][g] = XS_handler.get_double ("Chi "+ g_num);
        nu[m][g] = XS_handler.get_double ("Nu "+ g_num);
        for (unsigned int p = 0; p < n_prec; ++p)
        {
          std::string p_num = Utilities::int_to_string(p+1, Utilities::needed_digits(p+1));
          chi_d[m][p][g] = XS_handler.get_double ("Delayed Chi from Precursor " + p_num + " to Group " + g_num);
        }
      }
      
      XS_handler.leave_subsection ();
      
      
    }
  }
  
  // Get Remaining Transient Data
  XS_handler.enter_subsection ("Transient Data");
    
    for (unsigned int g = 0; g < n_groups; ++g)
      velocity[g] = XS_handler.get_double ("Velocity "+ Utilities::int_to_string(g+1,Utilities::needed_digits(g+1)));
    
    for (unsigned int p = 0; p < n_prec; ++p)
    {
      lambda[p] = XS_handler.get_double ("Lambda "+ Utilities::int_to_string(p+1,Utilities::needed_digits(p+1)));
      beta[p] = XS_handler.get_double ("Beta "+ Utilities::int_to_string(p+1,Utilities::needed_digits(p+1)));
    }
    
  XS_handler.leave_subsection ();
  
  return;
}
  

//======================================================================================================================
//
// Temperature Dependent Cross-Section Accessors
//
//======================================================================================================================
double MaterialData::get_diffusion_coefficient (const unsigned int group, const unsigned int material_id, const double temp) const
{
Assert (group < n_groups, ExcIndexRange (group, 0, n_groups));
Assert (material_id < n_materials, ExcIndexRange (material_id, 0, n_materials));

if (!temp_dep || material_temp_count.at(material_id) == 1)
  return diffusion.at(material_id).at(group)[0][1];
else
  return get_value (temp,diffusion.at(material_id).at(group));
}


//======================================================================================================================
double MaterialData::get_removal_XS (const unsigned int group, const unsigned int material_id, const double temp) const
{
Assert (group < n_groups, ExcIndexRange (group, 0, n_groups));
Assert (material_id < n_materials, ExcIndexRange (material_id, 0, n_materials));

if (!temp_dep || material_temp_count.at(material_id) == 1)
  return sigma_r.at(material_id).at(group)[0][1];
else
  return get_value (temp,sigma_r.at(material_id).at(group));
}

//======================================================================================================================
double
MaterialData::get_fission_XS (const unsigned int group, const unsigned int material_id, const double temp) const
{
Assert (group < n_groups, ExcIndexRange (group, 0, n_groups));
Assert (material_id < n_materials, ExcIndexRange (material_id, 0, n_materials));

if (!temp_dep || material_temp_count.at(material_id) == 1)
  return nu_sigma_f.at(material_id).at(group)[0][1];
else
  return get_value (temp,nu_sigma_f.at(material_id).at(group));
}


//======================================================================================================================
double MaterialData::get_scattering_XS (const unsigned int group_1, const unsigned int group_2, const unsigned int material_id, const double temp) const
{
Assert (group_1 < n_groups, ExcIndexRange (group_1, 0, n_groups));
Assert (group_2 < n_groups, ExcIndexRange (group_2, 0, n_groups));
Assert (material_id < n_materials, ExcIndexRange (material_id, 0, n_materials));

if (!temp_dep || material_temp_count.at(material_id) == 1)
  return sigma_s.at(material_id).at(group_1).at(group_2)[0][1];
else
  return get_value (temp,sigma_s.at(material_id).at(group_1).at(group_2));
}

//======================================================================================================================
double MaterialData::get_fission_dist_XS (const unsigned int group_1, const unsigned int group_2, const unsigned int material_id, const double temp) const
{
return get_fission_spectrum(group_1, material_id) * get_fission_XS(group_2, material_id, temp);
//return get_fission_spectrum(group_1, material_id) * get_fission_XS(group_2, material_id);
}

//======================================================================================================================
//
// Temperature Independent Cross-Section Accessors
//
//======================================================================================================================
double MaterialData::get_fission_spectrum (const unsigned int group, const unsigned int material_id) const
{
Assert (group < n_groups, ExcIndexRange (group, 0, n_groups));
Assert (material_id < n_materials, ExcIndexRange (material_id, 0, n_materials));

return chi[material_id][group];
}

//======================================================================================================================
double MaterialData::get_chi (const unsigned int group, const unsigned int material_id) const
{
Assert (group < n_groups, ExcIndexRange (group, 0, n_groups));
Assert (material_id < n_materials, ExcIndexRange (material_id, 0, n_materials));

return chi[material_id][group];
}

//======================================================================================================================
double MaterialData::get_nu (const unsigned int group, const unsigned int material_id) const
{
Assert (group < n_groups, ExcIndexRange (group, 0, n_groups));
Assert (material_id < n_materials, ExcIndexRange (material_id, 0, n_materials));

return nu[material_id][group];
}

//======================================================================================================================
// get_lambda
double MaterialData::get_lambda (const unsigned int prec) const
{
Assert (prec < n_prec, ExcIndexRange (prec, 0, n_prec));

return lambda[prec];
}

//======================================================================================================================
// get_beta
double MaterialData::get_beta (const unsigned int prec) const
{
Assert (prec < n_prec, ExcIndexRange (prec, 0, n_prec));

return beta[prec];
}

//======================================================================================================================
// get_velocity
double MaterialData::get_velocity (const unsigned int group) const
{
Assert (group < n_groups, ExcIndexRange (group, 0, n_groups));

return velocity[group];
}

//======================================================================================================================
// get_chi_d
double MaterialData::get_chi_d (const unsigned int group, const unsigned int prec, const unsigned int material_id) const
{
Assert (group < n_groups, ExcIndexRange (group, 0, n_groups));
Assert (prec < n_prec, ExcIndexRange (prec, 0, n_prec));
Assert (material_id < n_materials, ExcIndexRange (material_id, 0, n_materials));

return chi_d[material_id][prec][group];
}

//======================================================================================================================
// get_n_groups
unsigned int MaterialData::get_n_groups () const
{
return n_groups;
}

//======================================================================================================================
// get_n_materials
unsigned int MaterialData::get_n_materials () const
{
return n_materials;
}

//======================================================================================================================
// get_n_prec
unsigned int MaterialData::get_n_prec () const
{
return n_prec;
}

//======================================================================================================================
// get_value
double MaterialData::get_value (const double &val_in, const Table<2,double> &tab_in) const
{
  // Local Variables
  double val (0.0);
  unsigned int n = tab_in.size(0);
  
  if (n == 1) return tab_in[0][1];
  
  // Linearly interpolate if input value is less than lowest in table
  // or greater than highest in table
  if (val_in <= tab_in[0][0])
  {
    val = tab_in[0][1] - (tab_in[1][1] - tab_in[0][1])/(tab_in[1][0] - tab_in[0][0]) * (tab_in[0][0] - val_in);
    return val;
  }
  else if (val_in >= tab_in[n-1][0])
  {
    val = tab_in[n-1][1] + (tab_in[n-1][1] - tab_in[n-2][1])/(tab_in[n-1][0] - tab_in[n-2][0]) * (val_in - tab_in[n-1][0]);
    return val;
  }
  
  // Loop through elements in table
  for (unsigned int i = 1; i < n; ++i)
  {
    if (val_in < tab_in[i][0])
    {
      val = tab_in[i-1][1] + (tab_in[i][1] - tab_in[i-1][1])/(tab_in[i][0] - tab_in[i-1][0]) * (val_in - tab_in[i-1][0]);
      return val;
    }
  }
  
  return val;
}

//======================================================================================================================
// get_value
double MaterialData::get_derivative (const double &val_in, const Table<2,double> &tab_in) const
{
  // Local Variables
  double val (0.0);
  unsigned int n = tab_in.size(0);
  
  if (n == 1) return 0.0;
  
  // Linearly interpolate if input value is less than lowest in table
  // or greater than highest in table
  if (val_in <= tab_in[0][0])
  {
    val = (tab_in[1][1] - tab_in[0][1])/(tab_in[1][0] - tab_in[0][0]);
    return val;
  }
  else if (val_in >= tab_in[n-1][0])
  {
    val = (tab_in[n-1][1] - tab_in[n-2][1])/(tab_in[n-1][0] - tab_in[n-2][0]);
    return val;
  }
  
  // Loop through elements in table
  for (unsigned int i = 1; i < n; ++i)
  {
    if (val_in < tab_in[i][0])
    {
      val = (tab_in[i][1] - tab_in[i-1][1])/(tab_in[i][0] - tab_in[i-1][0]);
      return val;
    }
  }
  
  return val;
  
  return val;
}

//======================================================================================================================
#endif

