// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Driver function for a simple beam problem

// OOMPH-LIB includes
#include "generic.h"
#include "beam.h"
#include "meshes/one_d_lagrangian_mesh.h"

using namespace std;
using namespace oomph;

//========start_of_namespace========================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
  /// Non-dimensional thickness
  double H = 0.0;

  /// Non-dimensional coefficient (FSI)
  double I = 0.0;

  /// Angle between the two arms of the beam
  double Alpha = 0.0;


  // These are parameters that can be set from the command line:

  /// Aspect ratio: Don't change this on the fly; should only be assigned
  /// once before the mesh is generated
  /// first arm length = |q+0.5|, second arm length = |q-0.5|
  // double Q = 0.4;

  // Define the length of the beam in the GeomObejct
  double Stretch_ratio = 0.4;

  /// Initial value for theta_eq in the Newton solve
  double Initial_value_for_theta_eq = -1.57;

  /// Default value for desired ds
  double Ds_default = 1.0e-5;

  // To prevent large solution jumps in critical intervals for I, try to reduce
  // ds specifically in those areas for smoother and more precise results

  /// End point for the first I interval [Interval1_start,Interval1_end]
  double Interval1_start = 0.0;
  double Interval1_end = 0.01;

  /// Start and end points for the second I interval
  /// [Interval2_start,Interval2_end]
  double Interval2_start = 0.07;
  double Interval2_end = 0.08;

  // If the interval_ds is smaller than the default, it will automatically
  // revert to using the default interval

  /// Value of ds for first interval
  double Ds_interval1 = 10.0;

  /// Value of ds for second interval
  double Ds_interval2 = 10.0;


} // namespace Global_Physical_Variables


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=========================================================================
/// RigidBodyElement
//=========================================================================
class RigidBodyElement : public GeneralisedElement
{
public:
  /// Constructor: Pass initial values for rigid body parameters (pinned
  /// by default) and time stepper pointer
  RigidBodyElement(const double& X0,
                   const double& Y0,
                   const double& Theta_eq,
                   TimeStepper*& time_stepper_pt)
    : Time_stepper_pt(0)
  {
    // Pass time stepper pointer
    Time_stepper_pt = time_stepper_pt;

    // Create internal data which contains the "rigid body" parameters
    for (unsigned i = 0; i < 3; i++)
    {
      // Create data: One value, timedependence, free by default
      add_internal_data(new Data(Time_stepper_pt, 1));
    }

    // Give initial values t=0
    internal_data_pt(0)->set_value(0, 0, X0);
    internal_data_pt(1)->set_value(0, 0, Y0);
    internal_data_pt(2)->set_value(0, 0, Theta_eq);

    /* // These are just initial values so pin
    internal_data_pt(3)->set_value(0, X0);
    internal_data_pt(3)->pin(0);
    internal_data_pt(4)->set_value(0, Y0);
    internal_data_pt(4)->pin(0); */
  }


  /// Function that returns the Vector of pointers to the "rigid body"
  /// parameters
  Vector<Data*> rigid_body_parameters()
  {
    Vector<Data*> tmp_pt(3);
    for (unsigned i = 0; i < 3; i++)
    {
      tmp_pt[i] = internal_data_pt(i);
    }
    return tmp_pt;
  }


  /// Helper function to compute the meaningful parameter values
  /// from enumerated data
  void get_parameters(double& X0, double& Y0, double& Theta_eq)
  {
    X0 = internal_data_pt(0)->value(0);
    Y0 = internal_data_pt(1)->value(0);
    Theta_eq = internal_data_pt(2)->value(0);
  }


  /// Helper function to compute the first time derivative of the meaningful
  /// parameter values
  void get_first_time_derivative_of_parameters(double& dX0_dt,
                                               double& dY0_dt,
                                               double& dTheta_eq_dt)
  {
    dX0_dt = Time_stepper_pt->time_derivative(1, internal_data_pt(0), 0);
    dY0_dt = Time_stepper_pt->time_derivative(1, internal_data_pt(1), 0);
    dTheta_eq_dt = Time_stepper_pt->time_derivative(1, internal_data_pt(2), 0);
  }


  /// Pass pointer to the Mesh of HaoHermiteBeamElements
  /// and add their unknowns to be external data for this element
  void set_pointer_to_beam_meshes(const Vector<SolidMesh*>& beam_mesh_pt)
  {
    // Store the pointer for future reference
    Beam_mesh_pt = beam_mesh_pt;

    // Loop over the nodes in the all mesh and add them as external Data
    // because they affect the traction and therefore the total drag
    // and torque on the object
    unsigned npointer = beam_mesh_pt.size();
    for (unsigned i = 0; i < npointer; i++)
    {
      unsigned nnode = beam_mesh_pt[i]->nnode();
      for (unsigned j = 0; j < nnode; j++)
      {
        add_external_data(beam_mesh_pt[i]->node_pt(j)->variable_position_pt());
      }
    }
  }


  /// Compute the beam's centre of mass
  void compute_centre_of_mass(Vector<double>& sum_r_centre);


  /// Compute the drag and torque on the entire beam structure according
  /// to slender body theory
  void compute_drag_and_torque(Vector<double>& sum_total_drag,
                               double& sum_total_torque);


  /// Output the Theta_eq, Theta_eq_orientation (make comparision with paper's
  /// results), drag and torque on the entire beam structure
  void output(std::ostream& outfile)
  {
    Vector<double> sum_total_drag(2);
    double sum_total_torque = 0.0;

    // Compute the drag and torque on the entire beam structure
    compute_drag_and_torque(sum_total_drag, sum_total_torque);

    // Output Theta_eq
    double Theta_eq = internal_data_pt(2)->value(0);
    outfile << fmod(Theta_eq, 2.0 * acos(-1.0)) << "  ";

    // Make a transformation from Theta_eq to Theta_eq_orientation
    // Note that here Theta_eq_orientation is controlled in the range of
    // [-2*PI,2*PI]
    double Theta_eq_orientation = fmod(
      fmod(Theta_eq, 2.0 * acos(-1.0)) + acos(-1.0) / 2.0, 2.0 * acos(-1.0));

    // To escape the jump of the solutions
    if (fabs(Theta_eq_orientation) > 1.5 * acos(-1.0))
    {
      if (Theta_eq_orientation > 0)
      {
        outfile << Theta_eq_orientation - 2.0 * acos(-1.0) << "  ";
      }
      else
      {
        outfile << Theta_eq_orientation + 2.0 * acos(-1.0) << "  ";
      }
    }
    else
    {
      outfile << Theta_eq_orientation << "  ";
    }

    // Output X0,Y0
    double X0 = internal_data_pt(0)->value(0);
    double Y0 = internal_data_pt(1)->value(0);
    outfile << X0 << "  ";
    outfile << Y0 << "  ";

    // Output drag and torque on the entire beam structure
    outfile << sum_total_drag[0] << "  ";
    outfile << sum_total_drag[1] << "  ";
    outfile << sum_total_torque << "  ";
  }

protected:
  // Fill in contribution to residuals
  void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
    // oomph_info << "ndof in element: " << residuals.size() << std::endl;

    // Get current total drag and torque
    Vector<double> sum_total_drag(2);
    double sum_total_torque = 0.0;
    compute_drag_and_torque(sum_total_drag, sum_total_torque);

    // current time
    // double time = Time_stepper_pt->time();

    unsigned n_internal = ninternal_data();
    for (unsigned i = 0; i < n_internal; i++)
    {
      // Get the local equation number of the zeroth dof
      // associated with this internal Data object
      unsigned j = 0;
      int eqn_number = internal_local_eqn(i, j);

      // Is it an actual dof
      if (eqn_number >= 0)
      {
        if (i == 0)
        {
          // Eqn for V:
          residuals[eqn_number] = sum_total_drag[0];
          // residuals[eqn_number] = internal_data_pt(i)->value(j) - time;
          // residuals[eqn_number] = internal_data_pt(i)->value(j) - time;
        }
        else if (i == 1)
        {
          // Eqn for U0:
          residuals[eqn_number] = sum_total_drag[1];
          // residuals[eqn_number] = internal_data_pt(i)->value(j) - sin(time);
          // residuals[eqn_number] = internal_data_pt(i)->value(j) - time *
          // time;
        }
        else if (i == 2)
        {
          // Eqn for Theta_eq:
          residuals[eqn_number] = sum_total_torque;
          // residuals[eqn_number] =
          //  internal_data_pt(i)->value(j) - time * atan(1.0);
        }
        else
        {
          oomph_info << "Never get here\n";
          abort();
        }

        // std::cout << "internal data " << i << " is not pinned\n";
      }
      else
      {
        // std::cout << "internal data " << i << " is pinned\n";
      }
    }
  }

private:
  /// Pointer to the Mesh of HaoHermiteBeamElements
  Vector<SolidMesh*> Beam_mesh_pt;

  /// Time stepper pointer
  TimeStepper* Time_stepper_pt;
};


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


//=====================================================================
/// Upgraded Hermite Beam Element to incorporate slender body traction
//=====================================================================
class HaoHermiteBeamElement : public virtual HermiteBeamElement
{
public:
  /// Constructor: Initialise private member data
  HaoHermiteBeamElement()
    : Rigid_body_element_pt(0), I_pt(0), Theta_initial_pt(0)
  {
  }


  /// Pass pointer to RigidBodyElement that contains the rigid body parameters
  void set_pointer_to_rigid_body_element(
    RigidBodyElement* rigid_body_element_pt)
  {
    // Store the pointer for future reference
    Rigid_body_element_pt = rigid_body_element_pt;

    // Get the rigid body parameters
    Vector<Data*> rigid_body_data_pt =
      Rigid_body_element_pt->rigid_body_parameters();

#ifdef PARANOID
    if (rigid_body_data_pt.size() != 3)
    {
      std::ostringstream error_message;
      error_message << "rigid_body_data_pt should have size 3, not "
                    << rigid_body_data_pt.size() << std::endl;

      // loop over all entries
      for (unsigned i = 0; i < 3; i++)
      {
        if (rigid_body_data_pt[i]->nvalue() != 1)
        {
          error_message << "rigid_body_data_pt[" << i
                        << "] should have 1 value, not "
                        << rigid_body_data_pt[i]->nvalue() << std::endl;
        }
      }

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Add the rigid body parameters as the external data for this element
    for (unsigned i = 0; i < 3; i++)
    {
      add_external_data(rigid_body_data_pt[i]);
    }
  }


  /// Pointer to non-dimensional coefficient (FSI)
  double*& i_pt()
  {
    return I_pt;
  }


  /// Pointer to initial angle
  void theta_initial_pt(const double* theta_initial_pt)
  {
    Theta_initial_pt = theta_initial_pt;
  }

  /// Initial angle
  double theta_initial() const
  {
    if (Theta_initial_pt == 0)
    {
      return 0.0;
    }
    else
    {
      return *Theta_initial_pt;
    }
  }


  /// Compute the element's contribution to the (\int r ds) and length of beam
  void compute_contribution_to_int_r_and_length(Vector<double>& int_r,
                                                double& length)
  {
#ifdef PARANOID
    if (int_r.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "int_r should have size 2, not " << int_r.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Initialise
    int_r[0] = 0.0;
    int_r[1] = 0.0;
    length = 0.0;

    // Local coordinate (1D!)
    Vector<double> s(1);

    // Set # of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Return local coordinate s[j]  of i-th integration point.
      unsigned j = 0;
      s[j] = integral_pt()->knot(ipt, j);

      // Get position vector to and non-unit tangent vector on wall:
      // dr/ds. NOTE: This is before we apply the rigid body motion!
      // so in terms of the write-up the position vector is R_0
      Vector<double> R_0(2);
      Vector<double> drds(2);
      get_non_unit_tangent(s, R_0, drds);

      // Jacobian of mapping between local and global coordinates
      double J = sqrt(drds[0] * drds[0] + drds[1] * drds[1]);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Translate rigid body parameters into meaningful variables
      double X0 = 0.0;
      double Y0 = 0.0;
      double Theta_eq = 0.0;
      Rigid_body_element_pt->get_parameters(X0, Y0, Theta_eq);


      // hierher use Theta_initial everywhere whenever you're processing
      // Theta_eq

      // Apply rigid body translation and rotation to get the actual
      // shape of the deformed body in the fluid
      Vector<double> R(2);
      R[0] = cos(Theta_eq + theta_initial()) * R_0[0] -
             sin(Theta_eq + theta_initial()) * R_0[1] + X0;
      R[1] = sin(Theta_eq + theta_initial()) * R_0[0] +
             cos(Theta_eq + theta_initial()) * R_0[1] + Y0;

      // Add 'em.
      length += W;
      int_r[0] += R[0] * W;
      int_r[1] += R[1] * W;
    }
  }


  /// Compute the slender body traction acting on the actual beam onto the
  /// element at local coordinate s
  void compute_slender_body_traction_on_actual_beam(const Vector<double>& s,
                                                    Vector<double>& traction)
  {
#ifdef PARANOID
    if (traction.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "traction should have size 2, not " << traction.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get the Eulerian position and the unit normal.
    // NOTE: This is before we apply the rigid body motion!
    // so in terms of the write-up the position vector is R_0 and N_0
    Vector<double> R_0(2);
    Vector<double> N_0(2);
    get_normal(s, R_0, N_0);

    // Translate rigid body parameters into meaningful variables
    double X0 = 0.0;
    double Y0 = 0.0;
    double Theta_eq = 0.0;
    Rigid_body_element_pt->get_parameters(X0, Y0, Theta_eq);


    // Compute R which is after translation and rotation
    Vector<double> R(2);
    R[0] = cos(Theta_eq + theta_initial()) * R_0[0] -
           sin(Theta_eq + theta_initial()) * R_0[1] + X0;
    R[1] = sin(Theta_eq + theta_initial()) * R_0[0] +
           cos(Theta_eq + theta_initial()) * R_0[1] + Y0;

    // Compute normal N which is after translation and rotation
    Vector<double> N(2);
    N[0] = cos(Theta_eq + theta_initial()) * N_0[0] -
           sin(Theta_eq + theta_initial()) * N_0[1];
    N[1] = sin(Theta_eq + theta_initial()) * N_0[0] +
           cos(Theta_eq + theta_initial()) * N_0[1];


    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_dofs = nnodal_position_type();

    Vector<double> dR0_dt(n_dim);

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_dofs);

    // Get shape functions
    shape(s, psi);

    Vector<double> interpolated_xi(n_lagrangian);
    interpolated_xi[0] = 0.0;

    // Loop over coordinate directions/components of Vector
    for (unsigned i = 0; i < n_dim; i++)
    {
      // Initialise time derivative of R_0
      dR0_dt[i] = 0.0;
    }


    // Calculate spatial derivatives
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over positional dofs
      for (unsigned k = 0; k < n_position_dofs; k++)
      {
        // Loop over Lagrangian coordinate directions [xi_gen[] are the
        // the *gen*eralised Lagrangian coordinates: node, type, direction]
        for (unsigned i = 0; i < n_lagrangian; i++)
        {
          interpolated_xi[i] +=
            raw_lagrangian_position_gen(l, k, i) * psi(l, k);
        }

        // Loop over components of the deformed position Vector
        for (unsigned i = 0; i < n_dim; i++)
        {
          dR0_dt[i] += raw_dnodal_position_gen_dt(1, l, k, i) * psi(l, k);
        }
      }
    }

    // Translate rigid body parameters into meaningful variables (first time
    // derivative)
    double dX0_dt = 0.0;
    double dY0_dt = 0.0;
    double dTheta_eq_dt = 0.0;
    Rigid_body_element_pt->get_first_time_derivative_of_parameters(
      dX0_dt, dY0_dt, dTheta_eq_dt);

    // Particle velocity from kinematics
    Vector<double> U(2);
    U[0] = dX0_dt - dTheta_eq_dt * sin(Theta_eq) * R_0[0] +
           cos(Theta_eq) * dR0_dt[0] - dTheta_eq_dt * cos(Theta_eq) * R_0[1] -
           sin(Theta_eq) * dR0_dt[1];
    U[1] = dY0_dt + dTheta_eq_dt * cos(Theta_eq) * R_0[0] +
           sin(Theta_eq) * dR0_dt[0] - dTheta_eq_dt * sin(Theta_eq) * R_0[1] +
           cos(Theta_eq) * dR0_dt[1];

    // Compute the traction onto the element at local coordinate s
    traction[0] =
      (1.0 - 0.5 * N[1] * N[1]) * (R[1] - U[0]) - 0.5 * (N[1] * N[0] * U[1]);

    traction[1] =
      0.5 * N[1] * N[0] * (R[1] - U[0]) - (1.0 - 0.5 * N[0] * N[0]) * U[1];
  }


  /// Compute the slender body traction acting on the beam in the reference
  /// configuration (i.e. without rigid body motion!) at local coordinate s
  void compute_slender_body_traction_on_beam_in_reference_configuration(
    const Vector<double>& s, Vector<double>& traction_0)
  {
#ifdef PARANOID
    if (traction_0.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "traction_0 should have size 2, not "
                    << traction_0.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Translate rigid body parameters into meaningful variables
    double X0 = 0.0;
    double Y0 = 0.0;
    double Theta_eq = 0.0;
    Rigid_body_element_pt->get_parameters(X0, Y0, Theta_eq);

    // Compute the slender body traction acting on the actual beam onto the
    // element at local coordinate s
    Vector<double> traction(2);
    compute_slender_body_traction_on_actual_beam(s, traction);

    // Rotate the traction from the actual beam back to the reference
    // configuration.
    traction_0[0] = traction[0] * cos(Theta_eq + theta_initial()) +
                    traction[1] * sin(Theta_eq + theta_initial());
    traction_0[1] = -traction[0] * sin(Theta_eq + theta_initial()) +
                    traction[1] * cos(Theta_eq + theta_initial());
  }


  // overloaded load_vector to apply the computed traction_0 (i.e. the
  // traction acting on the beam before its rigid body motion is applied)
  // including the non-dimensional coefficient I (FSI)
  void load_vector(const unsigned& intpt,
                   const Vector<double>& xi,
                   const Vector<double>& x,
                   const Vector<double>& N,
                   Vector<double>& load)
  {
    /// Return local coordinate s[j] at the specified integration point.
    Vector<double> s(1);
    unsigned j = 0;
    s[j] = integral_pt()->knot(intpt, j);

    compute_slender_body_traction_on_beam_in_reference_configuration(s, load);
    load[0] = *(i_pt()) * load[0];
    load[1] = *(i_pt()) * load[1];
  }


  // Compute the element's contribution to the total drag and torque on
  // the entire beam structure according to slender body theory
  void compute_contribution_to_drag_and_torque(Vector<double>& drag,
                                               double& torque)
  {
#ifdef PARANOID
    if (drag.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "drag should have size 2, not " << drag.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Initialise
    drag[0] = 0.0;
    drag[1] = 0.0;
    torque = 0.0;

    // Compute the beam's positon of centre of mass
    Vector<double> sum_r_centre(2);
    Rigid_body_element_pt->compute_centre_of_mass(sum_r_centre);

    // Local coordinate (1D!)
    Vector<double> s(1);

    // Set # of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Translate rigid body parameters into meaningful variables
    double X0 = 0.0;
    double Y0 = 0.0;
    double Theta_eq = 0.0;
    Rigid_body_element_pt->get_parameters(X0, Y0, Theta_eq);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      /// Return local coordinate s[j] of i-th integration point.
      unsigned j = 0;
      s[j] = integral_pt()->knot(ipt, j);


      // Get position vector to and non-unit tangent vector on wall:
      // dr/ds
      // NOTE: This is before we apply the rigid body motion!
      // so in terms of the write-up the position vector is R_0
      Vector<double> R_0(2);
      Vector<double> drds(2);
      get_non_unit_tangent(s, R_0, drds);

      // Jacobian. Since Jacobian is the same for R, still use it here.
      double J = sqrt(drds[0] * drds[0] + drds[1] * drds[1]);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Compute the slender body traction on actual beam; note this is
      // inefficient since we've already computed certain quantities that will
      // be needed in this function
      Vector<double> traction(2);
      compute_slender_body_traction_on_actual_beam(s, traction);

      // Compute R (after translation and rotation)
      Vector<double> R(2);
      R[0] = cos(Theta_eq + theta_initial()) * R_0[0] -
             sin(Theta_eq + theta_initial()) * R_0[1] + X0;
      R[1] = sin(Theta_eq + theta_initial()) * R_0[0] +
             cos(Theta_eq + theta_initial()) * R_0[1] + Y0;

      // calculate the contribution to torque
      double local_torque = (R[0] - sum_r_centre[0]) * traction[1] -
                            (R[1] - sum_r_centre[1]) * traction[0];

      // Add 'em
      drag[0] += traction[0] * W;
      drag[1] += traction[1] * W;
      torque += local_torque * W;
    }
  }


  /// Overloaded output function
  void output(std::ostream& outfile, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << std::endl;

    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_dofs = nnodal_position_type();

    Vector<double> R_0(n_dim);

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_dofs);

    // Loop over element plot points
    for (unsigned l1 = 0; l1 < n_plot; l1++)
    {
      s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

      // Get shape functions
      shape(s, psi);

      Vector<double> interpolated_xi(n_lagrangian);
      interpolated_xi[0] = 0.0;

      // Loop over coordinate directions/components of Vector
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Initialise
        R_0[i] = 0.0;
      }

      // Calculate positions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_dofs; k++)
        {
          // Loop over Lagrangian coordinate directions [xi_gen[] are the
          // the *gen*eralised Lagrangian coordinates: node, type, direction]
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] +=
              raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over components of the deformed position Vector
          for (unsigned i = 0; i < n_dim; i++)
          {
            R_0[i] += raw_dnodal_position_gen_dt(0, l, k, i) * psi(l, k);
          }
        }
      }

      // Get the normal vector N0 at each plotted point
      Vector<double> N_0(n_dim);
      get_normal(s, N_0);

      // Compute slender body traction acting on the actual beam
      Vector<double> traction(n_dim);
      compute_slender_body_traction_on_actual_beam(s, traction);

      // Compute slender body traction acting on the beam in the reference
      // configuration
      Vector<double> traction_0(n_dim);
      compute_slender_body_traction_on_beam_in_reference_configuration(
        s, traction_0);

      // Translate rigid body parameters into meaningful variables
      double X0 = 0.0;
      double Y0 = 0.0;
      double Theta_eq = 0.0;
      Rigid_body_element_pt->get_parameters(X0, Y0, Theta_eq);

      // Compute R after translation and rotation
      Vector<double> R(n_dim);
      R[0] = cos(Theta_eq + theta_initial()) * R_0[0] -
             sin(Theta_eq + theta_initial()) * R_0[1] + X0;
      R[1] = sin(Theta_eq + theta_initial()) * R_0[0] +
             cos(Theta_eq + theta_initial()) * R_0[1] + Y0;

      // Compute normal N after translation and rotation
      Vector<double> N(n_dim);
      N[0] = cos(Theta_eq + theta_initial()) * N_0[0] -
             sin(Theta_eq + theta_initial()) * N_0[1];
      N[1] = sin(Theta_eq + theta_initial()) * N_0[0] +
             cos(Theta_eq + theta_initial()) * N_0[1];

      // Output R0 which is clamped at the origin
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << R_0[i] << " ";
      }

      // Output R which is after translation and rotation
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << R[i] << " ";
      }

      // Output unit normal N0
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << N_0[i] << " ";
      }

      // Output unit normal N which is after translation and rotation
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << N[i] << " ";
      }

      // Output traction acting on the beam in the reference configuration
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << traction_0[i] << " ";
      }

      // Output traction acting on the actual beam
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << traction[i] << " ";
      }

      // Output the velocity of the background
      outfile << R[1] << "  " << 0;
      outfile << std::endl;
    }
  }

private:
  /// Pointer to element that controls the rigid body motion
  RigidBodyElement* Rigid_body_element_pt;

  /// Pointer to non-dimensional coefficient (FSI)
  double* I_pt;

  /// Pointer to initial rotation of the element when it's in its (otherwise)
  /// undeformed configuration
  const double* Theta_initial_pt;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=============================================================================
/// Compute the beam's centre of mass (defined outside class to avoid
/// forward references)
//=============================================================================
void RigidBodyElement::compute_centre_of_mass(Vector<double>& sum_r_centre)
{
#ifdef PARANOID
  if (r_centre.size() != 2)
  {
    std::ostringstream error_message;
    error_message << "r_centre should have size 2, not " << r_centre.size()
                  << std::endl;

    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Initialise
  sum_r_centre[0] = 0.0;
  sum_r_centre[1] = 0.0;
  Vector<double> int_r(2);
  double length = 0.0;

  // Find number of beam meshes
  unsigned npointer = Beam_mesh_pt.size();

  // Loop over the beam meshes to compute the centre of mass of the entire beam
  for (unsigned i = 0; i < npointer; i++)
  {
    // Initialise
    Vector<double> total_int_r(2);
    double total_length = 0.0;

    // Find number of elements in the mesh
    unsigned n_element = Beam_mesh_pt[i]->nelement();

    // Loop over the elements to compute the sum of elements' contribution to
    // the (\int r ds) and the length of beam
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast to the specific element type
      HaoHermiteBeamElement* elem_pt =
        dynamic_cast<HaoHermiteBeamElement*>(Beam_mesh_pt[i]->element_pt(e));

      // Compute contribution to the the (\int r ds) and length of beam within
      // the e-th element
      elem_pt->compute_contribution_to_int_r_and_length(int_r, length);

      // Sum the elements' contribution to the (\int r ds) and length of beam
      total_int_r[0] += int_r[0];
      total_int_r[1] += int_r[1];
      total_length += length;
    } // end of loop over elements

    // Assemble the (\int r ds) and beam length to get the centre of mass for
    // one arm
    Vector<double> r_centre(2);
    r_centre[0] = (1.0 / total_length) * total_int_r[0];
    r_centre[1] = (1.0 / total_length) * total_int_r[1];

    // Compute the centre of mass of the entire beam
    sum_r_centre[0] = sum_r_centre[0] + r_centre[0];
    sum_r_centre[1] = sum_r_centre[1] + r_centre[1];
  }
}


//=============================================================================
/// Compute the drag and torque on the entire beam structure according to
/// slender body theory (Type=0: first arm, Type=1: second arm.)
//=============================================================================
void RigidBodyElement::compute_drag_and_torque(Vector<double>& sum_total_drag,
                                               double& sum_total_torque)
{
#ifdef PARANOID
  if (total_drag.size() != 2)
  {
    std::ostringstream error_message;
    error_message << "total_drag should have size 2, not " << total_drag.size()
                  << std::endl;

    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Initialise
  sum_total_drag[0] = 0.0;
  sum_total_drag[1] = 0.0;
  sum_total_torque = 0.0;

  Vector<double> drag(2);
  double torque = 0.0;

  // Find number of beam meshes
  unsigned npointer = Beam_mesh_pt.size();

  // Loop over the beam meshes to compute the drag and torque of the entire beam
  for (unsigned i = 0; i < npointer; i++)
  {
    // Initialise
    Vector<double> total_drag(2);
    double total_torque = 0.0;

    // Find number of elements in the mesh
    unsigned n_element = Beam_mesh_pt[i]->nelement();

    // Loop over the elements to compute the sum of elements' contribution to
    // the drag and torque on the entire beam structure
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast to the specific element type
      HaoHermiteBeamElement* elem_pt =
        dynamic_cast<HaoHermiteBeamElement*>(Beam_mesh_pt[i]->element_pt(e));

      // Compute contribution to the drag and torque within the e-th element
      elem_pt->compute_contribution_to_drag_and_torque(drag, torque);

      // Sum the elements' contribution to the drag and torque
      total_drag[0] += drag[0];
      total_drag[1] += drag[1];
      total_torque += torque;
    } // end of loop over elements

    // Compute the drag and torque of the entire beam
    sum_total_drag[0] = sum_total_drag[0] + total_drag[0];
    sum_total_drag[1] = sum_total_drag[1] + total_drag[1];
    sum_total_torque = sum_total_torque + total_torque;
  }
}


//======start_of_problem_class==========================================
/// Beam problem object
//======================================================================
class ElasticBeamProblem : public Problem
{
public:
  /// Constructor: The arguments are the number of elements and the parameter to
  /// determine the length of the beam
  ElasticBeamProblem(const unsigned& n_elem1, const unsigned& n_elem2);

  /// Conduct a parameter study
  // void parameter_study();

  /// No actions need to be performed after a solve
  void actions_after_newton_solve() {}

  /// No actions need to be performed before a solve
  void actions_before_newton_solve() {}

  /// Update the problem specs after solve (empty)
  void actions_after_implicit_timestep() {}

  /// Update the problem specs before next timestep:
  /// Set Dirchlet boundary conditions from exact solution.
  void actions_before_implicit_timestep() {}

  /// Set initial condition (incl previous timesteps) according
  /// to specified function.
  void set_initial_condition();

  /// Doc the solution
  void doc_solution(DocInfo& doc_info, ofstream& trace_file);

  /// Dump problem data to allow for later restart
  void dump_it(ofstream& dump_file)
  {
    // Current value of the FSI parameter
    dump_file << Global_Physical_Variables::I << " # FSI parameter"
              << std::endl;

    // hierher maybe add q and alpha but then issue warning
    // if it differs from the one specified on the command line

    // Dump the refinement pattern and the generic problem data
    Problem::dump(dump_file);
  }

  /// Read problem data for restart
  void restart(ifstream& restart_file)
  {
    // Read line up to termination sign
    string input_string;
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Read in FSI parameter
    Global_Physical_Variables::I = double(atof(input_string.c_str()));

    // Refine the mesh and read in the generic problem data
    Problem::read(restart_file);
  }

private:
  /// Pointer to geometric object that represents the beam's undeformed shape
  GeomObject* Undef_beam_pt1;

  GeomObject* Undef_beam_pt2;

  /// Pointer to RigidBodyElement that actually contains the rigid body data
  RigidBodyElement* Rigid_body_element_pt;

  /// Pointer to beam mesh (first arm)
  OneDLagrangianMesh<HaoHermiteBeamElement>* Beam_mesh_first_arm_pt;

  /// Pointer to beam mesh (second arm)
  OneDLagrangianMesh<HaoHermiteBeamElement>* Beam_mesh_second_arm_pt;

  /// Pointer to mesh containing the rigid body element
  Mesh* Rigid_body_element_mesh_pt;

}; // end of problem class


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=========================================================================
/// Steady, straight 1D line in 2D space
///  \f[ x = 0.0 \f]
///  \f[ y = \zeta \f]
//=========================================================================
class StraightLineVertical : public GeomObject
{
public:
  /// Constructor derives from GeomObject(1, 2)
  StraightLineVertical() : GeomObject(1, 2) {}

  /// Broken copy constructor
  StraightLineVertical(const StraightLineVertical& dummy) = delete;

  /// Broken assignment operator
  void operator=(const StraightLineVertical&) = delete;

  /// Position Vector at Lagrangian coordinate zeta
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    r[0] = 0.0;
    r[1] = zeta[0];
  }


  /// Derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// Evaluated at current time.
  virtual void dposition(const Vector<double>& zeta,
                         DenseMatrix<double>& drdzeta) const
  {
    // Tangent vector
    drdzeta(0, 0) = 0.0;
    drdzeta(0, 1) = 1.0;
  }


  /// 2nd derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i). Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }


  /// Posn Vector and its  1st & 2nd derivatives
  /// w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i).
  /// Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          Vector<double>& r,
                          DenseMatrix<double>& drdzeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Position Vector
    r[0] = 0.0;
    r[1] = zeta[0];

    // Tangent vector
    drdzeta(0, 0) = 0.0;
    drdzeta(0, 1) = 1.0;

    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }
};

//=========================================================================
/// Steady, straight 1D line in 2D space with stretch_ratio
///  \f[ x = 0.0 \f]
///  \f[ y = \zeta*stretch_ratio  \f]
//=========================================================================
class StraightLineVertical_new : public GeomObject
{
public:
  /// Constructor derives from GeomObject(1, 2)
  /// Constructor: Pass stretch_ratio
  StraightLineVertical_new(const double& stretch_ratio) : GeomObject(1, 2)
  {
    Stretch_ratio = stretch_ratio;
  }

  /// Broken copy constructor
  StraightLineVertical_new(const StraightLineVertical_new& dummy) = delete;

  /// Broken assignment operator
  void operator=(const StraightLineVertical_new&) = delete;

  /// Position Vector at Lagrangian coordinate zeta
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    r[0] = 0.0;
    r[1] = zeta[0] * Stretch_ratio;
  }


  /// Derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// Evaluated at current time.
  virtual void dposition(const Vector<double>& zeta,
                         DenseMatrix<double>& drdzeta) const
  {
    // Tangent vector
    drdzeta(0, 0) = 0.0;
    drdzeta(0, 1) = Stretch_ratio;
  }


  /// 2nd derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i). Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }


  /// Posn Vector and its  1st & 2nd derivatives
  /// w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i).
  /// Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          Vector<double>& r,
                          DenseMatrix<double>& drdzeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Position Vector
    r[0] = 0.0;
    r[1] = zeta[0] * Stretch_ratio;

    // Tangent vector
    drdzeta(0, 0) = 0.0;
    drdzeta(0, 1) = Stretch_ratio;

    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }

private:
  /// Define the length of the beam
  double Stretch_ratio;
};


//=============start_of_constructor=====================================
/// Constructor for elastic beam problem
//======================================================================
ElasticBeamProblem::ElasticBeamProblem(const unsigned& n_elem1,
                                       const unsigned& n_elem2)
{
  // Allocate the timestepper -- this constructs the Problem's
  // time object with a sufficient amount of storage to store the
  // previous timsteps.
  add_time_stepper_pt(new BDF<1>);

  // x position of clamped point
  double x0 = 0.0;

  // y position of clamped point
  double y0 = 0.0;

  // Beam's inclination
  double theta_eq = Global_Physical_Variables::Initial_value_for_theta_eq;

  // Make the RigidBodyElement that stores the parameters for the rigid body
  // motion
  Rigid_body_element_pt =
    new RigidBodyElement(x0, y0, theta_eq, time_stepper_pt());

  // Add the rigid body element to its own mesh
  Rigid_body_element_mesh_pt = new Mesh;
  Rigid_body_element_mesh_pt->add_element_pt(Rigid_body_element_pt);

  // Set the undeformed beam shape for two arms (in the reference orientation
  // before applying the rigid body motion)
  // Undef_beam_pt = new StraightLineVertical();

  // Aspect ratio to determine the length of the beam
  // first arm length = |q+0.5|, second arm length = |q-0.5|
  // double* q_pt = &Global_Physical_Variables::Q;

  // Still use the same expression of q as before to represent the length of the
  // two arms
  double* stretch_ratio_pt = &Global_Physical_Variables::Stretch_ratio;
  Undef_beam_pt1 = new StraightLineVertical_new(fabs(*stretch_ratio_pt + 0.5));
  Undef_beam_pt2 = new StraightLineVertical_new(fabs(*stretch_ratio_pt - 0.5));

  // Create the (Lagrangian!) mesh, using the StraightLineVertical object
  // Undef_beam_pt to specify the initial (Eulerian) position of the
  // nodes. (first arm)
  // double length_1 = fabs(*q_pt + 0.5);
  double length_1 = 1.0;
  Beam_mesh_first_arm_pt = new OneDLagrangianMesh<HaoHermiteBeamElement>(
    n_elem1, length_1, Undef_beam_pt1, time_stepper_pt());

  // Create the (Lagrangian!) mesh, using the StraightLineVertical object
  // Undef_beam_pt to specify the initial (Eulerian) position of the
  // nodes. (second arm)
  // double length_2 = fabs(*q_pt - 0.5);
  double length_2 = 1.0;
  Beam_mesh_second_arm_pt = new OneDLagrangianMesh<HaoHermiteBeamElement>(
    n_elem2, length_2, Undef_beam_pt2, time_stepper_pt());

  // Pass the pointer of the mesh to the RigidBodyElement class
  // so it can work out the drag and torque on the entire structure
  Vector<SolidMesh*> Beam_mesh_pt(2);
  Beam_mesh_pt[0] = Beam_mesh_first_arm_pt;
  Beam_mesh_pt[1] = Beam_mesh_second_arm_pt;
  Rigid_body_element_pt->set_pointer_to_beam_meshes(Beam_mesh_pt);

  // Build the problem's global mesh
  add_sub_mesh(Beam_mesh_first_arm_pt);
  add_sub_mesh(Beam_mesh_second_arm_pt);
  add_sub_mesh(Rigid_body_element_mesh_pt);
  build_global_mesh();

  // Set the boundary conditions: One end of the beam is clamped in space
  // Pin displacements in both x and y directions, and pin the derivative of
  // position Vector w.r.t. to coordinates in x direction. (first arm)
  Beam_mesh_first_arm_pt->boundary_node_pt(0, 0)->pin_position(0);
  Beam_mesh_first_arm_pt->boundary_node_pt(0, 0)->pin_position(1);
  Beam_mesh_first_arm_pt->boundary_node_pt(0, 0)->pin_position(1, 0);

  // Find number of elements in the mesh (first arm)
  unsigned n_element = Beam_mesh_first_arm_pt->nelement();

  // Loop over the elements to set physical parameters etc. (first arm)
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast to the specific element type
    HaoHermiteBeamElement* elem_pt = dynamic_cast<HaoHermiteBeamElement*>(
      Beam_mesh_first_arm_pt->element_pt(e));

    // Pass the pointer of RigidBodyElement to the each element
    // so we can work out the rigid body motion
    elem_pt->set_pointer_to_rigid_body_element(Rigid_body_element_pt);

    // Set physical parameters for each element:
    elem_pt->h_pt() = &Global_Physical_Variables::H;
    elem_pt->i_pt() = &Global_Physical_Variables::I;

    // Note: no rotation!

    // Set the undeformed shape for each element
    elem_pt->undeformed_beam_pt() = Undef_beam_pt1;

  } // end of loop over elements


  // Set the boundary conditions: One end of the beam is clamped in space
  // Pin displacements in both x and y directions, and pin the derivative of
  // position Vector w.r.t. to coordinates in x direction. (second arm)
  Beam_mesh_second_arm_pt->boundary_node_pt(0, 0)->pin_position(0);
  Beam_mesh_second_arm_pt->boundary_node_pt(0, 0)->pin_position(1);
  Beam_mesh_second_arm_pt->boundary_node_pt(0, 0)->pin_position(1, 0);

  // Find number of elements in the mesh (second arm)
  n_element = Beam_mesh_second_arm_pt->nelement();

  // Loop over the elements to set physical parameters etc. (second arm)
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast to the specific element type
    HaoHermiteBeamElement* elem_pt = dynamic_cast<HaoHermiteBeamElement*>(
      Beam_mesh_second_arm_pt->element_pt(e));

    // Pass the pointer of RigidBodyElement to the each element
    // so we can work out the rigid body motion
    elem_pt->set_pointer_to_rigid_body_element(Rigid_body_element_pt);

    // Set physical parameters for each element:
    elem_pt->h_pt() = &Global_Physical_Variables::H;
    elem_pt->i_pt() = &Global_Physical_Variables::I;

    // Rotate by opening angle
    elem_pt->theta_initial_pt(&Global_Physical_Variables::Alpha);

    // Set the undeformed shape for each element
    elem_pt->undeformed_beam_pt() = Undef_beam_pt2;

  } // end of loop over elements

  // Assign the global and local equation numbers
  cout << "# of dofs " << assign_eqn_numbers() << std::endl;

} // end of constructor

//=====start_of_set_ic=====================================================
/// Setup initial conditions -- either restart from solution
/// specified via command line or impulsive start.
//=========================================================================
void ElasticBeamProblem::set_initial_condition()
{
  // Assign impulsive start
  assign_initial_values_impulsive();
}


//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
void ElasticBeamProblem::doc_solution(DocInfo& doc_info, ofstream& trace_file)
{
  ofstream some_file1;
  char filename[100];

  // Number of plot points
  unsigned npts;
  npts = 5;


  cout << std::endl;
  cout << "=================================================" << std::endl;
  cout << "Docing solution for t=" << time_pt()->time() << std::endl;
  cout << "=================================================" << std::endl;


  // Output solution
  //-----------------
  sprintf(filename,
          "%s/soln_first_arm_%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file1.open(filename);
  Beam_mesh_first_arm_pt->output(some_file1, npts);
  some_file1.close();

  ofstream some_file2;
  sprintf(filename,
          "%s/soln_second_arm_%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file2.open(filename);
  Beam_mesh_second_arm_pt->output(some_file2, npts);
  some_file2.close();

  trace_file << Global_Physical_Variables::I << "  ";
  Rigid_body_element_pt->output(trace_file);
  trace_file << time_pt()->time() << std::endl;

} // end of doc_solution


/* //=======start_of_parameter_study==========================================
/// Solver loop to perform parameter study
//=========================================================================
void ElasticBeamProblem::parameter_study()
{
  // Over-ride the default maximum value for the residuals
  // Problem::Max_residuals = 1.0e10;
  // Problem::Max_newton_iterations = 20;
  Problem::Always_take_one_newton_step = true;
  Problem::Scale_arc_length = false;
  Problem::Theta_squared = 0.2;

  // Create label for output
  DocInfo doc_info;

  // Set output directory -- this function checks if the output
  // directory exists and issues a warning if it doesn't.
  doc_info.set_directory("RESLT");

  // Output file stream used for writing results
  ofstream file;

  // String used for the filename
  char filename[100];

  // Write the file name
  sprintf(filename,
          "RESLT/elastic_beam_I_theta_q_%.3f_alpha_%.3fpi_initial_%.2f.dat",
          Global_Physical_Variables::Stretch_ratio,
          Global_Physical_Variables::Alpha / acos(-1.0),
          Global_Physical_Variables::Initial_value_for_theta_eq);
  file.open(filename);

  // Counter to record the iterations for the while loop
  unsigned counter = 0;

  // Initialize the value of backup for dofs
  DoubleVector dofs_backup;

  // ////////////////////////////////////////////////////////////////////////
  // // TEST!
  // ofstream file1;
  // ofstream file2;
  // // Document the solution (first arm)
  // sprintf(filename,
  //         "RESLT/beam_first_arm_initial_%.2f_%d.dat",
  //         Global_Physical_Variables::Initial_value_for_theta_eq,
  //         counter);
  // file1.open(filename);
  // Beam_mesh_first_arm_pt->output(file1, 5);
  // file1.close();

  // // Document the solution (second arm)
  // sprintf(filename,
  //         "RESLT/beam_second_arm_initial_%.2f_%d.dat",
  //         Global_Physical_Variables::Initial_value_for_theta_eq,
  //         counter);
  // file2.open(filename);
  // Beam_mesh_second_arm_pt->output(file2, 5);
  // file2.close();
  // ///////////////////////////////////////////////////////////////////////////


  // exit(0);

  // Loop over different values for Non-dimensional coefficient (FSI) I by
  // using arclength increment
  // The loop stops when I becomes negative since only positive values of I are
  // considered.
  double I_backup = 0.0;
  double ds = 0.0;
  while (Global_Physical_Variables::I >= 0.0)
  {
    // Get the dofs
    Problem::get_dofs(dofs_backup);

    try
    {
      if (counter == 0)
      {
        // Solve the system
        newton_solve();
      }
      else
      {
        // Backup the FSI coefficient I
        I_backup = Global_Physical_Variables::I;

        // To prevent large solution jumps in critical intervals for I, try to
        // reduce ds specifically in those areas for smoother and more precise
        // results

        // First I interval [Interval1_start,Interval1_end]
        if (Global_Physical_Variables::Ds_default >
              Global_Physical_Variables::Ds_interval1 &&
            Global_Physical_Variables::I >=
              Global_Physical_Variables::Interval1_start &&
            Global_Physical_Variables::I <=
              Global_Physical_Variables::Interval1_end)
        {
          ds = Global_Physical_Variables::Ds_interval1;
        }
        // Second I interval [Interval2_start,Interval2_end]
        else if (Global_Physical_Variables::Ds_default >
                   Global_Physical_Variables::Ds_interval2 &&
                 Global_Physical_Variables::I >=
                   Global_Physical_Variables::Interval2_start &&
                 Global_Physical_Variables::I <=
                   Global_Physical_Variables::Interval2_end)
        {
          ds = Global_Physical_Variables::Ds_interval2;
        }
        else
        {
          // Use the default one
          ds = Global_Physical_Variables::Ds_default;
        }

        /// Use the arclength solve
        ds = arc_length_step_solve(&Global_Physical_Variables::I, ds);
      }

      // Document I
      file << Global_Physical_Variables::I << "  ";

      // Document the solution of Theta_eq, Theta_eq_orientation
      Rigid_body_element_pt->output(file);

      // Document maximum residuals at start and after each newton iteration
      file << Problem::Max_res[0] << "  ";

      // Document actual number of Newton iterations taken during the most
      // recent iteration
      file << Problem::Nnewton_iter_taken << "  ";

      // Step label
      file << counter << std::endl;

      // Output file stream used for writing results
      ofstream file1;
      ofstream file2;
      ofstream file3;

      // Document the solution (first arm)
      sprintf(filename,
              "RESLT/beam_first_arm_initial_%.2f_%d.dat",
              Global_Physical_Variables::Initial_value_for_theta_eq,
              counter);
      file1.open(filename);
      Beam_mesh_first_arm_pt->output(file1, 5);
      file1.close();

      // Document the solution (second arm)
      sprintf(filename,
              "RESLT/beam_second_arm_initial_%.2f_%d.dat",
              Global_Physical_Variables::Initial_value_for_theta_eq,
              counter);
      file2.open(filename);
      Beam_mesh_second_arm_pt->output(file2, 5);
      file2.close();

      // Write restart file
      sprintf(filename, "RESLT/restart%i.dat", counter);
      file2.open(filename);
      dump_it(file2);
      file2.close();

      // Bump counter for output
      counter = counter + 1;
    }
    catch (...)
    {
      // If the initial values are not appropriate, the newton method cannot
      // converge at the starting point I=0.0
      if (counter == 0)
      {
        oomph_info
          << "Initial values for I=0.0 are not appropriate. Please modify them!"
          << std::endl;
        break;
      }
      else
      {
        // Check if default ds is small enough
        if (Global_Physical_Variables::Ds_default < 1.0e-8)
        {
          // Since the default is not small enough, assign a smaller value
          Global_Physical_Variables::Ds_default = 1.0e-8;

          // Since the newton method cannot converge in this turn, overwrite the
          // failed I value
          Global_Physical_Variables::I = I_backup;
        }
        else
        {
          // If default ds is already small enough, it means there is no
          // solution any more
          break;
        }
      }

      // Reset the dofs
      Problem::set_dofs(dofs_backup);
    }
  }
  file.close();


} // end of parameter study */

//========start_of_main================================================
/// Driver for beam (string under tension) test problem
//=====================================================================
int main(int argc, char** argv)
{
  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Stretch_ratio
  CommandLineArgs::specify_command_line_flag(
    "--q", &Global_Physical_Variables::Stretch_ratio);

  // Aspect ratio
  // CommandLineArgs::specify_command_line_flag("--q",
  //&Global_Physical_Variables::Q);

  // Opening angle in degrees
  double alpha_in_degrees = 45.0;
  CommandLineArgs::specify_command_line_flag("--alpha_in_degrees",
                                             &alpha_in_degrees);

  // Initial value for theta_eq in the Newton solve
  CommandLineArgs::specify_command_line_flag(
    "--Initial_value_for_theta_eq",
    &Global_Physical_Variables::Initial_value_for_theta_eq);

  // Initial value for theta_eq in the Newton solve
  CommandLineArgs::specify_command_line_flag(
    "--ds_default", &Global_Physical_Variables::Ds_default);

  // End point for the first I interval [Interval1_start,Interval1_end]
  CommandLineArgs::specify_command_line_flag(
    "--interval1_start", &Global_Physical_Variables::Interval1_start);

  CommandLineArgs::specify_command_line_flag(
    "--interval1_end", &Global_Physical_Variables::Interval1_end);

  // Start and end points for the second I interval
  // [Interval2_start,Interval2_end]
  CommandLineArgs::specify_command_line_flag(
    "--interval2_start", &Global_Physical_Variables::Interval2_start);

  CommandLineArgs::specify_command_line_flag(
    "--interval2_end", &Global_Physical_Variables::Interval2_end);

  // Value of ds for first interval
  CommandLineArgs::specify_command_line_flag(
    "--ds_interval1", &Global_Physical_Variables::Ds_interval1);

  // Value of ds for second interval
  CommandLineArgs::specify_command_line_flag(
    "--ds_interval2", &Global_Physical_Variables::Ds_interval2);

  // Restart file
  std::string restart_file;
  CommandLineArgs::specify_command_line_flag("--restart_file", &restart_file);

  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  // Now that we've read the opening angle in degrees, update the value in
  // radians
  Global_Physical_Variables::Alpha = 4.0 * atan(1.0) / 180.0 * alpha_in_degrees;

  // Set the non-dimensional thickness
  Global_Physical_Variables::H = 0.01;

  // Number of elements (choose an even number if you want the control point
  // to be located at the centre of the beam)
  unsigned n_element1 = 20;
  unsigned n_element2 = 20;

  Global_Physical_Variables::I = 0.0001;

  // Construct the problem
  ElasticBeamProblem problem(n_element1, n_element2);

  // Setup labels for output
  DocInfo doc_info;

  // Output directory
  doc_info.set_directory("RESLT");

  // Output number
  doc_info.number() = 0;

  // Open a trace file
  ofstream trace_file;
  char filename[100];
  sprintf(filename, "%s/trace.dat", doc_info.directory().c_str());
  trace_file.open(filename);

  // Choose simulation interval and timestep
  double t_max = 500.0;
  double dt = 1;

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  problem.initialise_dt(dt);

  // Set IC
  problem.set_initial_condition();

  // Output initial condition
  problem.doc_solution(doc_info, trace_file);

  // Increment counter for solutions
  doc_info.number()++;

  // Find number of steps
  unsigned nstep = unsigned(t_max / dt);

  // Timestepping loop
  for (unsigned istep = 0; istep < nstep; istep++)
  {
    cout << " Timestep " << istep << std::endl;

    // Take timestep
    problem.unsteady_newton_solve(dt);

    // Output solution
    problem.doc_solution(doc_info, trace_file);

    // Increment counter for solutions
    doc_info.number()++;
  }

  // Close trace file
  trace_file.close();

  // Do the restart?
  if (CommandLineArgs::command_line_flag_has_been_set("--restart_file"))
  {
    // Open/read restart file
    std::ifstream file2;
    file2.open(restart_file.c_str());
    problem.restart(file2);
    file2.close();
  }

  // Conduct parameter study
  // problem.parameter_study();

} // end of main
