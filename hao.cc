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
// Driver function for a simple beam proble

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

  /// 2nd Piola Kirchhoff pre-stress
  double Sigma0 = 0.0;

  /// Pressure load
  double P_ext = 0.0;

  /// Non-dimensional coefficeient
  double Scale = 0.0;

  /// Shear rate
  double Gamma_dot = 0.0;

  /// Initial drift speed and accerlation of horiztontal motion
  double V = 0.0;

  /// Initial speed of horiztontal motion
  double U0 = 0.0;

  /// Initial Beam's orientation
  double Theta_eq = -0.3;

  /// Initial x position of clamped point
  double X0 = 0.8;

  /// Initial y position of clamped point
  double Y0 = 0.5;

  /// Load function: Apply a constant external pressure to the beam
  void load(const Vector<double>& xi,
            const Vector<double>& x,
            const Vector<double>& N,
            Vector<double>& load)
  {
    for (unsigned i = 0; i < 2; i++)
    {
      load[i] = -P_ext * N[i];
    }
    // load[0] = Scale*Gamma_dot*(x[1]-0.5*N[1]*N[1]*x[1]);
    // load[1] = Scale*Gamma_dot*0.5*N[0]*N[1]*x[1];
    // load[0] = 1.0e-7;
    // load[1] = 1.0e-7;
  }

} // namespace Global_Physical_Variables


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=========================================================================
/// hierher; come up with better name
//=========================================================================
class HaoElement : public GeneralisedElement
{
public:
  /// Constructor:  Pass initial values for rigid body parameters (pinned
  /// by default)
  HaoElement(const double& V,
             const double& U0,
             const double& Theta_eq,
             const double& X0,
             const double& Y0)
  {
    // Create internal data which contains the "rigid body" parameters
    for (unsigned i = 0; i < 5; i++)
    {
      // Create data: One value, no timedependence, free by default
      add_internal_data(new Data(1));

      // Pin the data
      internal_data_pt(i)->pin(0);
    }

    // Give them a value:
    internal_data_pt(0)->set_value(0, V);
    internal_data_pt(1)->set_value(0, U0);
    internal_data_pt(2)->set_value(0, Theta_eq);
    internal_data_pt(3)->set_value(0, X0);
    internal_data_pt(4)->set_value(0, Y0);
  }


  /// Function that returns the Vector of pointers to the "rigid body"
  /// parameters
  Vector<Data*> rigid_body_parameters()
  {
    Vector<Data*> tmp_pt(5);
    for (unsigned i = 0; i < 5; i++)
    {
      tmp_pt[i] = internal_data_pt(i);
    }
    return tmp_pt;
  }


  /// Helper function to compute the meaningful parameter values
  /// from enumerated data
  void get_parameters(
    double& V, double& U0, double& Theta_eq, double& X0, double& Y0)
  {
    V = internal_data_pt(0)->value(0);
    U0 = internal_data_pt(1)->value(0);
    Theta_eq = internal_data_pt(2)->value(0);
    X0 = internal_data_pt(3)->value(0);
    Y0 = internal_data_pt(4)->value(0);
  }


  /// Pass pointer to the Mesh of HaoHermiteBeamElements
  void set_pointer_to_mesh(Mesh* mesh_pt)
  {
    // Store the pointer for future reference
    Beam_mesh_pt = mesh_pt;
  }


  /// hierher implement this: Compute the sum of the elements' contribution
  /// to the drag and torque on the entire beam structure according to
  /// slender body theory
  void compute_drag_and_torque(Vector<double>& sum_drag, double& sum_torque);
  //{
  //  HIERHER DO THIS LOOP OVER ELEMENTS; THIS NEEDS THE POINTER
  //  TO THE BEAM MESH SO WRITE A SET FUNCTION.


  /// Output the sum of drag and torque on the entire beam structure
  void output(std::ostream& outfile)
  {
    Vector<double> sum_drag(2);
    double sum_torque = 0;

    // Compute the sum of drag and torque on the entire beam structure
    compute_drag_and_torque(sum_drag, sum_torque);

    // Output the sum of drag and torque
    outfile << sum_drag[0] << "  ";
    outfile << sum_drag[1] << "  ";
    outfile << sum_torque << std::endl;
  }

private:
  /// Pointer to the Mesh of HaoHermiteBeamElements
  Mesh* Beam_mesh_pt;
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
  /// Pass pointer to HaoElement that contains the rigid body parameters
  void set_hao_element_hierher_rename_me(HaoElement* hao_element_pt)
  {
    // Store the pointer for future reference
    Hao_element_pt = hao_element_pt;

    // Pass the rigid body parameters to HaoHermiteBeamElement
    Vector<Data*> rigid_body_data_pt = Hao_element_pt->rigid_body_parameters();

#ifdef PARANOID
    if (rigid_body_data_pt.size() != 5)
    {
      std::ostringstream error_message;
      error_message << "rigid_body_data_pt should have size 5, not "
                    << rigid_body_data_pt.size() << std::endl;

      // hierher loop over all entries
      for (unsigned i = 0; i < 5; i++)
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

    // add the rigid body parameters as the external data for each element
    for (unsigned i = 0; i < 5; i++)
    {
      add_external_data(rigid_body_data_pt[i]);
    }
  }


  /// Compute the local traction vector acting on the element
  /// at the local coordinate s according to slender body theory
  void compute_local_slender_body_traction(const Vector<double>& s,
                                           Vector<double>& drag)
  {
    abort();
    // drag[0]=V-0.5;
    // drag[1]=U0-0.3;
    // torque=Theta_eq-(0.3);
  }


  /// Compute the slender body traction onto the element at
  /// local coordinate s
  void compute_slender_body_traction(const Vector<double>& s,
                                     Vector<double>& traction)
  {
    // Get the Eulerian position and the unit normal
    Vector<double> posn(2);
    Vector<double> N(2);
    get_normal(s, posn, N);

    // hierher don't refer to namespace directly (but fine for now...)
    traction[0] = Global_Physical_Variables::Scale *
                  Global_Physical_Variables::Gamma_dot *
                  (posn[1] - 0.5 * N[1] * N[1] * posn[1]);

    traction[1] = Global_Physical_Variables::Scale *
                  Global_Physical_Variables::Gamma_dot * 0.5 * N[0] * N[1] *
                  posn[1];

    // traction[0] = 1.0e-7;
    // traction[1] = 1.0e-7;
  }


  /// Compute the element's contribution to the drag and torque on
  /// the entire beam structure according to slender body theory
  void compute_integrated_drag_and_torque(Vector<double>& drag, double& torque)
  {
// hierher paranoid check size of drag vector!
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

    // Local coordinate (1D!)
    Vector<double> s(1);

    // Set # of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      /// Return local coordinate s[j]  of i-th integration point.
      unsigned j = 0;
      s[j] = integral_pt()->knot(ipt, j);


      /// Get position vector to and non-unit tangent vector on wall:
      /// dr/ds
      Vector<double> posn(2);
      Vector<double> drds(2);
      get_non_unit_tangent(s, posn, drds);

      // Jacobian
      double J = sqrt(drds[0] * drds[0] + drds[1] * drds[1]);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Compute the slender body traction
      Vector<double> traction(2);
      compute_slender_body_traction(s, traction);

      // Translate parameters into meaningful variables do this elsewhere too
      double V = 0.0;
      double U0 = 0.0;
      double Theta_eq = 0.0;
      double X0 = 0.0;
      double Y0 = 0.0;
      Hao_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

      // calculate the contribution to torque
      double local_torque = 0.0;
      local_torque =
        (posn[0] - X0) * traction[1] - (posn[1] - Y0) * traction[0];

      // Add 'em
      drag[0] += traction[0] * W;
      drag[1] += traction[1] * W;
      torque += local_torque * W;
    }
  }


  /// Overloaded output function
  void output(std::ostream& outfile, const unsigned& n_plot)
  {
    // oomph_info << "Hello world --  this should work!" << std::endl;


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

    Vector<double> posn(n_dim);

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
        // Initialise acclerations and veloc
        posn[i] = 0.0;
      }

      // Calculate displacements, accelerations and spatial derivatives
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
            posn[i] += raw_dnodal_position_gen_dt(0, l, k, i) * psi(l, k);
          }
        }
      }

      // Get the normal vector at each plotted point
      Vector<double> N(n_dim);
      get_normal(s, N);

      double scalar_N = 0.0;
      Vector<double> load(n_dim);

      // Compute slender body traction
      compute_slender_body_traction(s, load);

      // Output position
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << posn[i] << " ";
      }

      // Output unit normal
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << N[i] << " ";
        // Verify the norm of the unit normal vector
        scalar_N += pow(N[i], 2);
      }
      outfile << sqrt(scalar_N) << " ";

      // Output U_inf
      outfile << Global_Physical_Variables::Gamma_dot * posn[1] << " ";
      outfile << 0 << " ";

      // Output traction
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << load[i] << " ";
      }

      // Test. Try to move beam back to the original position
      outfile << (posn[0] - Global_Physical_Variables::X0) *
                     cos(Global_Physical_Variables::Theta_eq) +
                   (posn[1] - Global_Physical_Variables::Y0) *
                     sin(Global_Physical_Variables::Theta_eq)
              << "  ";
      outfile << -(posn[0] - Global_Physical_Variables::X0) *
                     sin(Global_Physical_Variables::Theta_eq) +
                   (posn[1] - Global_Physical_Variables::Y0) *
                     cos(Global_Physical_Variables::Theta_eq);
      // outfile <<
      // posn[0]*cos(-Global_Physical_Variables::Theta_eq)-posn[1]*sin(-Global_Physical_Variables::Theta_eq)-Global_Physical_Variables::X0<<"
      // "; outfile <<
      // posn[0]*sin(-Global_Physical_Variables::Theta_eq)+posn[1]*cos(-Global_Physical_Variables::Theta_eq)-Global_Physical_Variables::Y0;
      outfile << std::endl;
    }
  }

private:
  /// Pointer to element that controls the rigid body motion
  HaoElement* Hao_element_pt;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=============================================================================
/// Compute the sum of the elements' contribution to the drag and torque on the
/// entire beam structure according to slender body theory
//=============================================================================
void HaoElement::compute_drag_and_torque(Vector<double>& sum_drag,
                                         double& sum_torque)
{
#ifdef PARANOID
  if (sum_drag.size() != 2)
  {
    std::ostringstream error_message;
    error_message << "sum_drag should have size 2, not " << sum_drag.size()
                  << std::endl;

    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Find number of elements in the mesh
  unsigned n_element = Beam_mesh_pt->nelement();

  Vector<double> drag(2);
  double torque = 0;

  // Loop over the elements to compute the sum of elements' contribution to the
  // drag and torque on the entire beam structure
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast to the specific element type
    HaoHermiteBeamElement* elem_pt =
      dynamic_cast<HaoHermiteBeamElement*>(Beam_mesh_pt->element_pt(e));

    // Compute contribution to the drag and torque within the e-th element
    elem_pt->compute_integrated_drag_and_torque(drag, torque);

    // Sum the elements' contribution to the drag and torque
    sum_drag[0] += drag[0];
    sum_drag[1] += drag[1];
    sum_torque += torque;
  } // end of loop over elements
}


//======start_of_problem_class==========================================
/// Beam problem object
//======================================================================
class ElasticBeamProblem : public Problem
{
public:
  /// Constructor: The arguments are the number of elements,
  /// the length of domain
  ElasticBeamProblem(const unsigned& n_elem, const double& length);

  /// Conduct a parameter study
  void parameter_study();

  /// Return pointer to the mesh
  OneDLagrangianMesh<HaoHermiteBeamElement>* mesh_pt()
  {
    return dynamic_cast<OneDLagrangianMesh<HaoHermiteBeamElement>*>(
      Problem::mesh_pt());
  }

  /// No actions need to be performed after a solve
  void actions_after_newton_solve() {}

  /// No actions need to be performed before a solve
  void actions_before_newton_solve() {}

private:
  /// Pointer to the node whose displacement is documented
  Node* Doc_node_pt;

  /// Length of domain (in terms of the Lagrangian coordinates)
  double Length;

  /// Pointer to geometric object that represents the beam's undeformed shape
  GeomObject* Undef_beam_pt;

  /// Pointer to HaoElement that actually contains the rigid body data
  HaoElement* Hao_element_pt;

}; // end of problem class


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=========================================================================
/// Steady, straight 1D line in 2D space
///  \f[ x = \zeta \f]  hierher clean up
///  \f[ y = H \f]
//=========================================================================
class StraightLineVertical : public GeomObject
{
public:
  /// Constructor: Pass in pointer to HaoElement that contains the
  /// rigid body Data
  StraightLineVertical(HaoElement* hao_element_pt)
    : GeomObject(1, 2), Hao_element_pt(hao_element_pt)
  {
    // Pass the rigid body parameters to this class
    Vector<Data*> rigid_body_data_pt = Hao_element_pt->rigid_body_parameters();

#ifdef PARANOID
    if (rigid_body_data_pt.size() != 5)
    {
      std::ostringstream error_message;
      error_message << "rigid_body_data_pt should have size 5, not "
                    << rigid_body_data_pt.size() << std::endl;

      // hierher loop over all entries
      for (unsigned i = 0; i < 5; i++)
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

    // V, U0, X0, Y0, theta_eq they're all members of this data object
    Geom_data_pt.resize(5);
    for (unsigned i = 0; i < 5; i++)
    {
      Geom_data_pt[i] = rigid_body_data_pt[i];
    }

    // Data has been created externally: Must not clean up
    Must_clean_up = false;
  }

  /// Constructor:  Pass parameters (pinned by default)
  StraightLineVertical(const double& V,
                       const double& U0,
                       const double& Theta_eq,
                       const double& X0,
                       const double& Y0)
    : GeomObject(1, 2)
  {
    // Create Data for straightVertical-line object
    Geom_data_pt.resize(5);
    for (unsigned i = 0; i < 5; i++)
    {
      // Create data: One value, no timedependence, free by default
      Geom_data_pt[i] = new Data(1);

      // Pin the data
      Geom_data_pt[i]->pin(0);
    }

    // I've created the data, I need to clean up
    Must_clean_up = true;

    // Give them values:
    Geom_data_pt[0]->set_value(0, V);
    Geom_data_pt[1]->set_value(0, U0);
    Geom_data_pt[2]->set_value(0, Theta_eq);
    Geom_data_pt[3]->set_value(0, X0);
    Geom_data_pt[4]->set_value(0, Y0);
  }


  /// Broken copy constructor
  StraightLineVertical(const StraightLineVertical& dummy) = delete;

  /// Broken assignment operator
  void operator=(const StraightLineVertical&) = delete;

  /// Destructor:  Clean up if necessary
  ~StraightLineVertical()
  {
    // Do I need to clean up?
    if (Must_clean_up)
    {
      for (unsigned i = 0; i < 5; i++)
      {
        delete Geom_data_pt[i];
        Geom_data_pt[i] = 0;
      }
    }
  }


  /// Position Vector at Lagrangian coordinate zeta
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    // assign a value to t
    double t = 0.0;

    // Translate parameters into meaningful variables do this elsewhere too
    double V = 0.0;
    double U0 = 0.0;
    double Theta_eq = 0.0;
    double X0 = 0.0;
    double Y0 = 0.0;
    Hao_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

    // Position Vector: for t=0 so it depends on X_0, Y_0 and theta_eq
    // later: we'll provide a properly time dependent version which takes
    // the effect of V and U0 into account. // hierher
    r[0] = -zeta[0] * sin(Theta_eq) + 0.5 * V * t * t + U0 * t + X0;
    r[1] = zeta[0] * cos(Theta_eq) + V * t + Y0;
  }


  /// Parametrised position on object: r(zeta). Evaluated at
  /// previous timestep. t=0: current time; t>0: previous
  /// timestep.
  void position(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& r) const
  {
#ifdef PARANOID
    if (t > Geom_data_pt[0]->time_stepper_pt()->nprev_values())
    {
      std::ostringstream error_message;
      error_message << "t > nprev_values() " << t << " "
                    << Geom_data_pt[0]->time_stepper_pt()->nprev_values()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif


    // hierher wrong -- t is the discrete previous time level.
    oomph_info << "wrong -- t is the discrete previous time level;"
               << "implement properly and/or provide a version that"
               << "depends on continuous time." << std::endl;
    abort();

    // // Position Vector at time level t
    // r[0] =
    // -zeta[0]*sin(Geom_data_pt[2]->value(0))+0.5*(Geom_data_pt[0]->value(0))*t*t+(Geom_data_pt[1]->value(0))*t+Geom_data_pt[3]->value(0);
    // r[1] =
    // zeta[0]*cos(Geom_data_pt[2]->value(0))+(Geom_data_pt[0]->value(0))*t+Geom_data_pt[4]->value(0);
  }


  /// Derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// Evaluated at current time.
  virtual void dposition(const Vector<double>& zeta,
                         DenseMatrix<double>& drdzeta) const
  {
    // Translate parameters into meaningful variables do this elsewhere too
    double V = 0.0;
    double U0 = 0.0;
    double Theta_eq = 0.0;
    double X0 = 0.0;
    double Y0 = 0.0;
    Hao_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

    // Tangent vector
    drdzeta(0, 0) = -sin(Theta_eq);
    drdzeta(0, 1) = cos(Theta_eq);
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
    // assign a value to t
    double t = 0.0;

    // Translate parameters into meaningful variables do this elsewhere too
    double V = 0.0;
    double U0 = 0.0;
    double Theta_eq = 0.0;
    double X0 = 0.0;
    double Y0 = 0.0;
    Hao_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

    // Position Vector
    r[0] = -zeta[0] * sin(Theta_eq) + 0.5 * V * t * t + U0 * t + X0;
    r[1] = zeta[0] * cos(Theta_eq) + V * t + Y0;

    // Tangent vector
    drdzeta(0, 0) = -sin(Theta_eq);
    drdzeta(0, 1) = cos(Theta_eq);

    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }


  /// How many items of Data does the shape of the object depend on?
  unsigned ngeom_data() const
  {
    return Geom_data_pt.size();
  }

  /// Return pointer to the j-th Data item that the object's
  /// shape depends on
  Data* geom_data_pt(const unsigned& j)
  {
    return Geom_data_pt[j];
  }

private:
  /// Vector of pointers to Data items that affects the object's shape
  Vector<Data*> Geom_data_pt;

  /// Do I need to clean up?
  bool Must_clean_up;

  /// Pointer to HaoElement that actually contains all the data
  HaoElement* Hao_element_pt;
};


//=============start_of_constructor=====================================
/// Constructor for elastic beam problem
//======================================================================
ElasticBeamProblem::ElasticBeamProblem(const unsigned& n_elem,
                                       const double& length)
  : Length(length)
{
  // // Set the undeformed beam shape
  // Undef_beam_pt=new StraightLineVertical(Global_Physical_Variables::V,
  //                                        Global_Physical_Variables::U0,
  //                                        Global_Physical_Variables::Theta_eq,
  //                                        Global_Physical_Variables::X0,
  //                                        Global_Physical_Variables::Y0);


  // Make the HaoElement that stores the parameters for the rigid body motion
  Hao_element_pt = new HaoElement(Global_Physical_Variables::V,
                                  Global_Physical_Variables::U0,
                                  Global_Physical_Variables::Theta_eq,
                                  Global_Physical_Variables::X0,
                                  Global_Physical_Variables::Y0);


  // Set the undeformed beam shape
  Undef_beam_pt = new StraightLineVertical(Hao_element_pt);

  // Undef_beam_pt=new StraightLineVertical(0.0,0.0,0.0,0.0,0.0);

  // Create the (Lagrangian!) mesh, using the geometric object
  // Undef_beam_pt to specify the initial (Eulerian) position of the
  // nodes.
  Problem::mesh_pt() = new OneDLagrangianMesh<HaoHermiteBeamElement>(
    n_elem, length, Undef_beam_pt);

  // Pass the pointer of the mesh to the HaoElement class
  Hao_element_pt->set_pointer_to_mesh(mesh_pt());

  // Set the boundary conditions: Each end of the beam is fixed in space
  // Loop over the boundaries (ends of the beam)
  for (unsigned b = 0; b < 1; b++)
  {
    // Pin displacements in both x and y directions, and pin the derivative of
    // position Vector w.r.t. to coordinates in x direction. [Note: The
    // mesh_pt() function has been overloaded
    //  to return a pointer to the actual mesh, rather than
    //  a pointer to the Mesh base class. The current mesh is derived
    //  from the SolidMesh class. In such meshes, all access functions
    //  to the nodes, such as boundary_node_pt(...), are overloaded
    //  to return pointers to SolidNodes (whose position can be
    //  pinned) rather than "normal" Nodes.]
    mesh_pt()->boundary_node_pt(b, 0)->pin_position(0);
    mesh_pt()->boundary_node_pt(b, 0)->pin_position(1);
    mesh_pt()->boundary_node_pt(b, 0)->pin_position(1, 0);
  }

  // Find number of elements in the mesh
  unsigned n_element = mesh_pt()->nelement();

  // Loop over the elements to set physical parameters etc.
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast to the specific element type
    HaoHermiteBeamElement* elem_pt =
      dynamic_cast<HaoHermiteBeamElement*>(mesh_pt()->element_pt(e));

    // Pass the pointer of HaoElement to the each element
    elem_pt->set_hao_element_hierher_rename_me(Hao_element_pt);

    // Set physical parameters for each element:
    elem_pt->sigma0_pt() = &Global_Physical_Variables::Sigma0;
    elem_pt->h_pt() = &Global_Physical_Variables::H;

    // Set the load Vector for each element
    elem_pt->load_vector_fct_pt() = &Global_Physical_Variables::load;

    // Set the undeformed shape for each element
    elem_pt->undeformed_beam_pt() = Undef_beam_pt;
  } // end of loop over elements

  // Choose node at which displacement is documented (halfway along -- provided
  // we have an odd number of nodes; complain if this is not the
  // case because the comparison with the exact solution will be wrong
  // otherwise!)
  unsigned n_nod = mesh_pt()->nnode();
  if (n_nod % 2 != 1)
  {
    cout << "Warning: Even number of nodes " << n_nod << std::endl;
    cout << "Comparison with exact solution will be misleading..." << std::endl;
  }
  Doc_node_pt = mesh_pt()->node_pt((n_nod + 1) / 2 - 1);

  // Assign the global and local equation numbers
  cout << "# of dofs " << assign_eqn_numbers() << std::endl;

} // end of constructor


//=======start_of_parameter_study==========================================
/// Solver loop to perform parameter study
//=========================================================================
void ElasticBeamProblem::parameter_study()
{
  // Over-ride the default maximum value for the residuals
  Problem::Max_residuals = 1.0e10;

  // Set the increments in control parameters
  // double pext_increment = 0.001;
  double pext_increment = 1.0e-6;

  // Set initial values for control parameters
  Global_Physical_Variables::P_ext = 0.0 - pext_increment;

  Global_Physical_Variables::Scale = 1.0e-4;
  // Global_Physical_Variables::Scale = 0.0;

  // Create label for output
  DocInfo doc_info;

  // Set output directory -- this function checks if the output
  // directory exists and issues a warning if it doesn't.
  doc_info.set_directory("RESLT");

  // Open a trace file
  ofstream trace("RESLT/trace_beam.dat");

  // Write a header for the trace file
  trace << "VARIABLES=\"p_e_x_t\",\"d\""
        << ", \"p_e_x_t_(_e_x_a_c_t_)\"" << std::endl;

  // Output file stream used for writing results
  ofstream file;
  // String used for the filename
  char filename[100];

  // Loop over parameter increments
  unsigned nstep = 10;
  for (unsigned i = 1; i <= nstep; i++)
  {
    // Increment pressure
    Global_Physical_Variables::P_ext += pext_increment;
    Global_Physical_Variables::Gamma_dot = (1.0e-2) * i;

    // Solve the system
    newton_solve();

    // Calculate exact solution for `string under tension' (applicable for
    // small wall thickness and pinned ends)

    // The tangent of the angle beta
    double tanbeta = -2.0 * Doc_node_pt->x(1) / Length;

    double exact_pressure = 0.0;
    // If the beam has deformed, calculate the pressure required
    if (tanbeta != 0)
    {
      // Calculate the opening angle alpha
      double alpha = 2.0 * atan(2.0 * tanbeta / (1.0 - tanbeta * tanbeta));

      // Jump back onto the main branch if alpha>180 degrees
      if (alpha < 0) alpha += 2.0 * MathematicalConstants::Pi;

      // Green strain:
      double gamma =
        0.5 *
        (0.25 * alpha * alpha / (sin(0.5 * alpha) * sin(0.5 * alpha)) - 1.0);

      // Calculate the exact pressure
      exact_pressure = Global_Physical_Variables::H *
                       (Global_Physical_Variables::Sigma0 + gamma) * alpha /
                       Length;
    }

    // Document the solution
    sprintf(filename, "RESLT/beam%i.dat", i);
    file.open(filename);
    mesh_pt()->output(file, 5);
    file.close();

    // Write trace file: Pressure, displacement and exact solution
    // (for string under tension)
    trace << Global_Physical_Variables::P_ext << " " << abs(Doc_node_pt->x(1))
          << " " << exact_pressure << std::endl;

    // Output file stream used for writing results
    ofstream file2;

    // String used for the filename
    char filename2[100];

    // Open the file and give it a name
    sprintf(filename2, "RESLT/sum_drag_and_torque%d.dat", i);
    file2.open(filename2);

    // Document the sum of drag and torque on the entire beam structure
    Hao_element_pt->output(file2);
    file2.close();
  }


} // end of parameter study

//========start_of_main================================================
/// Driver for beam (string under tension) test problem
//=====================================================================
int main()
{
  // Set the non-dimensional thickness
  // Global_Physical_Variables::H=0.01;
  Global_Physical_Variables::H = 0.25;

  // Set the 2nd Piola Kirchhoff prestress
  // Global_Physical_Variables::Sigma0=0.1;
  Global_Physical_Variables::Sigma0 = 0.1;

  // Set the length of domain
  double L = 10.0;

  // Number of elements (choose an even number if you want the control point
  // to be located at the centre of the beam)
  unsigned n_element = 100;

  // Construst the problem
  ElasticBeamProblem problem(n_element, L);

  // Check that we're ready to go:
  cout << "\n\n\nProblem self-test ";
  if (problem.self_test() == 0)
  {
    cout << "passed: Problem can be solved." << std::endl;
  }
  else
  {
    throw OomphLibError(
      "Self test failed", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  // Conduct parameter study
  problem.parameter_study();

} // end of main
