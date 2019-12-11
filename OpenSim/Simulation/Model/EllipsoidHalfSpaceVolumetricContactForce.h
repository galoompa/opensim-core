#ifndef _OPENSIM_ELLIPSOID_HALF_SPACE_VOLUMETRIC_CONTACT_FORCE_H_
#define _OPENSIM_ELLIPSOID_HALF_SPACE_VOLUMETRIC_CONTACT_FORCE_H_
/* -------------------------------------------------------------------------- *
 *             EllipsoidHalfSpaceVolumetricContactForce.h                     *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Peter Brown                                                     *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied    *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

//=============================================================================
// INCLUDES
//=============================================================================
#include "OpenSim/Simulation/osimSimulationDLL.h"
#include "OpenSim/Simulation/Model/Force.h"
#include "OpenSim/Simulation/Model/PhysicalFrame.h"
#include "OpenSim/Simulation/Model/ContactEllipsoid.h"

// TODO: fill in. For example of comments, see https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1Joint.html#details
// for an example with decent documentation
/**
This class implements a volumetric contact model, which is based on elastic
foundation or "bed of springs" contact model.

Volumetric contact relies on an analytical derivation of the volume of
penetration (and some other geometry) between two surfaces, which can result
in single set of equations describing the contact forces between two surfaces,
without needing to discretise the surfaces. Volumetric was first proposed in [1].

The ellipsoid-plane implementation of volumetric contact, along with the velocity-based
friction model, is described in more detail in [2,3].

[1] Gonthier, Y., McPhee, J., Lange, C., & Piedbœuf, J.-C. (2004). A Regularized Contact Model With Asymmetric Damping and Dwell-Time Dependent Friction. Multibody System Dynamics, 11(3), 209–233.
[2] Brown, P., & McPhee, J. (2017). Volumetric Contact Model of Ellipsoid-Plane Geometries. ECCOMAS Thematic Conference on Multibody Dynamics, 365–374.
[3] Brown, P. (2017). Contact Modelling for Forward Dynamics of Human Motion (University of Waterloo).

<b>C++ example</b>
\code{.cpp}
// Define contact between an ellipsoidal body and the ground plane
EllipsoidHalfSpaceVolumetricContactForce *eContact2 = new EllipsoidHalfSpaceVolumetricContactForce("ellipsoid_plane_contact",
        *ellipsoidBody, Transform(),
        model.getGround(), Transform(Rotation(-Pi/2, CoordinateAxis::XCoordinateAxis())), // transformed Z-axis should point up (plane normal)
        1e3, // volumetric stiffness
        -1, // volumetric damping
        0.51, // static friction
        0.5, // dynamic friction
        1e-4, // friction transition velocity
        1e-4, // friction transition angular velocity
        Vec3(a, b, c)); // ellipsoid dimensions
    model.addForce(eContact2);
\endcode

@author Peter Brown
*/
namespace OpenSim {
class OSIMSIMULATION_API EllipsoidHalfSpaceVolumetricContactForce 
    : public Force
{
    OpenSim_DECLARE_CONCRETE_OBJECT(EllipsoidHalfSpaceVolumetricContactForce, Force);
public:
//=============================================================================
// PROPERTIES
//=============================================================================
    OpenSim_DECLARE_LIST_PROPERTY(frames, PhysicalFrame,
        "Offset frames for internal use");
    OpenSim_DECLARE_PROPERTY(volumetricStiffness, double,
        "Volumetric stiffness (N/m^3)");
    OpenSim_DECLARE_PROPERTY(volumetricDamping, double,
        "Volumetric damping");
    OpenSim_DECLARE_PROPERTY(staticFrictionCoefficient, double,
        "Static friction");
    OpenSim_DECLARE_PROPERTY(dynamicFrictionCoefficient, double,
        "Dynamic friction");
    OpenSim_DECLARE_PROPERTY(frictionTransitionVelocity, double,
        "Friction transition velocity (m/s)");
    OpenSim_DECLARE_PROPERTY(frictionTransitionAngularVelocity, double,
        "Friction angular transition velocity (rad/s)");
    OpenSim_DECLARE_PROPERTY(ellipsoidDimensions, SimTK::Vec3,
        "Ellipsoid dimensions (x, y, z axis) (m)");

    //OpenSim_DECLARE_SOCKET(contactEllipsoid, ContactEllipsoid,
    //    "The contact ellipsoid geometry.");

//==============================================================================
// SOCKETS
//==============================================================================
    OpenSim_DECLARE_SOCKET(ellipsoidFrame, PhysicalFrame,
        "Frame connected to the centre of the ellipsoid geometry");
    OpenSim_DECLARE_SOCKET(halfspaceFrame, PhysicalFrame,
        "Frame connected to the halfspace (plane), on the surface, with the positive Z-axis pointing away from the half-space.");

//=============================================================================
// METHODS
//=============================================================================
public:
    // Default Constructor
    EllipsoidHalfSpaceVolumetricContactForce();

    /** 
    Create an ellipsoid-half-space volumetric contact force between a body with ellipsoid geometry and a plane.
    @name                                   the name of this contact
    @ellipsoidFrame                         the frame that the ellipsoid geometry is attached to
    @transformInEllipsoidFrame              the transform from ellipsoidFrame to the centre of the ellipsoid
                                            including the rotation between the ellipsoid primary axis and the frame axis
    @halfSpaceFrame                         the frame that the plane/half-space geometry is attached to
    @transformInHalfSpaceFrame              the transform from halfSpaceFrame to a point on the plane,
                                            including the rotation (note that the positive Z-axis points away from the plane).
    @volumetricStiffness                    the volumetric stiffness parameter of the contact (N/m^3)
    @volumetricDamping                      the volumetric damping parameter of the contact
    @staticFrictionCoefficient              the static coefficient of friction of the contact
    @dynamicFrictionCoefficient             the dynamic coefficient of friction of the contact
    @frictionTransitionVelocity             the transition velocity for linear friction (m/s)
    @frictionTransitionAngularVelocity      the transition velocity for spinning friction (rad/s)
    @ellipsoidDimensions                    the radii of the ellipsoid geometry (along the local x-, y-, and z-axis, respectively)
    */
    EllipsoidHalfSpaceVolumetricContactForce(const std::string& name,
        const PhysicalFrame& ellipsoidFrame, const SimTK::Transform& transformInEllipsoidFrame,
        const PhysicalFrame& halfSpaceFrame, const SimTK::Transform& transformInHalfSpaceFrame,
        const double& volumetricStiffness,
        const double& volumetricDamping,
        const double& staticFrictionCoefficient,
        const double& dynamicFrictionCoefficient,
        const double& frictionTransitionVelocity,
        const double& frictionTransitionAngularVelocity,
        const SimTK::Vec3& ellipsoidDimensions);

    //--------------------------------------------------------------------------
    // COMPUTATION
    //--------------------------------------------------------------------------
    /** 
    Calculates the force between the two bodies
      */
    void computeForce(const SimTK::State& s, 
                      SimTK::Vector_<SimTK::SpatialVec>& bodyForces, 
                      SimTK::Vector& generalizedForces) const override;

    
    // TODO: consider using ContactHalfSpace and make a ContactEllipsoid

    //-----------------------------------------------------------------------------
    // Reporting
    //-----------------------------------------------------------------------------
    /** 
     * Provide name(s) of the quantities (column labels) of the force value(s) to be reported
     */
    OpenSim::Array<std::string> getRecordLabels() const override;
    /**
    *  Provide the value(s) to be reported that correspond to the labels
    */
    OpenSim::Array<double> getRecordValues(const SimTK::State& state) const override;

    void generateDecorations(bool fixed, const ModelDisplayHints& hints,
        const SimTK::State& s, SimTK::Array_<SimTK::DecorativeGeometry>& geometry) const;

private:
    void setNull();
    void constructProperties();

    /**
    Calculates the value for mu based on the current velocity, the transition velocity, and the
    values for the static and dynamic coefficients of friction.
    */
    double findMu(const double &v, const double &vt) const;

    
//=============================================================================
};  // END of class EllipsoidHalfSpaceVolumetricContactForce
//=============================================================================

} // end of namespace OpenSim

#endif // _OPENSIM_ELLIPSOID_HALF_SPACE_VOLUMETRIC_CONTACT_FORCE_H_


