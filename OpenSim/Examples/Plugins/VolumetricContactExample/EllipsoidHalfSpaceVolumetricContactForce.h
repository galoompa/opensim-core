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
#include <string>
#include "osimPluginDLL.h" // TODO: remove this
#include <OpenSim/Simulation/Model/Force.h>
#include <OpenSim/Simulation/Model/PhysicalFrame.h>

// TODO: fill in. For example of comments, see https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1Joint.html#details
// for an example with decent documentation
/**
This class implements a volumetric contact model, which is based on elastic
foundation or "bed of springs" contact model.

Volumetric contact relies on an analytical derivation of the volume of
penetration (and some other geometry) between two surfaces, which can result
in single set of equations describing the contact forces between two surfaces,
without needing to discretise the surfaces.

TODO: fill in
stomJoint). For more
details on how the underlying formulation supports coupled curvilinear
joints, see "Minimal formulation of joint motion for biomechanisms", 2010
A Seth, M Sherman, P Eastman, S Delp; Nonlinear dynamics 62 (1), 291-303

<b>C++ example</b>
\code{.cpp}
// Define a pin joint that attaches pendulum (an OpenSim::Body) to ground.
PinJoint* myPin = new PinJoint("pendulumToGround", myModel.getGround(),
                               pendulum);
\endcode

@author Peter Brown
*/
namespace OpenSim {
    // TODO: OSIMPLUGIN_API will change to OSEMSIMULATION_API
class OSIMPLUGIN_API EllipsoidHalfSpaceVolumetricContactForce 
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

    // TODO: comments / documentation
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
      */
    void computeForce(const SimTK::State& s, 
                      SimTK::Vector_<SimTK::SpatialVec>& bodyForces, 
                      SimTK::Vector& generalizedForces) const override;

    // TODO: consider implementing denerateDecorations (eg. see contactSphere)

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


private:
    void setNull();
    void constructProperties();

    double findMu(const double &v, const double &vt) const;

    
//=============================================================================
};  // END of class EllipsoidHalfSpaceVolumetricContactForce
//=============================================================================

} // end of namespace OpenSim

#endif // _OPENSIM_ELLIPSOID_HALF_SPACE_VOLUMETRIC_CONTACT_FORCE_H_


