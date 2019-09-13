#ifndef _OPENSIM_VOLUMETRIC_CONTACT_FORCE_PLUGIN_H_
#define _OPENSIM_VOLUMETRIC_CONTACT_FORCE_PLUGIN_H_
/* -------------------------------------------------------------------------- *
 *                        VolumetricContactForce.h                            *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s):                                                                 *
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
// Headers define the various property types that OpenSim objects can read 
#include <string>
#include "osimPluginDLL.h"
#include <OpenSim/Simulation/Model/Force.h>
#include <OpenSim/Simulation/Model/TwoFrameLinker.h>


//=============================================================================
//=============================================================================
/*
 * A class template for writing a custom Force plugin. 
 * Applies ...
 * @author         
 * @version    
 */
namespace OpenSim {
    // TODO: OSIMPLUGIN_API will change to OSEMSIMULATION_API
class OSIMPLUGIN_API VolumetricContactForce 
    : public TwoFrameLinker<Force, PhysicalFrame>
{
    OpenSim_DECLARE_CONCRETE_OBJECT(VolumetricContactForce, TwoFrameLinker);
public:
//=============================================================================
// PROPERTIES
//=============================================================================
    // TODO: look for an osimsimulation library registered type, and attach it to that library instead of a separate one?

    OpenSim_DECLARE_PROPERTY(k_V, double,
        "Volumetric stiffness (N/m^3)");
    OpenSim_DECLARE_PROPERTY(a_V, double,
        "Volumetric damping");
    OpenSim_DECLARE_PROPERTY(mu_s, double,
        "Static friction");
    OpenSim_DECLARE_PROPERTY(mu_d, double,
        "Dynamic friction");
    OpenSim_DECLARE_PROPERTY(v_t, double,
        "Friction transition velocity (m/s)");
    OpenSim_DECLARE_PROPERTY(w_t, double,
        "Friction angular transition velocity (rad/s)");
    OpenSim_DECLARE_PROPERTY(ellipsoid_dims, SimTK::Vec3,
        "Ellipsoid dimensions (x, y, z axis) (m)");

//=============================================================================
// METHODS
//=============================================================================
public:
    // Default Constructor
    VolumetricContactForce();

    VolumetricContactForce(const std::string& name,
        const PhysicalFrame& ellipsoidFrame, const SimTK::Transform& transformInEllipsoidFrame,
        const PhysicalFrame& planeFrame, const SimTK::Transform& transformInPlaneFrame,
        const double& k_V,
        const double& a_V,
        const double& mu_s,
        const double& mu_d,
        const double& v_t,
        const double& w_t,
        const SimTK::Vec3& ellipsoid_dims);

    //--------------------------------------------------------------------------
    // COMPUTATION
    //--------------------------------------------------------------------------
    /** 
      */
    void computeForce(const SimTK::State& s, 
                      SimTK::Vector_<SimTK::SpatialVec>& bodyForces, 
                      SimTK::Vector& generalizedForces) const override;

    /** This may be more difficult for volumetric... */
    double computePotentialEnergy(const SimTK::State& s) const override;

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

protected:
    // TODO: is there any reason this should be implemented? It was in BodyDragForce
    //virtual void connectToModel(Model& aModel);


private:
    void setNull();
    void constructProperties();

    double findMu(const double &v, const double &vt) const;

    
//=============================================================================
};  // END of class VolumetricContactForce
//=============================================================================
//=============================================================================

} // end of namespace OpenSim

#endif // _OPENSIM_VOLUMETRIC_CONTACT_FORCE_PLUGIN_H_


