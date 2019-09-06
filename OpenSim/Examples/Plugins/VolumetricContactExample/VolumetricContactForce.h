#ifndef _OPENSIM_VOLUMETRIC_CONTACT_FORCE_PLUGIN_H_
#define _OPENSIM_VOLUMETRIC_CONTACT_FORCE_PLUGIN_H_
/* -------------------------------------------------------------------------- *
 *                             VolumetricContactForce.h                                *
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


//=============================================================================
//=============================================================================
/*
 * A class template for writing a custom Force plugin. 
 * Applies ...
 * @author         
 * @version    
 */
namespace OpenSim {

class OSIMPLUGIN_API VolumetricContactForce : public Force  
{
OpenSim_DECLARE_CONCRETE_OBJECT(VolumetricContactForce, Force);
public:
//=============================================================================
// PROPERTIES
//=============================================================================

    /** String property containing the name of the ellipsoid body*/
    OpenSim_DECLARE_PROPERTY(ellipsoid_name, std::string,
        "Ellipsoid body name.");

    /** String property containing the name of the plane body (positive Z-axis is the plane normal)*/
    OpenSim_DECLARE_PROPERTY(plane_name, std::string,
        "Plane body name.");

    // TODO: perhaps the above body names should be replaced with a socket for the bodies

    /** Double property containing the volumetric stiffness (k_V)*/
    OpenSim_DECLARE_PROPERTY(k_V, double,
        "Volumetric stiffness");

    /** Double property containing the volumetric damping (a_V)*/
    OpenSim_DECLARE_PROPERTY(a_V, double,
        "Volumetric damping");

    /** Double property containing the static friction (mu_s)*/
    OpenSim_DECLARE_PROPERTY(mu_s, double,
        "Static friction");

    /** Double property containing the dynamic friction (mu_d)*/
    OpenSim_DECLARE_PROPERTY(mu_d, double,
        "Dynamic friction");

    /** Double property containing the linear transition velocity (for friction) (v_t)*/
    OpenSim_DECLARE_PROPERTY(v_t, double,
        "Friction transition velocity");

    /** Double property containing the angular transition velocity (for friction) (w_t)*/
    OpenSim_DECLARE_PROPERTY(w_t, double,
        "Friction angular transition velocity");

    /** Vector property containing the ellipsoid dimensions */
    OpenSim_DECLARE_PROPERTY(ellipsoid_dims, SimTK::Vec3,
        "Ellipsoid dimensions (x, y, z axis)");

    // Here are some examples of other scalar property types.
    // Uncomment them as you need them.
    // ------------------------------------------------------
    //// My string property
    //OpenSim_DECLARE_PROPERTY(string_property, std::string, 
    //"My string property."); 

    //// My int property
    //OpenSim_DECLARE_PROPERTY(int_property, int, 
    //"My int property."); 

    //// My bool property
    //OpenSim_DECLARE_PROPERTY(bool_property, bool, 
    //"My bool property."); 

    //// My double property
    //OpenSim_DECLARE_PROPERTY(double_property, double, 
    //"My double property."); 


//=============================================================================
// METHODS
//=============================================================================
public:
    // Default Constructor
    VolumetricContactForce();


    //--------------------------------------------------------------------------
    // COMPUTATION
    //--------------------------------------------------------------------------
    /** Compute the bushing force contribution to the system and add in to appropriate
      * bodyForce and/or system generalizedForce. The bushing force is [K]*dq + [D]*dqdot
      * where, [K] is the spatial 6dof stiffness matrix between the two frames 
               dq is the deflection in body spatial coordinates with rotations in Euler angles
      *        [D] is the spatial 6dof damping matrix opposing the velocity between the frames
      *        dqdot is the relative spatial velocity of the two frames
      * VolumetricContactForce implementation based SimTK::Force::LinearBushing
      * developed and implemented by Michael Sherman.
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
    OpenSim::Array<double>
    getRecordValues(const SimTK::State& state) const override;

protected:
    virtual void connectToModel(Model& aModel);


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


