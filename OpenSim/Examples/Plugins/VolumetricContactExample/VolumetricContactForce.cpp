/* -------------------------------------------------------------------------- *
 *                           VolumetricContactForce.cpp                                *
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
#include <OpenSim/Simulation/Model/BodySet.h>
#include <OpenSim/Simulation/Model/Model.h>
#include "VolumetricContactForce.h"

//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace SimTK;
using namespace OpenSim;

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/**
 * Default constructor
 */
VolumetricContactForce::VolumetricContactForce() : Force()
{
    setNull();
    constructProperties();
}



//=============================================================================
// CONSTRUCTION
//=============================================================================
//_____________________________________________________________________________
/**
 * Set the data members of this VolumetricContactForce to their null values.
 */
void VolumetricContactForce::setNull()
{
    // no internal data members to initialize.
}

//_____________________________________________________________________________
/**
 * Connect properties to local pointers.
 */
void VolumetricContactForce::constructProperties()
{
    constructProperty_ellipsoid_name("Unassigned");
    constructProperty_plane_name("Unassigned");
    constructProperty_k_V(1e4);
    constructProperty_a_V(-1);
    constructProperty_ellipsoid_dims(Vec3(1));
    constructProperty_mu_d(0.5);
    constructProperty_mu_s(0.5);
    constructProperty_v_t(1e-6);
    constructProperty_w_t(1e-6);

    // Here are some examples of other constructing other scalar property types.
    // Uncomment them as you need them.
    // ------------------------------------------------------
    //constructProperty_string_property("defaultString");
    //constructProperty_int_property(10);
    //constructProperty_bool_property(true);
    //constructProperty_double_property(1.5);
}

//_____________________________________________________________________________
/**
 * Perform some set up functions that happen after the
 * object has been deserialized or copied.
 *
 * @param aModel OpenSim model containing this VolumetricContactForce.
 */
void VolumetricContactForce::connectToModel(Model& aModel)
{
    string errorMessage;

    // Base class
    Super::connectToModel(aModel);

    // Look up the body and report an error if it is not found 
    if (!aModel.updBodySet().contains(get_ellipsoid_name())) {
        errorMessage = "Invalid bodyName (" + get_ellipsoid_name() + ") specified in Force " + getName();
        throw (Exception(errorMessage.c_str()));
    }
    if (!aModel.updBodySet().contains(get_plane_name())) {
        errorMessage = "Invalid bodyName (" + get_plane_name() + ") specified in Force " + getName();
        throw (Exception(errorMessage.c_str()));
    }
}


//=============================================================================
// COMPUTATION
//=============================================================================

void VolumetricContactForce::computeForce(const SimTK::State& s, 
                              SimTK::Vector_<SimTK::SpatialVec>& bodyForces, 
                              SimTK::Vector& generalizedForces) const
{
    if(!_model) return;     // some minor error checking

    const BodySet& bs = _model->updBodySet();         // get body set
    const Body& eBody = bs.get(get_ellipsoid_name()); // The ellipsoid body
    const Body& pBody = bs.get(get_plane_name()); // The plane body
    const Ground& ground = _model->getGround(); // reference to ground

    // Find and name variables to match ellipsoid equations exported from Maple
    SimTK::Rotation eR = eBody.getTransformInGround(s).R();
    eR = eBody.findTransformBetween(s, pBody).R();
    Vec3 pp_e = pBody.findStationLocationInAnotherFrame(s, Vec3(0), eBody);
    double R1 = eR.x()[0];
    double R2 = eR.y()[0];
    double R3 = eR.z()[0];
    double R4 = eR.x()[1];
    double R5 = eR.y()[1];
    double R6 = eR.z()[1];
    double R7 = eR.x()[2];
    double R8 = eR.y()[2];
    double R9 = eR.z()[2];
    double xp = pp_e[0];
    double yp = pp_e[1];
    double zp = pp_e[2];
    double a = get_ellipsoid_dims()[0];
    double b = get_ellipsoid_dims()[1];
    double c = get_ellipsoid_dims()[2];

    // Volumetric equations for ellipsoids (derived and optimised in Maple)
    double V = (0.0e0 < -0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0.4e1 / 0.3e1 * 0.3141592654e1 * a * b * c : (0.0e0 <= 0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? -0.3141592654e1 * pow(0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.2e1) * (-0.2e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1)) * a * b * c / 0.3e1 : 0));
    double cenLx = (0.0e0 < -0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0 : (0.0e0 <= 0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0.3e1 * a * a * pow(-0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.2e1) / (-0.8e1 + 0.4e1 * (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1)) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) * R7 : 0));
    double cenLy = (0.0e0 < -0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0 : (0.0e0 <= 0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0.3e1 * b * b * pow(-0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.2e1) / (-0.8e1 + 0.4e1 * (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1)) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) * R8 : 0));
    double cenLz = (0.0e0 < -0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0 : (0.0e0 <= 0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0.3e1 * c * c * pow(-0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.2e1) / (-0.8e1 + 0.4e1 * (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1)) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) * R9 : 0));
    double Ji = (0.0e0 < -0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0.4e1 / 0.15e2 * 0.3141592654e1 : (0.0e0 <= 0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0.3141592654e1 * pow(0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.3e1) * (0.3e1 * pow(0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.2e1) + 0.5e1 - 0.15e2 * (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1)) / 0.60e2 : 0));
    double Jp_xx = (0.0e0 <= 0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0.3141592654e1 * pow(0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.3e1) * (0.3e1 * pow(0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.2e1) + 0.5e1 - 0.15e2 * (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1)) * pow(((R1 * R1 + R2 * R2 - 0.1e1) * a * a - c * c * R1 * R1) * b * b - R2 * R2 * a * a * c * c, 0.2e1) * b * a * c / (((pow(R2, 0.4e1) + (R1 * R1 - R4 * R4 + R5 * R5 - 0.1e1) * R2 * R2 + 0.2e1 * R1 * R2 * R4 * R5 + R4 * R4) * a * a - ((R1 * R1 - R4 * R4 + 0.1e1) * R2 * R2 + 0.2e1 * R1 * R2 * R4 * R5 + R4 * R4 + R5 * R5 - 0.1e1) * c * c) * pow(b, 0.4e1) + ((R1 * R1 + R4 * R4 - 0.1e1) * (R1 * R1 + R2 * R2 - 0.1e1) * pow(a, 0.4e1) - 0.2e1 * (pow(R2, 0.4e1) + (R1 * R1 + R5 * R5 - 0.1e1) * R2 * R2 + R1 * R2 * R4 * R5 + R1 * R1 * (R1 * R1 + R4 * R4 - 0.1e1)) * c * c * a * a + ((R1 * R1 - R4 * R4 + 0.1e1) * R2 * R2 + 0.2e1 * R1 * R2 * R4 * R5 + pow(R1, 0.4e1) + R1 * R1 * R4 * R4 + R4 * R4 + R5 * R5 - 0.1e1) * pow(c, 0.4e1)) * b * b - R2 * R2 * ((R1 * R1 + R4 * R4 - 0.1e1) * a * a - c * c * (R1 * R1 + R2 * R2 + R4 * R4 + R5 * R5 - 0.1e1)) * a * a * c * c) / 0.60e2 : 0);
    double Jp_yy = (0.0e0 <= 0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0.3141592654e1 * pow(0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.3e1) * (0.3e1 * pow(0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.2e1) + 0.5e1 - 0.15e2 * (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1)) * ((pow(R1 * R4 + R2 * R5, 0.2e1) * pow(a, 0.4e1) - 0.2e1 * (R1 * R1 * R4 * R4 + R1 * R2 * R4 * R5 + R2 * R2 / 0.2e1 + R5 * R5 / 0.2e1 - 0.1e1 / 0.2e1) * c * c * a * a + pow(c, 0.4e1) * R1 * R1 * R4 * R4) * pow(b, 0.4e1) - ((0.2e1 * R1 * R2 * R4 * R5 + 0.2e1 * R2 * R2 * R5 * R5 + R1 * R1 + R4 * R4 - 0.1e1) * a * a - c * c * (0.2e1 * R1 * R2 * R4 * R5 + R1 * R1 + R2 * R2 + R4 * R4 + R5 * R5 - 0.1e1)) * a * a * c * c * b * b + pow(a, 0.4e1) * pow(c, 0.4e1) * R2 * R2 * R5 * R5) * b * a * c / (((pow(R2, 0.4e1) + (R1 * R1 - R4 * R4 + R5 * R5 - 0.1e1) * R2 * R2 + 0.2e1 * R1 * R2 * R4 * R5 + R4 * R4) * a * a - ((R1 * R1 - R4 * R4 + 0.1e1) * R2 * R2 + 0.2e1 * R1 * R2 * R4 * R5 + R4 * R4 + R5 * R5 - 0.1e1) * c * c) * pow(b, 0.4e1) + ((R1 * R1 + R4 * R4 - 0.1e1) * (R1 * R1 + R2 * R2 - 0.1e1) * pow(a, 0.4e1) - 0.2e1 * (pow(R2, 0.4e1) + (R1 * R1 + R5 * R5 - 0.1e1) * R2 * R2 + R1 * R2 * R4 * R5 + R1 * R1 * (R1 * R1 + R4 * R4 - 0.1e1)) * c * c * a * a + ((R1 * R1 - R4 * R4 + 0.1e1) * R2 * R2 + 0.2e1 * R1 * R2 * R4 * R5 + pow(R1, 0.4e1) + R1 * R1 * R4 * R4 + R4 * R4 + R5 * R5 - 0.1e1) * pow(c, 0.4e1)) * b * b - R2 * R2 * ((R1 * R1 + R4 * R4 - 0.1e1) * a * a - c * c * (R1 * R1 + R2 * R2 + R4 * R4 + R5 * R5 - 0.1e1)) * a * a * c * c) / 0.60e2 : 0);
    double Jp_xy = (0.0e0 <= 0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1) ? 0.3141592654e1 * pow(0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.3e1) * (0.3e1 * pow(0.1e1 + (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1), 0.2e1) + 0.5e1 - 0.15e2 * (R7 * xp + R8 * yp + R9 * zp) * pow(a * a * R7 * R7 + b * b * R8 * R8 + c * c * R9 * R9, -0.1e1 / 0.2e1)) * (((R1 * R1 + R2 * R2 - 0.1e1) * a * a - c * c * R1 * R1) * b * b - R2 * R2 * a * a * c * c) * b * (((R1 * R4 + R2 * R5) * a * a - c * c * R1 * R4) * b * b - a * a * c * c * R2 * R5) * a * c / (((pow(R2, 0.4e1) + (R1 * R1 - R4 * R4 + R5 * R5 - 0.1e1) * R2 * R2 + 0.2e1 * R1 * R2 * R4 * R5 + R4 * R4) * a * a - ((R1 * R1 - R4 * R4 + 0.1e1) * R2 * R2 + 0.2e1 * R1 * R2 * R4 * R5 + R4 * R4 + R5 * R5 - 0.1e1) * c * c) * pow(b, 0.4e1) + ((R1 * R1 + R4 * R4 - 0.1e1) * (R1 * R1 + R2 * R2 - 0.1e1) * pow(a, 0.4e1) - 0.2e1 * (pow(R2, 0.4e1) + (R1 * R1 + R5 * R5 - 0.1e1) * R2 * R2 + R1 * R2 * R4 * R5 + R1 * R1 * (R1 * R1 + R4 * R4 - 0.1e1)) * c * c * a * a + ((R1 * R1 - R4 * R4 + 0.1e1) * R2 * R2 + 0.2e1 * R1 * R2 * R4 * R5 + pow(R1, 0.4e1) + R1 * R1 * R4 * R4 + R4 * R4 + R5 * R5 - 0.1e1) * pow(c, 0.4e1)) * b * b - R2 * R2 * ((R1 * R1 + R4 * R4 - 0.1e1) * a * a - c * c * (R1 * R1 + R2 * R2 + R4 * R4 + R5 * R5 - 0.1e1)) * a * a * c * c) / 0.60e2 : 0);

    // Centroid position and velocity (universal equations)
    Vec3 centroidPosE = Vec3(cenLx, cenLy, cenLz); // centroid expressed in ellipsoid frame
    Vec3 centroidPos = eBody.expressVectorInGround(s, Vec3(cenLx, cenLy, cenLz)); // centroid expressed in ground frame
    Vec3 centroidPosP = eBody.expressVectorInAnotherFrame(s, Vec3(cenLx, cenLy, cenLz), pBody); // centroid expressed in plane frame

    // Relative linear and angular velocity (at centroid) expressed in plane frame
    Vec3 v_E_cen_G = eBody.findStationVelocityInGround(s, centroidPosE); // velocity of ellipsoid at centroid, expressed in ground frame
    Vec3 v_P_cen_G = pBody.findStationVelocityInGround(s, centroidPosP); // velocity of plane at centroid, expressed in ground frame
    Vec3 v_cen_P = ground.expressVectorInAnotherFrame(s, v_E_cen_G - v_P_cen_G, pBody); // relative velocity (of ellipsoid) at centroid, expressed in plane frame
    Vec3 w_P = ground.expressVectorInAnotherFrame(s, // relative angular velocity, expressed in plane frame
        eBody.getAngularVelocityInGround(s) - pBody.getAngularVelocityInGround(s), pBody);
    
    // Normal force equation
    double forceZP = get_k_V() * V * (1 + get_a_V() * v_cen_P[2]);

    // Rolling resistance equations (moment in tangential directions)
    double torqueXP = get_k_V() * get_a_V() * (Jp_xx * w_P[0] + Jp_xy * w_P[1]);
    double torqueYP = get_k_V() * get_a_V() * (Jp_xy * w_P[0] + Jp_yy * w_P[1]);

    // Tangential friction equations (force in tangential directions)
    double vt_cen_P = pow(pow(v_cen_P[0], 2) + pow(v_cen_P[1], 2), 0.5); // magnitude of relative velocity at centroid
    double forceXP, forceYP;
    if (vt_cen_P == 0) // avoid divide by 0 errors
    {
        forceXP = 0;
        forceYP = 0;
    }
    else
    {
        double mu = findMu(vt_cen_P, get_v_t());
        forceXP = -abs(forceZP) * mu * v_cen_P[0] / vt_cen_P;
        forceYP = -abs(forceZP) * mu * v_cen_P[1] / vt_cen_P;
    }
    
    // Spinning friction (torque in normal direction)
    double torqueZP;
    if (V == 0) // avoid divide by 0 errors
    {
        torqueZP = 0;
    }
    else
    {
        double mu = findMu(w_P[2], get_w_t());
        torqueZP = -abs(forceZP) / V * mu * (Jp_xx + Jp_yy); // note: Jp_zz = (Jp_xx + Jp_yy)
    }
    
    // Convert force and torque from plane frame to ground frame
    Vec3 forceG = pBody.expressVectorInGround(s, Vec3(forceXP, forceYP, forceZP)); // force vector in ground frame
    Vec3 torqueG = pBody.expressVectorInGround(s, Vec3(torqueXP, torqueYP, torqueZP)); // torque vector in ground frame
    
    // TODO: figure out correct convention...
    // The previous comments (see bodyDragExample) may be completely wrong... the point should be in the body frame, and the force in the ground frame
    applyForceToPoint(s, eBody, centroidPosE, forceG, bodyForces);
    applyForceToPoint(s, pBody, centroidPosP, -forceG, bodyForces);

    applyTorque(s, eBody, torqueG, bodyForces);
    applyTorque(s, pBody, -torqueG, bodyForces);

    // Debuging info
    // --------------
    int deb = 0;
    if (deb)
    {
        cout << "Time = " << s.getTime() << endl;
        printf("cenL %f %f %f\n", cenLx, cenLy, cenLz);
        printf("F %f %f %f\n", forceXP, forceYP, forceZP);
        printf("T %f %f %f\n", torqueXP, torqueYP, torqueZP);
        /*for (int i = 0; i < 3; ++i)
        {
            cout << eR.x()[i] << "\t" << eR.y()[i] << "\t" << eR.z()[i] << endl;
        }*/
        auto res = system("pause");
    }

    return;
}

double VolumetricContactForce::findMu(const double & v, const double & vt) const
{
    return get_mu_d() * tanh(4 * (v / vt)) + (get_mu_s() - get_mu_d()) * (v / vt) / pow(0.25 * pow((v / vt), 2) + 0.75, 2);
}

/** Potential energy function */
double VolumetricContactForce::computePotentialEnergy(const SimTK::State& s) const
{
    // TODO: implement, if necessary
    return 0;
}


//=============================================================================
// REPORTING
//=============================================================================
/** 
 * Provide names of the quantities (column labels) of the force value(s) reported
 * 
 */
OpenSim::Array<std::string> VolumetricContactForce::getRecordLabels() const 
{
    OpenSim::Array<std::string> labels("");
    //labels.append(getName()+"."+get_body_name()+".force.X");
    return labels;
}
/**
 * Provide the value(s) to be reported that correspond to the labels
 */
OpenSim::Array<double> VolumetricContactForce::getRecordValues(const SimTK::State& s) const 
{
    OpenSim::Array<double> values(3);
    // need to assign values...
    return values;
}
