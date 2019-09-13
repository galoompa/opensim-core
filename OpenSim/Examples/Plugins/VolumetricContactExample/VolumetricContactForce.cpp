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
VolumetricContactForce::VolumetricContactForce()
    : Super()
//    : TwoFrameLinker<Force, PhysicalFrame>()
{
    setNull();
    constructProperties();
}

VolumetricContactForce::VolumetricContactForce(const std::string& name,
    const PhysicalFrame& ellipsoidFrame, const SimTK::Transform& transformInEllipsoidFrame,
    const PhysicalFrame& planeFrame, const SimTK::Transform& transformInPlaneFrame,
    const double& k_V,
    const double& a_V,
    const double& mu_s,
    const double& mu_d,
    const double& v_t,
    const double& w_t,
    const SimTK::Vec3& ellipsoid_dims)
    : Super(name, ellipsoidFrame, transformInEllipsoidFrame, planeFrame, transformInPlaneFrame)
{
    setNull();
    constructProperties();
    set_k_V(k_V);
    set_a_V(a_V);
    set_mu_s(mu_s);
    set_mu_d(mu_d);
    set_v_t(v_t);
    set_w_t(w_t);
    set_ellipsoid_dims(ellipsoid_dims);
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
    constructProperty_k_V(1e4);
    constructProperty_a_V(-1);
    constructProperty_ellipsoid_dims(Vec3(1));
    constructProperty_mu_d(0.5);
    constructProperty_mu_s(0.5);
    constructProperty_v_t(1e-6);
    constructProperty_w_t(1e-6);
}

//_____________________________________________________________________________
/**
 * Perform some set up functions that happen after the
 * object has been deserialized or copied.
 *
 * @param aModel OpenSim model containing this VolumetricContactForce.
 */
#if 0
void VolumetricContactForce::connectToModel(Model& aModel)
{
    string errorMessage;

    // Base class
    Super::connectToModel(aModel);

    // TODO: probably remove below (or this whole function)

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
#endif


//=============================================================================
// COMPUTATION
//=============================================================================

void VolumetricContactForce::computeForce(const SimTK::State& s, 
                              SimTK::Vector_<SimTK::SpatialVec>& bodyForces, 
                              SimTK::Vector& generalizedForces) const
{
    // Some references to simplify equations
    const Ground& ground = getModel().getGround(); // ground reference
    const PhysicalFrame& eFrame = getFrame1(); // ellipsoid frame reference
    const PhysicalFrame& pFrame = getFrame2(); // plane body reference

    // Find and name variables to match ellipsoid equations exported from Maple
    SimTK::Rotation eR = eFrame.getTransformInGround(s).R();
    eR = eFrame.findTransformBetween(s, pFrame).R();
    Vec3 pp_e = pFrame.findStationLocationInAnotherFrame(s, Vec3(0), eFrame);
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
    //Vec3 centroidPos = eFrame.expressVectorInGround(s, Vec3(cenLx, cenLy, cenLz)); // centroid expressed in ground frame
    Vec3 centroidPosP = eFrame.expressVectorInAnotherFrame(s, Vec3(cenLx, cenLy, cenLz), pFrame); // centroid expressed in plane frame

    // NOTE: The rest of this function could be generalised for any volumetric contact.

    // Relative linear and angular velocity (at centroid) expressed in plane frame
    Vec3 v_E_cen_G = eFrame.findStationVelocityInGround(s, centroidPosE); // velocity of ellipsoid at centroid, expressed in ground frame
    Vec3 v_P_cen_G = pFrame.findStationVelocityInGround(s, centroidPosP); // velocity of plane at centroid, expressed in ground frame
    Vec3 v_cen_P = ground.expressVectorInAnotherFrame(s, v_E_cen_G - v_P_cen_G, pFrame); // relative velocity (of ellipsoid) at centroid, expressed in plane frame
    Vec3 w_P = ground.expressVectorInAnotherFrame(s, // relative angular velocity, expressed in plane frame
        eFrame.getAngularVelocityInGround(s) - pFrame.getAngularVelocityInGround(s), pFrame);
    
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
    Vec3 forceG = pFrame.expressVectorInGround(s, Vec3(forceXP, forceYP, forceZP)); // force vector in ground frame
    Vec3 torqueG = pFrame.expressVectorInGround(s, Vec3(torqueXP, torqueYP, torqueZP)); // torque vector in ground frame
    
    // Apply the forces and torques to each body:
    applyForceToPoint(s, eFrame, centroidPosE, forceG, bodyForces);
    applyForceToPoint(s, pFrame, centroidPosP, -forceG, bodyForces);
    applyTorque(s, eFrame, torqueG, bodyForces);
    applyTorque(s, pFrame, -torqueG, bodyForces);

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
