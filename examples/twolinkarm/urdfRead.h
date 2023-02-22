/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin@fysx.org>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <cmath>
#include <rbdl/rbdl.h>
#include <rbdl/rbdl_utils.h>

#ifndef RBDL_BUILD_ADDON_URDFREADER
#error "Error: RBDL addon URDFReader not enabled."
#endif

#include <rbdl/addons/urdfreader/urdfreader.h>

namespace flapping_model{
	const double p = 1.1839;
	const double wing_area = 0.046;
	const double tail_area = 0.018;
	const double C_l0 = 1.0;
	const double D_l0 = 1.0;
	const double D_l1 = 1.0;
    const double v0 = 5.0;
    const double target_height = 10.0;
};

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

void CalF(Model& model, const VectorNd& Q, const VectorNd& QDot, int body_id,
			std::vector<SpatialVector>& fext, const Matrix3d& rotate_y){
	auto V_inertial = CalcPointVelocity(model, Q, QDot, body_id, {0, 0, 0}, true);
	auto MatWorld2Body = CalcBodyWorldOrientation(model, Q, body_id, false);
	auto V_local = MatWorld2Body * V_inertial;// CalcBodyWorldOrientation
	Vector3d V_local_proj{V_local[0], 0, V_local[2]};
	double angle_of_attack;
	if(V_local[0] == 0){
		double angle_of_attack = atan(INFINITY);
	}
	else{
		double angle_of_attack = atan(V_local[2] / V_local[0]);
	}
	double C_l = flapping_model::C_l0 * sin(2 * angle_of_attack);
	double C_d = flapping_model::D_l0 - flapping_model::D_l1 * cos(2 * angle_of_attack);
	double area = flapping_model::wing_area;
	if(body_id == 7){
		area = flapping_model::tail_area;
	}
	double F_l = 2.0 / 3.0 * C_l * flapping_model::p * area * (V_local_proj.norm() * V_local_proj.norm());
	double F_d = 2.0 / 3.0 * C_d * flapping_model::p * area * (V_local_proj.norm() * V_local_proj.norm());

	auto pos = CalcBodyToBaseCoordinates(model, Q, body_id, {0, 0, 0}, false);
	// std::cout << "pos: " << pos << std::endl;
	Matrix3d D;
	D << 0, -pos[2], pos[1], pos[2], 0, -pos[0], -pos[1], pos[0], 0;

	auto F = F_l * MatWorld2Body.transpose() * (rotate_y * V_local_proj.normalized())
			- F_d * (MatWorld2Body.transpose() * V_local_proj).normalized();
	// std::cout << "F: " << F << std::endl;
	auto T = D * F;
	// std::cout << "T: " << T << std::endl;
	fext[body_id][0] = T[0];
	fext[body_id][1] = T[1];
	fext[body_id][2] = T[2];
	fext[body_id][3] = F[0];
	fext[body_id][4] = F[1];
	fext[body_id][5] = F[2];
}

void urdfRead (vector<adouble> x, adouble u, VectorNd& QDDot) {
	rbdl_check_api_version (RBDL_API_VERSION);
	Model* model = new Model();
	if (!Addons::URDFReadFromFile ("/home/lj/new_model/urdf/model.urdf", model,true, true)) {
		std::cerr << "Error loading model " << std::endl;
		abort();
	}

	// std::cout << "Degree of freedom overview:" << std::endl;
	// std::cout << Utils::GetModelDOFOverview(*model);

	// std::cout << "Model Hierarchy:" << std::endl;
	// std::cout << Utils::GetModelHierarchy(*model);

	VectorNd Q = VectorNd::Zero (model->q_size);
    Q << double(x[1]), 0.0, double(x[2]), 0.0, double(x[0]), 0.0, 0.0, 0.0, 0.0, 0.0, double(x[3]);

	VectorNd QDot = VectorNd::Zero (model->qdot_size);
	QDot << double(x[5]), 0.0, double(x[6]), 0.0, double(x[4]), 0.0, 0.0, 0.0, 0.0, 0.0, double(x[7]);

	VectorNd Tau = VectorNd::Zero (model->qdot_size);
    Tau[10] = double(u);

	QDDot = VectorNd::Zero (model->qdot_size);
	std::vector<RigidBodyDynamics::Math::SpatialVector> fext;
	fext.resize(model->mBodies.size());
	for (unsigned int i = 0; i < fext.size(); ++i) {
		fext[i]=RigidBodyDynamics::Math::SpatialVector::Zero();
	}
	Matrix3d rotate_y;
	rotate_y << 0, 0, -1, 0, 1, 0, 1, 0, 0;
	CalF(*model, Q, QDot, model->GetBodyId("left_wing"), fext, rotate_y);
	CalF(*model, Q, QDot, model->GetBodyId("right_wing"), fext, rotate_y);
	CalF(*model, Q, QDot, model->GetBodyId("tail"), fext, rotate_y);

	// std::cout << "body size: " << model->mBodies.size() << std::endl;
	// std::cout << "left: " << model->GetBodyId("left") << std::endl;
	// std::cout << "left_wing: " << model->GetBodyId("left_wing") << std::endl;
	// std::cout << "right: " << model->GetBodyId("right") << std::endl;
	// std::cout << "right_wing: " << model->GetBodyId("right_wing") << std::endl;
	// std::cout << "tail: " << model->GetBodyId("tail") << std::endl;

	ForwardDynamics (*model, Q, QDot, Tau, QDDot, &fext);
	// InverseDynamics (*model, Q, QDot, QDDot, Tau);
	// std::cout << QDDot.transpose() << std::endl;

	delete model;

}

