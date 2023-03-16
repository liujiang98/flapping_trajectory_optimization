/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin@fysx.org>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */
#pragma once

#include <iostream>
#include <cmath>
#include <rbdl/rbdl.h>
#include <rbdl/rbdl_utils.h>

#ifndef RBDL_BUILD_ADDON_URDFREADER
#error "Error: RBDL addon URDFReader not enabled."
#endif

#include <rbdl/addons/urdfreader/urdfreader.h>

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
namespace flapping_model{
	const double p = 1.1839;
	const double wing_area = 0.046;
	const double tail_area = 0.018;
	const double C_l0 = 1.332;
	const double D_l0 = 1.713;
	const double D_l1 = -1.639;
	const Vector3d wing_pos{-0.03238, 0.0, 0.0};
	const Vector3d tail_pos{-0.10905, 0.0, 0.02466};
	const double v0 = 6.0;
    const double target_height = 0.8;
	const double f = 10.0; // 扑翼频率
};

void CalF(Model& model, const VectorNd& Q, const VectorNd& QDot, const char* body_name,
			std::vector<SpatialVector>& fext){
	int body_id = model.GetBodyId(body_name);
	double area = flapping_model::wing_area;
	Vector3d body_pos = flapping_model::wing_pos;
	if(body_id == 7){ // body_id of tail is 7
		area = flapping_model::tail_area;
		body_pos = flapping_model::tail_pos;
	}

	// 惯性系下的速度
	Vector3d V_inertial = CalcPointVelocity(model, Q, QDot, body_id, body_pos, true);
	// std::cout << "V_inertial: " << V_inertial.transpose() << std::endl;

	Matrix3d MatWorld2Body = CalcBodyWorldOrientation(model, Q, body_id, true);
	Eigen::Quaterniond q(MatWorld2Body);

	// 从惯性系到body(local)系下的旋转矩阵
	MatWorld2Body = q.normalized().toRotationMatrix();

	// body(local)系下的速度
	Vector3d V_local = MatWorld2Body * V_inertial;// CalcBodyWorldOrientation
	// std::cout << "V_local: " << V_local.transpose() << std::endl;

	Vector3d V_local_proj{V_local[0], 0, V_local[2]};
	double x = V_local[2] / V_local[0];
	double angle_of_attack;
	if(V_local[0] == 0.0){
		angle_of_attack = atan(INFINITY);
	}
	else{
		angle_of_attack = abs(atan(x));
	}
	Matrix3d rotate_y;
	rotate_y << 0, 0, -1, 0, 1, 0, 1, 0, 0;
	if(V_local[2] > 0.0){
		rotate_y(0, 2) = 1;
		rotate_y(2, 0) = -1;
	}
	double C_l = flapping_model::C_l0 * std::sin(2 * angle_of_attack);
	double C_d = flapping_model::D_l0 - flapping_model::D_l1 * std::cos(2 * angle_of_attack);
	double F_l = 2.0 / 3.0 * C_l * flapping_model::p * area * (V_local_proj.norm() * V_local_proj.norm());
	double F_d = 2.0 / 3.0 * C_d * flapping_model::p * area * (V_local_proj.norm() * V_local_proj.norm());
	// std::cout << "V: " << V_local_proj.norm() << std::endl;
	// std::cout << "F_l: " << F_l << std::endl;
	// std::cout << "F_d: " << F_d << std::endl;

	// body上的点在惯性系下的位置
	Vector3d pos = CalcBodyToBaseCoordinates(model, Q, body_id, body_pos, true);
	// std::cout << "pos: " << pos << std::endl;

	// 惯性系下的外力
	Vector3d F = F_l * MatWorld2Body.transpose() * (rotate_y * V_local_proj.normalized())
			- F_d * ((MatWorld2Body.transpose() * V_local_proj).normalized());
	// Vector3d F = F_l * (rotate_y * V_local_proj.normalized()) - F_d * V_local_proj.normalized();
	// std::cout << "1: " << F_d * (MatWorld2Body.transpose() * V_local_proj).normalized() << std::endl;
	// std::cout << "2: " << F_d * ((MatWorld2Body.transpose() * V_local_proj).normalized()) << std::endl;
	// std::cout << "F: " << F.transpose() << std::endl;

	// 惯性系下的力矩
	Vector3d T = pos.cross(F);
	// std::cout << "T: " << T << std::endl;
	fext[body_id][0] = T[0];
	fext[body_id][1] = T[1];
	fext[body_id][2] = T[2];
	fext[body_id][3] = F[0];
	fext[body_id][4] = F[1];
	fext[body_id][5] = F[2];
	// std::cout << "fext: " << fext[body_id].transpose() << std::endl;
}

void urdfRead (vector<adouble> x, adouble u, VectorNd& QDDot, adouble t) {
	rbdl_check_api_version (RBDL_API_VERSION);
	Model* model = new Model();
	if (!Addons::URDFReadFromFile ("/home/lj/new_model/src/urdf/model.urdf", model, true, false)) {
		std::cerr << "Error loading model " << std::endl;
		abort();
	}

	double time = t.value();
	double theta1 = 0.806 * std::sin(flapping_model::f * 2.0 * 3.1415 * time) + 0.241;
	double theta3 = -0.806 * std::sin(flapping_model::f * 2.0 * 3.1415 * time) - 0.241;
	double theta2 = 0.35 * std::sin(flapping_model::f * 2.0 * 3.1415 * time + 3.1415 / 2.0);
	double theta4 = theta2;
	double theta1_dot = flapping_model::f * 2.0 * 3.1415 * 0.806 * 
						std::cos(flapping_model::f * 2.0 * 3.1415 * time);
	double theta3_dot = -0.806 * flapping_model::f * 2.0 * 3.1415 * 
						std::cos(flapping_model::f * 2.0 * 3.1415 * time);
	double theta2_dot = flapping_model::f * 2.0 * 3.1415 * 0.35 *
						std::cos(flapping_model::f * 2.0 * 3.1415 * time + 3.1415 / 2.0);
	double theta4_dot = theta2_dot;
	VectorNd Q = VectorNd::Zero (model->q_size);
    Q << x[1].value(), 0.0, x[2].value(), 0.0, std::sin(x[0].value() / 2.0), 0.0, 
			theta1, theta2, theta3, theta4, x[3].value(), std::cos(x[0].value() / 2.0);
    // std::cout << "Q : " << Q.transpose() << std::endl;

	VectorNd QDot = VectorNd::Zero (model->qdot_size);
	QDot << x[5].value(), 0.0, x[6].value(), 0.0, x[4].value(), 0.0,
			 theta1_dot, theta2_dot, theta3_dot, theta4_dot, x[7].value();
    // std::cout << "QDOT : " << QDot.transpose() << std::endl;

	VectorNd Tau = VectorNd::Zero (model->qdot_size);
    Tau[10] = u.value();
	// std::cout << "tau: " << Tau[10] << std::endl;
	QDDot = VectorNd::Zero (model->qdot_size);
	std::vector<RigidBodyDynamics::Math::SpatialVector> fext;
	fext.resize(model->mBodies.size());
	for (unsigned int i = 0; i < fext.size(); ++i) {
		fext[i]=RigidBodyDynamics::Math::SpatialVector::Zero();
	}
	CalF(*model, Q, QDot, "left_wing", fext);
	CalF(*model, Q, QDot, "right_wing", fext);
	CalF(*model, Q, QDot, "tail", fext);
	// std::cout << "tail force: " << fext[model->GetBodyId("tail")].transpose() << std::endl;
	// std::cout << "right wing force: " << fext[model->GetBodyId("right_wing")].transpose() << std::endl;
	// std::cout << "left wing force: " << fext[model->GetBodyId("left_wing")].transpose() << std::endl;

	ForwardDynamics (*model, Q, QDot, Tau, QDDot, &fext);
	// ForwardDynamics (*model, Q, QDot, Tau, QDDot);
	// std::cout << "QDDot: " << QDDot.transpose() << std::endl;

	delete model;

}

