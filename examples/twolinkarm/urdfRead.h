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
	const double wing_area = 0.085; // 0.085 0.046
	const double tail_area = 0.027; //0.03845369   0.018  0.027
	const double C_l0 = 1.332;
	const double D_l0 = 1.713;
	const double D_l1 = -1.639;
	const Vector3d wing_pos{-0.0555, 0.0, 0.0}; //-0.03238, 0.0, 0.0
	const Vector3d tail_pos{-0.14187, 0.0, 0.0}; //-0.2635, 0, 0.0268   -0.10905, 0.0, 0.02466
	const double v0 = 6.0;
    const double target_height = 1.0;
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
	Matrix3d MatWorld2Body = CalcBodyWorldOrientation(model, Q, body_id, true);
	Eigen::Quaterniond q(MatWorld2Body);

	// 从惯性系到body(local)系下的旋转矩阵
	MatWorld2Body = q.normalized().toRotationMatrix();

	// body(local)系下的速度
	Vector3d V_local = MatWorld2Body * V_inertial;// CalcBodyWorldOrientation
	Vector3d V_local_proj{V_local[0], 0, V_local[2]};
	double x, angle_of_attack;
	if(V_local[0] == 0){
		angle_of_attack = atan(INFINITY);
	}
	else{
		x = V_local[2] / V_local[0];
	}
	angle_of_attack = atan(x);
	Matrix3d rotate_y;
	rotate_y << 0, 0, 1, 0, 1, 0, -1, 0, 0;
	double C_l = flapping_model::C_l0 * std::sin(2 * angle_of_attack);
	double C_d = flapping_model::D_l0 - flapping_model::D_l1 * std::cos(2 * angle_of_attack);
	double F_l = 2.0 / 3.0 * C_l * flapping_model::p * area * (V_local_proj.norm() * V_local_proj.norm());
	double F_d = 2.0 / 3.0 * C_d * flapping_model::p * area * (V_local_proj.norm() * V_local_proj.norm());

	// body上的点在惯性系下的位置
	Vector3d pos = CalcBodyToBaseCoordinates(model, Q, body_id, body_pos, true);

	// 惯性系下的外力
	Vector3d F = F_l * MatWorld2Body.transpose() * (rotate_y * V_local_proj.normalized())
			- F_d * ((MatWorld2Body.transpose() * V_local_proj).normalized());

	// 惯性系下的力矩
	Vector3d T = pos.cross(F);
	fext[body_id][0] = T[0];
	fext[body_id][1] = T[1];
	fext[body_id][2] = T[2];
	fext[body_id][3] = F[0];
	fext[body_id][4] = F[1];
	fext[body_id][5] = F[2];
}

void SetDynamicParam (vector<adouble> x, adouble u, VectorNd& QDDot, adouble t, Model* model) {

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

	VectorNd QDot = VectorNd::Zero (model->qdot_size);
	QDot << x[5].value(), 0.0, x[6].value(), 0.0, x[4].value(), 0.0,
			 theta1_dot, theta2_dot, theta3_dot, theta4_dot, x[7].value();

	VectorNd Tau = VectorNd::Zero (model->qdot_size);
    Tau[10] = u.value();
	QDDot = VectorNd::Zero (model->qdot_size);
	std::vector<RigidBodyDynamics::Math::SpatialVector> fext;
	fext.resize(model->mBodies.size());
	for (unsigned int i = 0; i < fext.size(); ++i) {
		fext[i]=RigidBodyDynamics::Math::SpatialVector::Zero();
	}
	CalF(*model, Q, QDot, "left_wing", fext);
	CalF(*model, Q, QDot, "right_wing", fext);
	CalF(*model, Q, QDot, "tail", fext);

	ForwardDynamics (*model, Q, QDot, Tau, QDDot, &fext);
}

