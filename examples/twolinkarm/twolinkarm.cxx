#include "psopt.h"
#include "urdfRead.h"
#include <cmath>

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    adouble u = controls[0];
    return  u * u;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    vector<adouble> x{states[0], states[1], states[2], 
                    states[3], states[4], states[5], states[6], states[7]};
    // std::cout << "t: " << time << std::endl;
    // std::cout << "states: " << states[0] << " " << states[1] << " " << states[2] 
    //             << " " << states[3] << " " << states[4] << " " << states[5]
    //             << " " << states[6] << " " << states[7] << std::endl;
    adouble u = controls[0];
    // std::cout << u.value() << std::endl;
    VectorNd QDDot;
    adouble t = time;
    urdfRead(x, u, QDDot, t);
    // std::cout << ++i << std::endl;
    derivatives[ 0 ] = x[4];
    derivatives[ 1 ] = x[5];
    derivatives[ 2 ] = x[6];
    derivatives[ 3 ] = x[7];
    derivatives[ 4 ] = QDDot[4];
    derivatives[ 5 ] = QDDot[0];
    derivatives[ 6 ] = QDDot[2];
    derivatives[ 7 ] = QDDot[10];
    // std::cout << "dea: " << QDDot[4] << " " << derivatives[4] << std::endl;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   adouble x10 = initial_states[ 0 ];
   adouble x20 = initial_states[ 1 ];
   adouble x30 = initial_states[ 2 ];
   adouble x40 = initial_states[ 3 ];
   adouble x50 = initial_states[ 4 ];
   adouble x60 = initial_states[ 5 ];
   adouble x70 = initial_states[ 6 ];
   adouble x80 = initial_states[ 7 ];
   adouble pzf = final_states[ 2 ];
   adouble pzf_dot = final_states[ 6 ];
    //    std::cout << "states: " << initial_states[0] << " " << initial_states[1] << " " << initial_states[2] 
    //                 << " " << initial_states[3] << " " << initial_states[4] << " " << initial_states[5]
    //                 << " " << initial_states[6] << " " << initial_states[7] << std::endl;
   e[ 0 ] = x10;
   e[ 1 ] = x20;
   e[ 2 ] = x30;
   e[ 3 ] = x40;
   e[ 4 ] = x50;
   e[ 5 ] = x60 - flapping_model::v0 * cos(x10);
   e[ 6 ] = x70 + flapping_model::v0 * sin(x10);
   e[ 7 ] = x80;
   e[ 8 ] = pzf;
   e[ 9 ] = pzf_dot;
//    std::cout << "cos(3.14): " << cos(x) << std::endl;
//    std::cout << "test: " << x10 << " " << cos(x10) << " " << 
//                 flapping_model::v0 * cos(x10) << " "
//                 <<  x60 - flapping_model::v0 * cos(x10) << " "<< e[4] << std::endl;
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////


int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        					= "flapping wing robot";

    problem.outfilename                 	= "flapping.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   						= 1;
    problem.nlinkages                   	= 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////
    int nodes_num = 300;
    problem.phases(1).nstates   				= 8;
    problem.phases(1).ncontrols 				= 1;
    problem.phases(1).nevents   				= 10;
    problem.phases(1).npath     				= 0;
    problem.phases(1).nodes                     << nodes_num;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;
    DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.phases(1).bounds.lower.states(0) = -PSOPT::pi / 2.0;
    problem.phases(1).bounds.lower.states(1) = 0.0;
    problem.phases(1).bounds.lower.states(2) = 0.0;
    problem.phases(1).bounds.lower.states(3) = -1.547;
    problem.phases(1).bounds.lower.states(4) = -6.0;
    problem.phases(1).bounds.lower.states(5) = 0.0;
    problem.phases(1).bounds.lower.states(6) = -5.0;
    problem.phases(1).bounds.lower.states(7) = -50.0;

    problem.phases(1).bounds.upper.states(0) = PSOPT::pi / 2.0;
    problem.phases(1).bounds.upper.states(1) = 20;
    problem.phases(1).bounds.upper.states(2) = 2 * flapping_model::target_height;
    problem.phases(1).bounds.upper.states(3) = 1.547;
    problem.phases(1).bounds.upper.states(4) = 6.0;
    problem.phases(1).bounds.upper.states(5) = flapping_model::v0;
    problem.phases(1).bounds.upper.states(6) = 5.0;
    problem.phases(1).bounds.upper.states(7) = 50.0;

    problem.phases(1).bounds.lower.controls(0) = -10.0;
    problem.phases(1).bounds.upper.controls(0) = 10.0;

    problem.phases(1).bounds.lower.events(0) = -PSOPT::pi / 2.0;
    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.lower.events(2) = 0.0;
    problem.phases(1).bounds.lower.events(3) = -0.0;
    problem.phases(1).bounds.lower.events(4) = 0.0;
    problem.phases(1).bounds.lower.events(5) = -0.0;
    problem.phases(1).bounds.lower.events(6) = -0.0;
    problem.phases(1).bounds.lower.events(7) = 0.0;
    problem.phases(1).bounds.lower.events(8) = flapping_model::target_height - 0.01;
    problem.phases(1).bounds.lower.events(9) = -0.01;
    
    problem.phases(1).bounds.upper.events(0) = 0.0;
    problem.phases(1).bounds.upper.events(1) = 0.0;
    problem.phases(1).bounds.upper.events(2) = 0.0;
    problem.phases(1).bounds.upper.events(3) = -0.0;
    problem.phases(1).bounds.upper.events(4) = 0.0;
    problem.phases(1).bounds.upper.events(5) = 0.0;
    problem.phases(1).bounds.upper.events(6) = 0.0;
    problem.phases(1).bounds.upper.events(7) = 0.0;
    problem.phases(1).bounds.upper.events(8) = flapping_model::target_height + 0.01;
    problem.phases(1).bounds.upper.events(9) = 0.01;



    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.0;
    problem.phases(1).bounds.upper.EndTime      = 3.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 				    = &integrand_cost;
    problem.endpoint_cost 					= &endpoint_cost;
    problem.dae 							= &dae;
    problem.events 							= &events;
    problem.linkages						= &linkages;



////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes = problem.phases(1).nodes(0);

    MatrixXd x0(8,nnodes);

    x0 <<  linspace(-0.35, -0.35, nnodes),
           linspace(0.0, 2.0, nnodes),
           linspace(0.0, flapping_model::target_height, nnodes),
           linspace(-0.0, 1.0, nnodes),
           linspace(0.0, 0.0, nnodes),
           linspace(0.5 * flapping_model::v0, 0.5 * flapping_model::v0, nnodes),
           linspace(0.5 * flapping_model::v0, 0.0, nnodes),
           linspace(0.0, 0.0, nnodes);

    problem.phases(1).guess.controls       = linspace(-0.0, 0.0, nnodes);
    problem.phases(1).guess.states         = x0;
    problem.phases(1).guess.time           = linspace(0.0, 3.0, nnodes);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    // algorithm.derivatives                 = "automatic";
    algorithm.derivatives                 = "numerical";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-2;
    algorithm.ode_tolerance               = 1.e-2;
    algorithm.collocation_method          = "Hermite-Simpson";
    // algorithm.mesh_refinement             = "automatic";
    algorithm.ipopt_max_cpu_time          = 72000.0;
    // algorithm.defect_scaling              = "jacobian-based";

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x 		= solution.get_states_in_phase(1);
    u 		= solution.get_controls_in_phase(1);
    t 		= solution.get_time_in_phase(1);
////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd new_x(3, x.cols());
    new_x << x.row(0), 
             x.row(2),
             x.row(3);

    plot(t, x.row(1), problem.name + ": states", "time (s)", "states", "x2");
    plot(t, new_x, problem.name + ": states", "time (s)", "states", "x1 x3 x4");
    plot(t, x.block(4, 0, 4, x.cols()), problem.name + ": states", "time (s)",
                                                         "states", "x5 x6 x7 x8");
    plot(t, u, problem.name + ": controls", "time (s)", "controls", "u");
    plot(t, x, problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4 x5 x6 x7 x8",
                                  "pdf", "flappingwing_states.pdf");
    plot(t, u, problem.name + ": controls", "time (s)", "controls", "u",
                              "pdf", "flappingwing_controls.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

