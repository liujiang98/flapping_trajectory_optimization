//////////////////////////////////////////////////////////////////////////
//////////////////        twolinkarm.cxx        //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:                 Two link arm problem      ////////////////
//////// Last modified:         04 January 2009           ////////////////
//////// Reference:             PROPT users guide         ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"
#include "urdfRead.h"
int i = 0;
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
    return  0.5*u * u;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    adouble xdot, ydot, vdot;
    vector<adouble> x{states[0], states[1], states[2], 
                    states[3], states[4], states[5], states[6], states[7]};

    adouble u = controls[0];
    // std::cout << u.value() << std::endl;
    VectorNd QDDot;
    adouble t = time;
    urdfRead(x, u, QDDot, t);
    ++i;
    // std::cout << ++i << std::endl;
    derivatives[ 0 ] = x[4];
    derivatives[ 1 ] = x[5];
    derivatives[ 2 ] = x[6];
    derivatives[ 3 ] = x[7];
    derivatives[ 4 ] = QDDot[4];
    derivatives[ 5 ] = QDDot[0];
    derivatives[ 6 ] = QDDot[2];
    derivatives[ 7 ] = QDDot[10];
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

   e[ 0 ] = x20;
   e[ 1 ] = x30;
   e[ 2 ] = x40;
   e[ 3 ] = x50;
   e[ 4 ] = x60 - flapping_model::v0 * cos(x10);
   e[ 5 ] = x70 - flapping_model::v0 * sin(x10);
   e[ 6 ] = x80;
   e[ 7 ] = pzf - flapping_model::target_height;
   e[ 8 ] = pzf_dot;

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

    problem.name        					= "Two link robotic arm";

    problem.outfilename                 	= "twolink.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   						= 1;
    problem.nlinkages                   	= 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   				= 8;
    problem.phases(1).ncontrols 				= 1;
    problem.phases(1).nevents   				= 9;
    problem.phases(1).npath     				= 0;
    problem.phases(1).nodes                     << 40;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;
    DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.phases(1).bounds.lower.states(0) = -3.14;
    problem.phases(1).bounds.lower.states(1) = 0.0;
    problem.phases(1).bounds.lower.states(2) = 0.0;
    problem.phases(1).bounds.lower.states(3) = -2.0;
    problem.phases(1).bounds.lower.states(4) = -10.0;
    problem.phases(1).bounds.lower.states(5) = -10.0;
    problem.phases(1).bounds.lower.states(6) = -10.0;
    problem.phases(1).bounds.lower.states(7) = -10.0;

    problem.phases(1).bounds.upper.states(0) = 3.14;
    problem.phases(1).bounds.upper.states(1) = 12.0;
    problem.phases(1).bounds.upper.states(2) = 12.0;
    problem.phases(1).bounds.upper.states(3) = 2.0;
    problem.phases(1).bounds.upper.states(4) = 10.0;
    problem.phases(1).bounds.upper.states(5) = 10.0;
    problem.phases(1).bounds.upper.states(6) = 10.0;
    problem.phases(1).bounds.upper.states(7) = 10.0;

    problem.phases(1).bounds.lower.controls(0) = -5.0;
    problem.phases(1).bounds.upper.controls(0) = 5.0;

    problem.phases(1).bounds.lower.events(0) = 0.0;
    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.lower.events(2) = 0.0;
    problem.phases(1).bounds.lower.events(3) = 0.0;
    problem.phases(1).bounds.lower.events(4) = 0.0;
    problem.phases(1).bounds.lower.events(5) = 0.0;
    problem.phases(1).bounds.lower.events(6) = 0.0;
    problem.phases(1).bounds.lower.events(7) = -0.5;
    problem.phases(1).bounds.lower.events(8) = -0.5;

    problem.phases(1).bounds.upper.events(0) = 0.0;
    problem.phases(1).bounds.upper.events(1) = 0.0;
    problem.phases(1).bounds.upper.events(2) = 0.0;
    problem.phases(1).bounds.upper.events(3) = 0.0;
    problem.phases(1).bounds.upper.events(4) = 0.0;
    problem.phases(1).bounds.upper.events(5) = 0.0;
    problem.phases(1).bounds.upper.events(6) = 0.0;
    problem.phases(1).bounds.upper.events(7) = 0.5;
    problem.phases(1).bounds.upper.events(8) = 0.5;



    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.0;
    problem.phases(1).bounds.upper.EndTime      = 10.0;


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


    MatrixXd x0(8,40);

    x0 <<  linspace(0.0, 0.0, 40),
           linspace(0.0, 0.0, 40),
           linspace(0.0, 0.0, 40),
           linspace(0.0, 0.0, 40),
           linspace(0.0, 0.0, 40),
           linspace(0.0, flapping_model::v0, 40),
           linspace(0.0, flapping_model::v0, 40),
           linspace(0.0, 0.0, 40);

    problem.phases(1).guess.controls       = zeros(1,40);
    problem.phases(1).guess.states         = x0;
    problem.phases(1).guess.time           = linspace(0.0, 10.0, 40);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.collocation_method          = "Hermite-Simpson";

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
    std::cout << "///////////////////////////use dae function " << i << std::endl;
////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4 x5 x6 x7 x8");

    plot(t,u,problem.name + ": controls", "time (s)", "controls", "u");


    plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4",
                                  "pdf", "twolinkarm_states.pdf");

    plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1 u2",
                              "pdf", "twolinkarm_controls.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

