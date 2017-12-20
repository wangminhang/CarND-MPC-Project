#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 15;
double dt = 0.05;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

double ref_cte = 0;
double ref_epsi = 0;
double ref_v = 50;
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;
const double w_state_cte = 300.0;
const double w_state_epsi = 100.0;
const double w_state_v = 1.0;
const double w_val_steering = 200.0;
const double w_val_throttle = 50.0;
const double w_seq_steering = 5000.0;
const double w_seq_throttle = 100.0;
const double max_angle = 25 * M_PI / 180;

class FG_eval
{
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars)
  {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    fg[0] = 0;

    for (size_t m = 0; m < N; m++)
    {
      fg[0] += w_state_cte * CppAD::pow(vars[cte_start + m] - ref_cte, 2);
      fg[0] += w_state_epsi * CppAD::pow(vars[epsi_start + m] - ref_epsi, 2);
      fg[0] += w_state_v * CppAD::pow(vars[v_start + m] - ref_v, 2);
    }

    for (size_t m = 0; m < N - 1; m++)
    {
      fg[0] += w_val_steering * CppAD::pow(vars[delta_start + m], 2);
      fg[0] += w_val_throttle * CppAD::pow(vars[a_start + m], 2);
    }

    for (size_t m = 0; m < N - 2; m++)
    {
      fg[0] += w_seq_steering * CppAD::pow(vars[delta_start + m + 1] - vars[delta_start + m], 2);
      fg[0] += w_seq_throttle * CppAD::pow(vars[a_start + m + 1] - vars[a_start + m], 2);
    }

    fg[x_start + 1] = vars[x_start];
    fg[y_start + 1] = vars[y_start];
    fg[psi_start + 1] = vars[psi_start];
    fg[v_start + 1] = vars[v_start];
    fg[cte_start + 1] = vars[cte_start];
    fg[epsi_start + 1] = vars[epsi_start];

    for (size_t m = 0; m < N - 1; m++)
    {
      AD<double> x1 = vars[x_start + m + 1];
      AD<double> y1 = vars[y_start + m + 1];
      AD<double> psi1 = vars[psi_start + m + 1];
      AD<double> v1 = vars[v_start + m + 1];
      AD<double> cte1 = vars[cte_start + m + 1];
      AD<double> epsi1 = vars[epsi_start + m + 1];

      AD<double> x0 = vars[x_start + m];
      AD<double> y0 = vars[y_start + m];
      AD<double> psi0 = vars[psi_start + m];
      AD<double> v0 = vars[v_start + m];
      AD<double> cte0 = vars[cte_start + m];
      AD<double> epsi0 = vars[epsi_start + m];

      AD<double> delta0 = vars[delta_start + m];
      AD<double> a0 = vars[a_start + m];

      // AD<double> f0 = coeffs[0] + coeffs[1] * x0;
      // AD<double> psides0 = CppAD::atan(coeffs[1]);
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0;
      AD<double> f0_deriv = coeffs[1]+ 2.0 * coeffs[2] * x0;
      AD<double> psides0 = CppAD::atan(f0_deriv);

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.

      // TODO: Setup the rest of the model constraints
      fg[x_start + m + 2] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[y_start + m + 2] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[psi_start + m + 2] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[v_start + m + 2] = v1 - (v0 + a0 * dt);
      fg[cte_start + m + 2] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs)
{
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0]; // posX
  double y = state[1]; // posY
  double psi = state[2]; // orientation
  double v = state[3]; // velocity
  double cte = state[4]; // cross track error cte
  double epsi = state[5]; // steering angle error epsi

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = N * 6 + (N - 1) * 2;
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int m = 0; m < n_vars; m++)
  {
    vars[m] = 0;
  }

  // initial values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  for (size_t m = 0; m < delta_start; m++)
  {
    vars_lowerbound[m] = -1 * numeric_limits<double>::max();
    vars_upperbound[m] = numeric_limits<double>::max();
  }

  for (size_t m = delta_start; m < a_start; m++)
  {
    vars_lowerbound[m] = -max_angle;
    vars_upperbound[m] = max_angle;
  }

  // Acceleration & braking
  for (size_t m = a_start; m < n_vars; m++)
  {
    vars_lowerbound[m] = -0.5;
    vars_upperbound[m] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);

  for (int m = 0; m < n_constraints; m++)
  {
    constraints_lowerbound[m] = 0;
    constraints_upperbound[m] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;
  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  this->mpc_x = {};
  this->mpc_y = {};

  for (size_t m = 0; m < N; m++)
  {
    this->mpc_x.push_back(solution.x[x_start+m]);
    this->mpc_y.push_back(solution.x[y_start+m]);
  }

  return {
    solution.x[x_start + 1],   solution.x[y_start + 1],
    solution.x[psi_start + 1], solution.x[v_start + 1],
    solution.x[cte_start + 1], solution.x[epsi_start + 1],
    solution.x[delta_start],   solution.x[a_start]
  };
}
