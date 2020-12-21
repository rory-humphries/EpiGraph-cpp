
#include <Eigen/Dense>
#include <EpiGraph/EpiGraph.hpp>
#include <EpiGraph/Models/SI.hpp>
#include <boost/numeric/odeint.hpp>
#include <iostream>

using namespace Eigen;
using namespace EpiGraph;
using namespace boost::numeric::odeint;

int main(int argc, char *argv[]) {
  print_banner2();

  int run_time = 1000;
  int time_step = 1;

  Vector2d x = {5000, 1};              // S, I
  Vector4d p = {0.3, 0.1, 0.01, 0.01}; // infection, recovery, birth, death

  runge_kutta4<Vector2d, double, Vector2d, double, vector_space_algebra> rk;

  auto sys = [p](const Vector2d &x, Vector2d &dxdt, double t) {
    si_ode(x, dxdt, p);
  };

  for (int t = 0; t < run_time; t += time_step) {
    rk.do_step(sys, x, t, time_step);
    std::cout << x.transpose() << std::endl;
  }
}
