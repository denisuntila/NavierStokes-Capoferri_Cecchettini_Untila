#include <iostream>

// TEST 01 3D

// Set the necessary macro variables before importing the file
// Make sure you set the variable DIM in the CMakeLkists.txt file
// otherwise it will not compile
#define NS_INPUT

#include "NavierStokes.hpp"


// Since we defined the macro variable NS_INPUT, the methods of inlet velocity
// are only declared, so we have to write their implementation in the main file.
static constexpr double U_m = 0.45;
static constexpr double H = 0.41;

void
NavierStokes::InletVelocity::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  for (unsigned int i = 0; i < dim + 1; ++i)
    values[i] = 0.0;
  
  values[0] = 16 * U_m * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);
  
}

double
NavierStokes::InletVelocity::value(const Point<dim> &p, const unsigned int component) const
{
  if (component == 0)
    return 16 * U_m * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);
  else
    return 0.0;
}


double
NavierStokes::InletVelocity::get_mean_vel()
{
  return 4.0 * U_m / 9.0;
}




int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string mesh_file_name("../../../../mesh/domain3D.msh");

  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;


  NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, 0.01, 2.0, 10); // last numbers: deltat, T, step
  problem.set_re_number(20);
  problem.setup();
  problem.compute_ordered_dofs_indices();
  problem.solve();
  return 0;
}
