#include <iostream>

// Set the necessary macro variables before importing the file
// Make sure you set the variable DIM in the CMakeLkists.txt file
// otherwise it will not compile
#define NS_INPUT

#include "NavierStokes.hpp"


// Since we defined the macro variable NS_INPUT, the methods of inlet velocity
// are only declared, so we have to write their implementation in the main file.
static constexpr double U_m = 10.0;
static constexpr double H = 20.0;

void
NavierStokes::InletVelocity::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  //values[0] = 1.0;
  for (unsigned int i = 0; i < dim + 1; ++i)
    values[i] = 0.0;
  
  values[0] = 4 * 1.5 * p[1] * (0.41 - p[1]) / (0.41 * 0.41);
  
}

double
NavierStokes::InletVelocity::value(const Point<dim> &p, const unsigned int component) const
{
  if (/*component == 1 ||*/ component == 0)
  {
    if (/*p[1] < 0.24*/1)
      return 4 * 1.5 * p[1] * (0.41 - p[1]) / (0.41 * 0.41);
    else
      return 0.0;
  }
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

#if DIM == 3
  const std::string mesh_file_name("../mesh/domain3D.msh");
#else
  const std::string mesh_file_name("../mesh/domain2D.msh");
#endif

  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;


  NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, 0.01, 0.2, 2); // last numbers: deltat, T, step
  //problem.set_re_number(20);
  problem.setup();
  problem.compute_ordered_dofs_indices();
  problem.solve();
  return 0;
}
