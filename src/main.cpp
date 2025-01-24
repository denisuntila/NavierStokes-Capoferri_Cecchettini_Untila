#include <iostream>
#define NS_INPUT
#include "NavierStokes.hpp"



void
NavierStokes::InletVelocity::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  //values[0] = 1.0;
  for (unsigned int i = 0; i < dim + 1; ++i)
    values[i] = 0.0;
  if (/*p[1] < 0.24*/ 1)
    values[0] = 4 * 1.5 * p[1] * (0.41 - p[1]) / (0.41 * 0.41);
  
  //values[0] = 10.0;
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

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string mesh_file_name("../mesh/domain2D.msh");
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;


  NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, 0.01, 5.6, 2);
  problem.setup();
  problem.compute_ordered_dofs_indices();
  problem.solve(550);
  return 0;
}
