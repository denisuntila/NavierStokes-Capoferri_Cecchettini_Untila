#include <iostream>
#include "NavierStokes.hpp"


int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

#ifdef DIM
  const std::string mesh_file_name("../mesh/domain3D.msh");
#else
  const std::string mesh_file_name("../mesh/domain2D.msh");
#endif

  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;


  NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, 0.01, 2.0, 2); // last numbers: deltat, T, step
  problem.set_re_number(20);
  problem.setup();
  problem.compute_ordered_dofs_indices();
  problem.solve(0);
  return 0;
}
