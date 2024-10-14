#include <iostream>
#include "NavierStokes.hpp"

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string mesh_file_name("../mesh/domain2D.msh");
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

  NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, 0.001, 0.003, 1);
  problem.setup();
  problem.compute_ordered_dofs_indices();
  //problem.solve();
  //problem.import_data();
  //problem.solve2();

  //problem.export_data_old();
  //problem.export_data();
  return 0;
}
