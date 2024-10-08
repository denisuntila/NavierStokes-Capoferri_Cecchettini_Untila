#include "NavierStokes.hpp"


void
NavierStokes::setup()
{
  // Import the mesh
  {
    pcout << "Importing the mesh from " << mesh_file_name << std::endl;

    Triangulation<dim> mesh_serial; 
    GridIn<dim> grid_in;

    grid_in.attach_triangulation(mesh_serial);

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::
      Utilities::create_description_from_triangulation(
      mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "\tNumber of elements = " << mesh.n_global_active_cells()
      << std::endl;
  }

  pcout << "---------------------------------------------------" << std::endl;

  // Initialize the finite element spaces for velocity and pressure
  {
    pcout << "Initializing the finite element spaces" << std::endl;

    const FE_Type<dim> fe_scalar_pressure(degree_pressure);
    const FE_Type<dim> fe_scalar_velocity(degree_velocity);

    fe = std::make_unique<FESystem<dim>>(
      fe_scalar_velocity, dim,
      fe_scalar_pressure, 1
    );

    pcout << "\tVelocity degree:\t\t" << fe_scalar_velocity.degree
      << std::endl; 
    pcout << "\tPressure degree:\t\t" << fe_scalar_pressure.degree
      << std::endl;
    pcout << "\tDoFs per cell:\t\t\t" << fe->dofs_per_cell
      << std::endl;

    quadrature = std::make_unique<Q_Type<dim>>(fe->degree + 1);
    pcout << "\tQuadrature points per cell:\t" << quadrature->size()
      << std::endl;

    quadrature_face = std::make_unique<Q_Type<dim - 1>>(fe->degree + 1);
    pcout << "\tQuadrature points per face:\t" << quadrature_face->size()
      << std::endl;
  }

  pcout << "---------------------------------------------------" << std::endl;

  // Initialize the DoF handler
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    block_owned_dofs.resize(2);
    block_relevant_dofs.resize(2);
    
    block_owned_dofs[0] = locally_owned_dofs.get_view(0, n_u);
    block_owned_dofs[1] = locally_owned_dofs.get_view(n_u, n_u + n_p);

    block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    pcout << "\tNumber of DoFs:" << std::endl;
    pcout << "\t\tvelocity:\t\t" << n_u << std::endl;
    pcout << "\t\tpressure:\t\t" << n_p << std::endl;
    pcout << "\t\ttotal:\t\t\t" << n_u + n_p << std::endl; 
  }

  pcout << "---------------------------------------------------" << std::endl;
  
  // Initialize the linear system
  {
    pcout << "Initializing the linear system" << std::endl;
    pcout << "\tInitializing the sparsity pattern" << std::endl;

    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    for (unsigned int c = 0; c < dim + 1; ++c)
    {
      for (unsigned int d = 0; d < dim + 1; ++d)
      {
        if (c == dim && d == dim)
          coupling[c][d] = DoFTools::none;
        else
          coupling[c][d] = DoFTools::always;
      }
    }

    TrilinosWrappers::BlockSparsityPattern sparsity(
      block_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
      coupling, sparsity);
    sparsity.compress();

    pcout << "\tInitializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    deltat_lumped_mass_inv.reinit(block_owned_dofs, MPI_COMM_WORLD);

    pcout << "\tInitializing the system right-hand side" << std::endl;
    system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "\tInitializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);

  }
  
}

void
NavierStokes::assemble(const double &time)
{

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();
  const unsigned int n_q_face      = quadrature_face->size();

  FEValues<dim> fe_values(*fe, *quadrature,
    update_values | update_gradients |
    update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_face_values(*fe, *quadrature_face,
    update_values | update_normal_vectors |
    update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_lumped_mass_inv(dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs    = 0.0;
  deltat_lumped_mass_inv = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  std::vector<Tensor<1, dim>> current_velocity_values(n_q);
  forcing_term.set_time(time);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    fe_values.reinit(cell);

    cell_matrix               = 0.0;
    cell_lumped_mass_inv      = 0.0;
    cell_rhs                  = 0.0;

    fe_values[velocity].get_function_values(solution, current_velocity_values);

    for (unsigned int q = 0; q < n_q; ++q)
    {
      Vector<double> forcing_term_loc(dim);
      forcing_term.vector_value(fe_values.quadrature_point(q),
        forcing_term_loc);
      Tensor<1, dim> forcing_term_tensor;
      for (unsigned int d = 0; d < dim; ++d)
        forcing_term_tensor[d] = forcing_term_loc[d];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        double temp = 0;
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          cell_matrix(i, j) +=
            scalar_product(fe_values[velocity].value(i, q),
              fe_values[velocity].value(j, q)) *
            fe_values.JxW(q) / deltat;

          // Viscosity term.
          cell_matrix(i, j) += nu *
            scalar_product(fe_values[velocity].gradient(i, q),
              fe_values[velocity].gradient(j, q)) *
            fe_values.JxW(q);

          // T1
          
          cell_matrix(i, j) += 
            fe_values[velocity].value(i, q) *
            fe_values[velocity].gradient(j, q) *
            current_velocity_values[q] *
            fe_values.JxW(q);
          

          // T2
          /*
          cell_matrix(i, j) += 
            current_velocity_values[q] * fe_values[velocity].gradient(j, q) *
            fe_values[velocity].value(i, q) *
            fe_values.JxW(q);
          */
          


          // Pressure term in the momentum equation.
          cell_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
            fe_values[pressure].value(j, q) *
            fe_values.JxW(q);

          // Pressure term in the continuity equation.
          cell_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
            fe_values[pressure].value(i, q) *
            fe_values.JxW(q);

          // Lumped mass matrix.
          temp += std::abs(
            fe_values[velocity].value(j, q) *
            fe_values[velocity].value(i, q) *
            fe_values.JxW(q)
          );
          
        }

        // Forcing term.
        cell_rhs(i) += scalar_product(forcing_term_tensor,
          fe_values[velocity].value(i, q)) *
          fe_values.JxW(q);
        
        cell_rhs(i) +=
          scalar_product(current_velocity_values[q],
            fe_values[velocity].value(i, q)) *
          fe_values.JxW(q) / deltat;
        
        // Used for aYosida preconditioner.
        //cell_lumped_mass_inv(i) = deltat / temp;
        cell_lumped_mass_inv(i) += temp;
      }
    }

    // Boundary integral for Neumann BCs.
    if (cell->at_boundary())
    {
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        if (cell->face(f)->at_boundary() &&
            cell->face(f)->boundary_id() == 1)
        {
          fe_face_values.reinit(cell, f);

          for (unsigned int q = 0; q < n_q_face; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) += -p_out *
                scalar_product(fe_face_values.normal_vector(q),
                  fe_face_values[velocity].value(i, q)) *
                fe_face_values.JxW(q);
            }
          }
        }
      }
    }

    cell->get_dof_indices(dof_indices);

    system_matrix.add(dof_indices, cell_matrix);
    system_rhs.add(dof_indices, cell_rhs);
    deltat_lumped_mass_inv.add(dof_indices, cell_lumped_mass_inv);
  }

  for (const auto &i : locally_owned_dofs)
  {
    deltat_lumped_mass_inv[i] = deltat / deltat_lumped_mass_inv[i];
  }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
  deltat_lumped_mass_inv.compress(VectorOperation::add);

  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    std::vector<bool> c_mask(dim, true);
    c_mask.push_back(false);

    // We interpolate first the inlet velocity condition alone, then the wall
    // condition alone, so that the latter "win" over the former where the two
    // boundaries touch.
    boundary_functions[3] = &inlet_velocity;
    VectorTools::interpolate_boundary_values(
      dof_handler, boundary_functions, boundary_values,
      ComponentMask(c_mask)
    );

    boundary_functions.clear();
    Functions::ZeroFunction<dim> zero_function(dim + 1);
    boundary_functions[0] = &zero_function;
    boundary_functions[2] = &zero_function;
    boundary_functions[4] = &zero_function;

    VectorTools::interpolate_boundary_values(
      dof_handler, boundary_functions, boundary_values,
      ComponentMask(c_mask)
    );

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution_owned, system_rhs, false
    );
  }
}


void
NavierStokes::solve_time_step()
{
  SolverControl solver_control(500, 1e-6 * system_rhs.l2_norm());

  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  //PreconditionIdentity preconditioner;

  /*
  PreconditionASIMPLE preconditioner;
  preconditioner.initialize(
    system_matrix.block(0, 0),
    system_matrix.block(1, 0),
    system_matrix.block(0, 1),
    solution_owned
  );
  */
  
  
  PreconditionAYosida preconditioner;
  preconditioner.initialize(
    system_matrix.block(0, 0),
    system_matrix.block(1, 0),
    system_matrix.block(0, 1),
    deltat_lumped_mass_inv,
    solution_owned
  );
  
  

  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
  pcout << "  " << solver_control.last_step() << " GMRES iterations" << std::endl;

  solution = solution_owned;

}


void
NavierStokes::output(const unsigned int &time_step) const
{
  TrilinosWrappers::MPI::BlockVector solution_ghost(block_owned_dofs,
    block_relevant_dofs,
    MPI_COMM_WORLD);
  
  solution_ghost = solution;       

  DataOut<dim> data_out;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(
    DataComponentInterpretation::component_is_scalar);
  std::vector<std::string> names(dim, "velocity");
  names.push_back("pressure");

  data_out.add_data_vector(dof_handler,
    solution_ghost,
    names,
    data_component_interpretation);

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-stokes";
  data_out.write_vtu_with_pvtu_record("./",
    output_file_name,
    time_step,
    MPI_COMM_WORLD);
}


void
NavierStokes::solve()
{
  //assemble_matrices();

  pcout << "===================================================" << std::endl;

  time = 0.0;

  // Initial condition
  {
    pcout << "Applying initial conditions" << std::endl;

    VectorTools::interpolate(dof_handler, initial_conditions, solution_owned);
    solution = solution_owned;

    output(0);

    pcout << "---------------------------------------------------" << std::endl;    
  }

  unsigned int time_step = 0;

  while (time < T - 0.5 * deltat)
  {
    time += deltat;
    ++time_step;

    pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
      << time << ":" << std::flush;
    
    assemble(time);
    solve_time_step();
    if (time_step % step == 0)
      output(time_step);
  }

}

void
NavierStokes::export_data()
{
  // For now it works in singlecore
  // I'll try to fix it using the ghost of the vector

  TrilinosWrappers::MPI::BlockVector solution_ghost(block_owned_dofs,
    block_relevant_dofs,
    MPI_COMM_WORLD);
  solution_ghost = solution;
  
  unsigned int local_size = solution.locally_owned_size();
  //std::cout << local_size << std::endl;
  unsigned int max_local_size = 0;
  MPI_Allreduce(&local_size, &max_local_size,
    1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  
  std::vector<int> local_indices(max_local_size);
  for (unsigned int i = 0; i < max_local_size; ++i)
  {
    local_indices.at(i) = ((i < local_size) ? 
      (locally_owned_dofs.nth_index_in_set(i)) : -1);
  }

  std::unique_ptr<int[]> temp;
  std::unique_ptr<double[]> rbuf;

  if (0 == mpi_rank)
  {
    temp = std::make_unique<int[]>(mpi_size * max_local_size);
    rbuf = std::make_unique<double[]>(solution.size());
  }

  for (unsigned int i = 0; i < max_local_size; ++i)
  {
    MPI_Gather(&local_indices.at(i), 1, MPI_UNSIGNED, temp.get() +
      i * mpi_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for (unsigned int i = 0; i < max_local_size; ++i)
  {
    double local_data = 0;
    int scount = 0;
    
    if (i < local_size)
    {
      local_data = solution_ghost[locally_owned_dofs.nth_index_in_set(i)];
      scount = 1;
    }

    if (0 == mpi_rank)
    {
      int *displs = temp.get() + i * mpi_size;
      std::unique_ptr<int[]> rcount = std::make_unique<int[]>(mpi_size);

      for (unsigned int j = 0; j < mpi_size; ++j)
      {
        rcount.get()[j] = ((-1 != displs[j]) ? 1 : 0);
      }
      MPI_Gatherv(&local_data, scount, MPI_DOUBLE, rbuf.get(),
        rcount.get(), displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else if (1 == scount)
      MPI_Gatherv(&local_data, scount, MPI_DOUBLE, nullptr, 
        nullptr, nullptr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0)
  {
    std::ofstream output_file("test4.bin", std::fstream::binary);
    output_file.write((char *)rbuf.get(), solution.size() * sizeof(double));
    output_file.close();
  }
}


void
NavierStokes::solve2()
{
  //assemble_matrices();

  pcout << "===================================================" << std::endl;

  time = 0.0;

  // Initial condition
  {
    pcout << "Applying initial conditions" << std::endl;

    //VectorTools::interpolate(dof_handler, initial_conditions, solution_owned);
    solution = solution_owned;

    output(0);

    pcout << "---------------------------------------------------" << std::endl;    
  }

  unsigned int time_step = 0;

  while (time < T - 0.5 * deltat)
  {
    time += deltat;
    ++time_step;

    pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
      << time << ":" << std::flush;
    
    assemble(time);
    solve_time_step();
    if (time_step % step == 0)
      output(time_step);
  }

}

void
NavierStokes::import_data()
{
  // For now it works only if we import a vector obtained by
  // executing the solver with the same number of MPI processes
  /*
  TrilinosWrappers::MPI::BlockVector solution_ghost(block_owned_dofs,
  block_relevant_dofs,
  MPI_COMM_WORLD);
  */
  //solution_ghost = solution;

  std::ifstream infile("test4.bin", std::fstream::binary);
  std::vector<unsigned char> inbuff(std::istreambuf_iterator<char>(infile), {});
  infile.close();
  
  
  for (const auto &i : locally_owned_dofs)
  {
    solution_owned[i] = *(double*)&inbuff[i * sizeof(double)];
  }
  //std::cout << "Buffer size = " << inbuff.size() << " Bytes" << std::endl; 
  //solution_owned = solution_ghost;
}
