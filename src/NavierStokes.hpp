#ifndef HH_NS_HH
#define HH_NS_HH

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>
#include <mpi.h>



using namespace dealii;


template<size_t Dim>
using FE_Type = FE_SimplexP<Dim>;

template<size_t Dim>
using Q_Type = QGaussSimplex<Dim>;

class NavierStokes
{
public:
  static constexpr unsigned int dim = 2;

  class ForcingTerm : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &/*p*/,
      const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity()
      : Function<dim>(dim + 1)
    {}

    virtual void
    vector_value(const Point<dim> &/*p*/, Vector<double> &values) const override
    {
      //values[0] = 1.0;
      for (unsigned int i = 0; i < dim + 1; ++i)
        values[i] = 0.0;
      if (/*p[1] < 0.24*/ 1)
        values[0] = 3.0;
      
      //values[0] = 10.0;
    }

    virtual double
    value(const Point<dim> &/*p*/, const unsigned int component = 0) const override
    {
      if (/*component == 1 ||*/ component == 0)
      {
        if (/*p[1] < 0.24*/1)
          return 3.0;
        else
          return 0.0;
      }
      else
        return 0.0;
    }
  };

  class FunctionH : public Function<dim>
  {
  public:
    // Constructor.
    FunctionH()
    {
    }

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int /*component*/) const override
    {
      return 0.;
    }
  };

  class InitialConditions : public Function<dim>
  {
  public:
    InitialConditions()
      : Function<dim>(dim + 1)
    {}

    virtual void
    vector_value(const Point<dim> &/*p*/, Vector<double> &values) const override
    {
      //values[0] = 1.0;
      for (unsigned int i = 0; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> &/*p*/, const unsigned int component = 0) const override
    {
      if (component == dim)
        return 0.0;
      else
        return 0.0;
    }
  };

  class PreconditionIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(TrilinosWrappers::MPI::BlockVector       &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      dst = src;
    }

  protected:
  };


  
  class PreconditionASIMPLE
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void initialize(const TrilinosWrappers::SparseMatrix &F_,
      const TrilinosWrappers::SparseMatrix &B_,
      const TrilinosWrappers::SparseMatrix &Bt_,
      const TrilinosWrappers::MPI::BlockVector &owned_solution)
    {
      F = &F_;
      B = &B_;
      Bt = &Bt_;
      
      // Creating D inverse
      // This is a vector since the matrix D is diagonal
      Di.reinit(owned_solution.block(0));
      for (unsigned int i : Di.locally_owned_elements())
      {
        double temp = F->diag_element(i);
        Di[i] = 1.0 / temp;
      }

      // Creating S
      B->mmult(S, *Bt, Di);

      preconditioner_F.initialize(*F);
      preconditioner_S.initialize(S);

      vec0.reinit(owned_solution.block(0));
      vec1.reinit(owned_solution.block(1));
    }

    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      
      const unsigned int maxit = 10000;
      const double tol = 1.e-2;
      // It didn't work because i applied the preconditioner
      // starting from the wrong matrix

      // Solve F on block 0    
      SolverControl solver_control_F(maxit, tol * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
      solver_F.solve(*F, vec0, src.block(0), preconditioner_F);
      B->vmult(vec1, vec0);
      vec1.sadd(-1.0, src.block(1));

      // Solve the system in S on block 1
      SolverControl solver_control_S(maxit, tol * vec1.l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
      solver_S.solve(S, dst.block(1), vec1, preconditioner_S);
      dst.block(1) *= - 1.0 / alpha;

      Bt->vmult(dst.block(0), dst.block(1));
      dst.block(0).scale(Di);
      dst.block(0).sadd(- 1.0, vec0);
      
    }

  protected:
    const double alpha = 0.5;

    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B;
    const TrilinosWrappers::SparseMatrix *Bt;
    TrilinosWrappers::SparseMatrix S;

    mutable TrilinosWrappers::MPI::Vector Di;
    mutable TrilinosWrappers::MPI::Vector vec0;
    mutable TrilinosWrappers::MPI::Vector vec1;

    TrilinosWrappers::PreconditionILU preconditioner_F;
    TrilinosWrappers::PreconditionILU preconditioner_S;

  };

  
  class PreconditionAYosida
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void initialize(const TrilinosWrappers::SparseMatrix &F_,
      const TrilinosWrappers::SparseMatrix &B_,
      const TrilinosWrappers::SparseMatrix &Bt_,
      const TrilinosWrappers::MPI::BlockVector &DTMli_,
      const TrilinosWrappers::MPI::BlockVector &owned_solution)
    {
      F = &F_;
      B = &B_;
      Bt = &Bt_;
      DTMli = &DTMli_;

      // Creating S
      B->mmult(S, *Bt, DTMli->block(0));

      preconditioner_F.initialize(*F);
      preconditioner_S.initialize(S);

      vec0.reinit(owned_solution.block(0));
      vec1.reinit(owned_solution.block(1));
    }

    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      
      const unsigned int maxit = 10000;
      const double tol = 1.e-2;
      // It didn't work because i applied the preconditioner
      // starting from the wrong matrix

      // Solve F on block 0    
      SolverControl solver_control_F(maxit, tol * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
      solver_F.solve(*F, vec0, src.block(0), preconditioner_F);
      B->vmult(vec1, vec0);
      vec1.sadd(1.0, -1.0, src.block(1));

      // Solve the system in S on block 1
      SolverControl solver_control_S(maxit, tol * vec1.l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
      solver_S.solve(S, dst.block(1), vec1, preconditioner_S);
      //dst.block(1) *= - 1.0;

      Bt->vmult(dst.block(0), dst.block(1));
      solver_F.solve(*F, dst.block(0), dst.block(0), preconditioner_F);
      dst.block(0).sadd(- 1.0, vec0);
      
    }

  protected:

    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B;
    const TrilinosWrappers::SparseMatrix *Bt;
    TrilinosWrappers::SparseMatrix S;

    const TrilinosWrappers::MPI::BlockVector *DTMli;
    mutable TrilinosWrappers::MPI::Vector vec0;
    mutable TrilinosWrappers::MPI::Vector vec1;

    TrilinosWrappers::PreconditionILU preconditioner_F;
    TrilinosWrappers::PreconditionILU preconditioner_S;

  };



  NavierStokes(const std::string &mesh_file_name_,
    const unsigned int &degree_velocity_,
    const unsigned int &degree_pressure_,
    const double &deltat_,
    const double &T_,
    const unsigned int &step_)
    : mesh_file_name(mesh_file_name_)
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
    , mesh(MPI_COMM_WORLD)
    , deltat(deltat_)
    , T(T_)
    , step(step_)
  {}

  void setup();

  void assemble_matrices();

  void assemble_rhs(const double &time);

  void assemble_static_matrices();

  void solve_time_step();

  void assemble(const double &time);

  void output(const unsigned int &time_step) const;

  void solve();


protected:
  const std::string mesh_file_name;

  const unsigned int mpi_size;
  const unsigned int mpi_rank;
  ConditionalOStream pcout;

  const unsigned int degree_velocity;
  const unsigned int degree_pressure;


  parallel::fullydistributed::Triangulation<dim> mesh;
  std::unique_ptr<FiniteElement<dim>> fe;
  std::unique_ptr<Quadrature<dim>> quadrature;
  std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

  DoFHandler<dim> dof_handler;

  IndexSet locally_owned_dofs;
  std::vector<IndexSet> block_owned_dofs;

  IndexSet locally_relevant_dofs;
  std::vector<IndexSet> block_relevant_dofs;

  TrilinosWrappers::BlockSparseMatrix system_matrix;
  TrilinosWrappers::MPI::BlockVector deltat_lumped_mass_inv;
  TrilinosWrappers::MPI::BlockVector system_rhs;
  TrilinosWrappers::MPI::BlockVector solution_owned;
  TrilinosWrappers::MPI::BlockVector solution;

  const double nu = 1.e-3;
  const double p_out = 10.0;

  const double deltat;

  double time;
  const double T;
  const unsigned int step;

  ForcingTerm forcing_term;
  InletVelocity inlet_velocity;
  InitialConditions initial_conditions;
  FunctionH function_h;

  double drag, lift;
};

#endif