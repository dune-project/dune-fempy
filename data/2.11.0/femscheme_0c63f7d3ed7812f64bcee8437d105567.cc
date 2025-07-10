#ifndef Guard_femscheme_0c63f7d3ed7812f64bcee8437d105567
#define Guard_femscheme_0c63f7d3ed7812f64bcee8437d105567

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/common/quadmath.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>
#include <dune/fem/schemes/femscheme.hh>
#include <dune/fem/solver/ldlsolver.hh>
#include <dune/fem/solver/spqrsolver.hh>
#include <dune/fem/solver/umfpacksolver.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/scheme.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>
#include <dune/vem/operator/diffusionmodel.hh>
#include <dune/vem/operator/vemdirichletconstraints.hh>
#include <dune/vem/operator/vemelliptic.hh>
#include <dune/vem/space/hk.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( femscheme_0c63f7d3ed7812f64bcee8437d105567, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::FemScheme< DirichletWrapperOperator< DifferentiableVEMEllipticOperator< Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::SparseRowMatrix<double,int>>,Dune::Fem::VirtualizedVemIntegrands< typename Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double >::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > >>,Dune::VemDirichletConstraints< Dune::Fem::VirtualizedVemIntegrands< typename Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double >::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > >,Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double >,true,2 >  >, Dune::Fem::UMFPACKInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >, typename Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::SparseRowMatrix<double,int>>::MatrixType > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Scheme",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::FemScheme< DirichletWrapperOperator< DifferentiableVEMEllipticOperator< Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::SparseRowMatrix<double,int>>,Dune::Fem::VirtualizedVemIntegrands< typename Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double >::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > >>,Dune::VemDirichletConstraints< Dune::Fem::VirtualizedVemIntegrands< typename Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double >::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > >,Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double >,true,2 >  >, Dune::Fem::UMFPACKInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >, typename Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::SparseRowMatrix<double,int>>::MatrixType > >"), Dune::Python::IncludeFiles{"dune/common/quadmath.hh","dune/fem/function/adaptivefunction.hh","dune/fem/gridpart/common/gridpart.hh","dune/fem/operator/linear/spoperator.hh","dune/fem/schemes/dirichletwrapper.hh","dune/fem/schemes/femscheme.hh","dune/fem/solver/ldlsolver.hh","dune/fem/solver/spqrsolver.hh","dune/fem/solver/umfpacksolver.hh","dune/fempy/parameter.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/scheme.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/misc/grid.hh","dune/vem/operator/diffusionmodel.hh","dune/vem/operator/vemdirichletconstraints.hh","dune/vem/operator/vemelliptic.hh","dune/vem/space/hk.hh"}).first;
    Dune::FemPy::registerScheme( cls0, cls );
  }
}
#endif