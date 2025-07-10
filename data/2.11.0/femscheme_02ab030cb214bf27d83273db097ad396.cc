#ifndef Guard_femscheme_02ab030cb214bf27d83273db097ad396
#define Guard_femscheme_02ab030cb214bf27d83273db097ad396

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/common/quadmath.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>
#include <dune/fem/schemes/femscheme.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
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

PYBIND11_MODULE( femscheme_02ab030cb214bf27d83273db097ad396, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::FemScheme< DirichletWrapperOperator< DifferentiableVEMEllipticOperator< Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::SparseRowMatrix<double,int>>,Dune::Fem::VirtualizedVemIntegrands< typename Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double >::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >>,Dune::VemDirichletConstraints< Dune::Fem::VirtualizedVemIntegrands< typename Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double >::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >,Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double >,true,3 >  >, Dune::Fem::KrylovInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double > > > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Scheme",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::FemScheme< DirichletWrapperOperator< DifferentiableVEMEllipticOperator< Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double > >,Dune::Fem::SparseRowMatrix<double,int>>,Dune::Fem::VirtualizedVemIntegrands< typename Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double >::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >>,Dune::VemDirichletConstraints< Dune::Fem::VirtualizedVemIntegrands< typename Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double >::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >,Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double >,true,3 >  >, Dune::Fem::KrylovInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double > > > >"), Dune::Python::IncludeFiles{"dune/common/quadmath.hh","dune/fem/function/adaptivefunction.hh","dune/fem/gridpart/common/gridpart.hh","dune/fem/operator/linear/spoperator.hh","dune/fem/schemes/dirichletwrapper.hh","dune/fem/schemes/femscheme.hh","dune/fem/solver/krylovinverseoperators.hh","dune/fempy/parameter.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/scheme.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/misc/grid.hh","dune/vem/operator/diffusionmodel.hh","dune/vem/operator/vemdirichletconstraints.hh","dune/vem/operator/vemelliptic.hh","dune/vem/space/hk.hh"}).first;
    Dune::FemPy::registerScheme( cls0, cls );
  }
}
#endif