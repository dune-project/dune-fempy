#ifndef Guard_femscheme_4901d05833ff7fe6f88ecd37363f35ba
#define Guard_femscheme_4901d05833ff7fe6f88ecd37363f35ba

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>
#include <dune/fem/schemes/molgalerkin.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/scheme.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femscheme_4901d05833ff7fe6f88ecd37363f35ba, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::GalerkinScheme< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::Fem::CodegenStorage>::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >, Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::Fem::CodegenStorage> >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::Fem::CodegenStorage> >,Dune::Fem::SparseRowMatrix<double,int>>, Dune::Fem::KrylovInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::Fem::CodegenStorage> > >, true >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Scheme",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::GalerkinScheme< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::Fem::CodegenStorage>::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >, Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::Fem::CodegenStorage> >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::Fem::CodegenStorage> >,Dune::Fem::SparseRowMatrix<double,int>>, Dune::Fem::KrylovInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::Fem::CodegenStorage> > >, true >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/spoperator.hh","dune/fem/schemes/dirichletwrapper.hh","dune/fem/schemes/molgalerkin.hh","dune/fem/solver/krylovinverseoperators.hh","dune/fem/space/lagrange.hh","dune/fempy/parameter.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/scheme.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/misc/grid.hh"}).first;
    Dune::FemPy::registerScheme( cls0, cls );
  }
}
#endif
