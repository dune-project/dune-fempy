#ifndef Guard_femscheme_9b7005253a7d0b23b665ce4b07198dde
#define Guard_femscheme_9b7005253a7d0b23b665ce4b07198dde

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>
#include <dune/fem/schemes/molgalerkin.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/scheme.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( femscheme_9b7005253a7d0b23b665ce4b07198dde, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::GalerkinScheme< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >, Dune::Fem::ISTLLinearOperator< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >,Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >>, Dune::Fem::ISTLInverseOperator< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> > >, false >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Scheme",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::GalerkinScheme< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >, Dune::Fem::ISTLLinearOperator< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >,Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >>, Dune::Fem::ISTLInverseOperator< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> > >, false >"), Dune::Python::IncludeFiles{"dune/fem/function/blockvectorfunction.hh","dune/fem/operator/linear/istloperator.hh","dune/fem/schemes/dirichletwrapper.hh","dune/fem/schemes/molgalerkin.hh","dune/fem/solver/istlinverseoperators.hh","dune/fem/space/lagrange.hh","dune/fempy/parameter.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/scheme.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerScheme( cls0, cls );
  }
}
#endif
