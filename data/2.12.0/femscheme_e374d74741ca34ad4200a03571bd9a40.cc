#ifndef Guard_femscheme_e374d74741ca34ad4200a03571bd9a40
#define Guard_femscheme_e374d74741ca34ad4200a03571bd9a40

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>
#include <dune/fem/schemes/molgalerkin.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
#include <dune/fem/space/combinedspace/tuplespace.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/scheme.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femscheme_e374d74741ca34ad4200a03571bd9a40, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::GalerkinScheme< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> >::GridPartType, std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 2 > >, std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 2 > > >, Dune::Fem::ISTLLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> > >>, Dune::Fem::ISTLInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> > > >, false >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Scheme",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::GalerkinScheme< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> >::GridPartType, std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 2 > >, std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 2 > > >, Dune::Fem::ISTLLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> > >>, Dune::Fem::ISTLInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> > > >, false >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/istloperator.hh","dune/fem/schemes/dirichletwrapper.hh","dune/fem/schemes/molgalerkin.hh","dune/fem/solver/istlinverseoperators.hh","dune/fem/space/combinedspace/tuplespace.hh","dune/fem/space/discontinuousgalerkin.hh","dune/fempy/parameter.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/scheme.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerScheme( cls0, cls );
  }
}
#endif
