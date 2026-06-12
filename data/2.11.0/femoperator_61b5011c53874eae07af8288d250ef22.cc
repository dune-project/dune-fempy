#ifndef Guard_femoperator_61b5011c53874eae07af8288d250ef22
#define Guard_femoperator_61b5011c53874eae07af8288d250ef22

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/schemes/galerkin.hh>
#include <dune/fem/space/combinedspace/tuplespace.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/hpdg/legendre.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/operator.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( femoperator_61b5011c53874eae07af8288d250ef22, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::GalerkinOperator< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 0 , Dune::Fem::SimpleStorage>::GridPartType, std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 2 > >, std::tuple< Dune::FieldVector< double, 1 > > >, Dune::Fem::TupleDiscreteFunction< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > >, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 0 , Dune::Fem::SimpleStorage> > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Operator",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::GalerkinOperator< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 0 , Dune::Fem::SimpleStorage>::GridPartType, std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 2 > >, std::tuple< Dune::FieldVector< double, 1 > > >, Dune::Fem::TupleDiscreteFunction< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > >, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 0 , Dune::Fem::SimpleStorage> > >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/fem/function/adaptivefunction.hh","dune/fem/function/tuplediscretefunction.hh","dune/fem/gridpart/adaptiveleafgridpart.hh","dune/fem/operator/linear/spoperator.hh","dune/fem/schemes/galerkin.hh","dune/fem/space/combinedspace/tuplespace.hh","dune/fem/space/finitevolume.hh","dune/fem/space/hpdg/legendre.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/grid/gridpart.hh","dune/fempy/py/gridview.hh","dune/fempy/py/operator.hh","dune/fempy/py/space.hh","dune/python/grid/gridview.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerOperator( cls0, cls );
    cls.def( "setCommunicate", []( DuneType &self, const bool communicate ) { self.setCommunicate( communicate ); } );
    cls.def( pybind11::init( [] ( const Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage>, Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >& dSpace, const Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 0 , Dune::Fem::SimpleStorage> &rSpace, Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 0 , Dune::Fem::SimpleStorage>::GridPartType, std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 2 > >, std::tuple< Dune::FieldVector< double, 1 > > > &integrands ) {return new DuneType( dSpace.gridPart(), integrands );
      } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
  }
}
#endif
