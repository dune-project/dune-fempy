#ifndef Guard_femoperator_4c44ccca08ea9eaf3bd959ef369eb801
#define Guard_femoperator_4c44ccca08ea9eaf3bd959ef369eb801

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem-dg/operator/common/orderreduction.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/combinedspace/tuplespace.hh>
#include <dune/fem/space/hpdg/legendre.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/operator.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( femoperator_4c44ccca08ea9eaf3bd959ef369eb801, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::OrderReduction< Dune::Fem::TupleDiscreteFunction< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Operator",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::OrderReduction< Dune::Fem::TupleDiscreteFunction< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > > >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/fem-dg/operator/common/orderreduction.hh","dune/fem/function/adaptivefunction.hh","dune/fem/function/tuplediscretefunction.hh","dune/fem/gridpart/adaptiveleafgridpart.hh","dune/fem/space/combinedspace/tuplespace.hh","dune/fem/space/hpdg/legendre.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/gridview.hh","dune/fempy/py/operator.hh","dune/fempy/py/space.hh","dune/python/grid/gridview.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerOperator( cls0, cls );
    cls.def( pybind11::init( [] ( const Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage>, Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > &dSpace ) {return new Dune::Fem::OrderReduction< Dune::Fem::TupleDiscreteFunction< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > > >(dSpace);
      } ), "dSpace"_a, pybind11::keep_alive< 1, 2 >() );
  }
}
#endif