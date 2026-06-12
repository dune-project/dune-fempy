#ifndef Guard_spaceadapt_28fe49a0f3c3ef85e10b2583f2c8f1d1
#define Guard_spaceadapt_28fe49a0f3c3ef85e10b2583f2c8f1d1

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/hpdg/legendre.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/fempy/py/spaceadapt.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( spaceadapt_28fe49a0f3c3ef85e10b2583f2c8f1d1, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::FemPy::SpaceAdaptation< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "SpaceAdaptation", Dune::Python::GenerateTypeName("Dune::FemPy::SpaceAdaptation< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/fem/function/adaptivefunction.hh","dune/fem/gridpart/adaptiveleafgridpart.hh","dune/fem/space/hpdg/legendre.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/gridview.hh","dune/fempy/py/space.hh","dune/fempy/py/spaceadapt.hh","dune/python/grid/gridview.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerSpaceAdaptation( cls0, cls );
  }
}
#endif
