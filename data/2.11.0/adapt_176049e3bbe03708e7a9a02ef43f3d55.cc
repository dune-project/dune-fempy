#ifndef Guard_adapt_176049e3bbe03708e7a9a02ef43f3d55
#define Guard_adapt_176049e3bbe03708e7a9a02ef43f3d55

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fempy/py/grid/adaptation.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( adapt_176049e3bbe03708e7a9a02ef43f3d55, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::FemPy::GridAdaptation< Dune::ALUGrid< 2, 2, Dune::cube > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "GridAdaptation", Dune::Python::GenerateTypeName("Dune::FemPy::GridAdaptation< Dune::ALUGrid< 2, 2, Dune::cube > >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/fempy/py/grid/adaptation.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerGridAdaptation( cls0, cls );
  }
}
#endif