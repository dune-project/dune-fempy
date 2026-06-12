#ifndef Guard_hierarchicalgrid_ecbddb274cb86041151a853587dceb93
#define Guard_hierarchicalgrid_ecbddb274cb86041151a853587dceb93

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( hierarchicalgrid_ecbddb274cb86041151a853587dceb93, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::ALUGrid< 2, 2, Dune::cube >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "HierarchicalGrid",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::ALUGrid< 2, 2, Dune::cube >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::Python::registerHierarchicalGrid( cls0, cls );
  }
}
#endif
