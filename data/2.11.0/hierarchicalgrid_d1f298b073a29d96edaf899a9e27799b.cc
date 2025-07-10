#ifndef Guard_hierarchicalgrid_d1f298b073a29d96edaf899a9e27799b
#define Guard_hierarchicalgrid_d1f298b073a29d96edaf899a9e27799b

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( hierarchicalgrid_d1f298b073a29d96edaf899a9e27799b, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Vem::Grid<2,2>;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "HierarchicalGrid",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Vem::Grid<2,2>"), Dune::Python::IncludeFiles{"dune/python/grid/hierarchical.hh","dune/vem/misc/grid.hh"}).first;
    Dune::Python::registerHierarchicalGrid( cls0, cls );
  }
}
#endif