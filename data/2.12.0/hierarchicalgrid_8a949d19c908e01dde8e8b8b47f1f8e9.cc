#ifndef Guard_hierarchicalgrid_8a949d19c908e01dde8e8b8b47f1f8e9
#define Guard_hierarchicalgrid_8a949d19c908e01dde8e8b8b47f1f8e9

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( hierarchicalgrid_8a949d19c908e01dde8e8b8b47f1f8e9, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Vem::CubeGrid<2,2>;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "HierarchicalGrid",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Vem::CubeGrid<2,2>"), Dune::Python::IncludeFiles{"dune/python/grid/hierarchical.hh","dune/vem/misc/grid.hh"}).first;
    Dune::Python::registerHierarchicalGrid( cls0, cls );
  }
}
#endif
