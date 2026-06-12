#ifndef Guard_hierarchicalgrid_157fac4e8cf46b26123765e475f12631
#define Guard_hierarchicalgrid_157fac4e8cf46b26123765e475f12631

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/grid/io/file/dgfparser/dgfoned.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( hierarchicalgrid_157fac4e8cf46b26123765e475f12631, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::OneDGrid;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "HierarchicalGrid",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::OneDGrid"), Dune::Python::IncludeFiles{"dune/grid/io/file/dgfparser/dgfoned.hh","dune/grid/onedgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::Python::registerHierarchicalGrid( cls0, cls );
  }
}
#endif
