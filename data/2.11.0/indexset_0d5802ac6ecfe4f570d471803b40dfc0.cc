#ifndef Guard_indexset_0d5802ac6ecfe4f570d471803b40dfc0
#define Guard_indexset_0d5802ac6ecfe4f570d471803b40dfc0

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/grid/indexset.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( indexset_0d5802ac6ecfe4f570d471803b40dfc0, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView::IndexSet;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "GridViewIndexSet", Dune::Python::GenerateTypeName("Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView::IndexSet"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/python/grid/hierarchical.hh","dune/python/grid/indexset.hh"}).first;
    Dune::Python::registerGridViewIndexSet( cls0, cls );
  }
}
#endif