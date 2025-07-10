#ifndef Guard_indexset_08aa6e5d32abc9657a050bbd75c2c7c6
#define Guard_indexset_08aa6e5d32abc9657a050bbd75c2c7c6

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/python/grid/hierarchical.hh>
#include <dune/python/grid/indexset.hh>
#include <dune/vem/misc/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( indexset_08aa6e5d32abc9657a050bbd75c2c7c6, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Vem::Grid<2,2>::LeafGridView::IndexSet;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "GridViewIndexSet", Dune::Python::GenerateTypeName("Dune::Vem::Grid<2,2>::LeafGridView::IndexSet"), Dune::Python::IncludeFiles{"dune/python/grid/hierarchical.hh","dune/python/grid/indexset.hh","dune/vem/misc/grid.hh"}).first;
    Dune::Python::registerGridViewIndexSet( cls0, cls );
  }
}
#endif