#ifndef Guard_view_2c4dc880cebe2f07a36e47627f499c79
#define Guard_view_2c4dc880cebe2f07a36e47627f499c79

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( view_2c4dc880cebe2f07a36e47627f499c79, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = typename Dune::ALUGrid< 2, 2, Dune::simplex >::LevelGridView;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "GridView", Dune::Python::GenerateTypeName("typename Dune::ALUGrid< 2, 2, Dune::simplex >::LevelGridView"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/python/grid/gridview.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::Python::registerGridView( cls0, cls );
  }
}
#endif