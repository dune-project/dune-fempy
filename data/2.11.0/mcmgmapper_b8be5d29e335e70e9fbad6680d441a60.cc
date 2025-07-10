#ifndef Guard_mcmgmapper_b8be5d29e335e70e9fbad6680d441a60
#define Guard_mcmgmapper_b8be5d29e335e70e9fbad6680d441a60

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/grid/mapper.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( mcmgmapper_b8be5d29e335e70e9fbad6680d441a60, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::MultipleCodimMultipleGeomTypeMapper< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "MultipleCodimMultipleGeomTypeMapper", Dune::Python::GenerateTypeName("Dune::MultipleCodimMultipleGeomTypeMapper< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/python/grid/hierarchical.hh","dune/python/grid/mapper.hh"}).first;
    Dune::Python::registerMultipleCodimMultipleGeomTypeMapper( cls0, cls );
  }
}
#endif