#ifndef Guard_hierarchicalgrid_a59f54d38d193cb8760732eb4a8edbfc
#define Guard_hierarchicalgrid_a59f54d38d193cb8760732eb4a8edbfc

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( hierarchicalgrid_a59f54d38d193cb8760732eb4a8edbfc, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "HierarchicalGrid",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >"), Dune::Python::IncludeFiles{"dune/grid/spgrid.hh","dune/grid/spgrid/dgfparser.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::Python::registerHierarchicalGrid( cls0, cls );
  }
}
#endif