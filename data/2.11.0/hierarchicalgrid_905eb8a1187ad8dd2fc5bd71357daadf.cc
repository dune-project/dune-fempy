#ifndef Guard_hierarchicalgrid_905eb8a1187ad8dd2fc5bd71357daadf
#define Guard_hierarchicalgrid_905eb8a1187ad8dd2fc5bd71357daadf

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/polygongrid/dgf.hh>
#include <dune/polygongrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( hierarchicalgrid_905eb8a1187ad8dd2fc5bd71357daadf, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::PolygonGrid< double >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "HierarchicalGrid",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::PolygonGrid< double >"), Dune::Python::IncludeFiles{"dune/polygongrid/dgf.hh","dune/polygongrid/grid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::Python::registerHierarchicalGrid( cls0, cls );
    cls.def( "dualGrid", []( DuneType &self ) { return self.dualGrid(); } );
    cls.def( "cachingStorage", []( DuneType &self ) { return false; } );
  }
}
#endif
