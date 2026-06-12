#ifndef Guard_hierarchicalgrid_2fbf00a9e3ca0098fcbe4e884ec704ee
#define Guard_hierarchicalgrid_2fbf00a9e3ca0098fcbe4e884ec704ee

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( hierarchicalgrid_2fbf00a9e3ca0098fcbe4e884ec704ee, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "HierarchicalGrid",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >"), Dune::Python::IncludeFiles{"dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::Python::registerHierarchicalGrid( cls0, cls );
    cls.def( pybind11::init( [] ( const Dune::EquidistantOffsetCoordinates< double, 1 >& coordinates, std::array<bool, 1> periodic, int overlap ) {std::bitset<1> periodic_;
        for (int i=0;i<1;++i) periodic_.set(i,periodic[i]);
        return new DuneType(coordinates,periodic_,overlap);
      } ), "coordinates"_a, "periodic"_a, "overlap"_a );
  }
}
#endif
