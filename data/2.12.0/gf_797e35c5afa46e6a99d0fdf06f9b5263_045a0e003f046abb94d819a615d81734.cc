#ifndef Guard_gf_797e35c5afa46e6a99d0fdf06f9b5263_045a0e003f046abb94d819a615d81734
#define Guard_gf_797e35c5afa46e6a99d0fdf06f9b5263_045a0e003f046abb94d819a615d81734

#include <config.h>

#define USING_DUNE_PYTHON 1


#include <cmath>
#include <dune/common/fvector.hh>
template <class GridView>
auto aTimesExact(double a) {
  return [a](const auto& en,const auto& xLocal) -> auto {
    auto x = en.geometry().global(xLocal);
    return a*(1./2.*(std::pow(x[0],2)+std::pow(x[1],2)) - 1./3.*(std::pow(x[0],3) - std::pow(x[1],3)) + 1.);
  };
}

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_797e35c5afa46e6a99d0fdf06f9b5263_045a0e003f046abb94d819a615d81734, module )
{
  module.def( "gf", [module] ( Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView &gv, double arg0 ) {
      auto callback=aTimesExact<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>( arg0);
      return Dune::Python::registerGridFunction<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView,decltype(callback)>(module,pybind11::cast(gv),"tmp",callback);
    },    pybind11::keep_alive<0,1>(),pybind11::keep_alive<0,2>());
}
#endif
