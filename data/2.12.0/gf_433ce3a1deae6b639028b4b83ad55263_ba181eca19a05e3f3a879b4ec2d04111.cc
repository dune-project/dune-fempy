#ifndef Guard_gf_433ce3a1deae6b639028b4b83ad55263_ba181eca19a05e3f3a879b4ec2d04111
#define Guard_gf_433ce3a1deae6b639028b4b83ad55263_ba181eca19a05e3f3a879b4ec2d04111

#include <config.h>

#define USING_DUNE_PYTHON 1


#include <cmath>
#include <dune/python/pybind11/numpy.h>
template <class GridView>
auto bnd(pybind11::array_t<double> &a, double h) {
  return [ptr=static_cast<double*>(a.request(true).ptr),h](const auto& en,const auto& xLocal) -> auto {
    auto x = en.geometry().global(xLocal);
    return sin(h*M_PI*(ptr[0]*x[0]+ptr[1]*x[1]));
  };
}

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/pybind11/numpy.h>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_433ce3a1deae6b639028b4b83ad55263_ba181eca19a05e3f3a879b4ec2d04111, module )
{
  module.def( "gf", [module] ( Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView &gv, pybind11::array_t<double> arg0, double arg1 ) {
      auto callback=bnd<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>( arg0,arg1);
      return Dune::Python::registerGridFunction<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView,decltype(callback)>(module,pybind11::cast(gv),"tmp",callback);
    },    pybind11::keep_alive<0,1>(),pybind11::keep_alive<0,2>(),pybind11::keep_alive<0,3>());
}
#endif
