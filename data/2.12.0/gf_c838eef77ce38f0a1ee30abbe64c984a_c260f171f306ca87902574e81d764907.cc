#ifndef Guard_gf_c838eef77ce38f0a1ee30abbe64c984a_c260f171f306ca87902574e81d764907
#define Guard_gf_c838eef77ce38f0a1ee30abbe64c984a_c260f171f306ca87902574e81d764907

#include <config.h>

#define USING_DUNE_PYTHON 1


#include <cmath>
#include <dune/common/fvector.hh>
template <class GridView, class GF>
auto aTimesExact(const GF &gf,Dune::FieldVector<double,1> &a) {
  return [lgf=localFunction(gf),&a] (const auto& en,const auto& xLocal) mutable -> auto {

    lgf.bind(en); // lambda must be mutable so that the non const function can be called
    return a[0]*lgf(xLocal);
  };
}

#include <dune/common/fvector.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <python/dune/generated/localfunction_23f5064f900e951408b9c650ab73109d_19659fe2f2f74bcfb5d532f1a0d38be9.cc>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_c838eef77ce38f0a1ee30abbe64c984a_c260f171f306ca87902574e81d764907, module )
{
  module.def( "gf", [module] ( Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView &gv, UFLLocalFunctions_23f5064f900e951408b9c650ab73109d::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > & arg0, Dune::FieldVector< double, 1 > & arg1 ) {
      auto callback=aTimesExact<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>( arg0,arg1);
      return Dune::Python::registerGridFunction<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView,decltype(callback)>(module,pybind11::cast(gv),"tmp",callback);
    },    pybind11::keep_alive<0,1>(),pybind11::keep_alive<0,2>(),pybind11::keep_alive<0,3>());
}
#endif
