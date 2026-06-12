#ifndef Guard_algorithm_ab2ef914b9e55b8b64a10df1cc8d44e7_30e6c9582f74f843ffa8db632cc01675
#define Guard_algorithm_ab2ef914b9e55b8b64a10df1cc8d44e7_30e6c9582f74f843ffa8db632cc01675

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>


    #include <dune/fem/misc/boundaryidprovider.hh>

    template <class GridView, class Intersection>
    int boundaryId( const GridView& gv, const Intersection& i )
    {
      return Dune::Fem::boundaryId( gv, i );
    }


PYBIND11_MODULE( algorithm_ab2ef914b9e55b8b64a10df1cc8d44e7_30e6c9582f74f843ffa8db632cc01675, module )
{
  module.def( "run", [] ( Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView & arg0, Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::Intersection & arg1 ) {
      return boundaryId( arg0, arg1 );
    } );
}
#endif
