#ifndef Guard_algorithm_18e0af3f7b117e06a5fbb3b9f03fca14_01e66ccfe515e42e55cf21194344523b
#define Guard_algorithm_18e0af3f7b117e06a5fbb3b9f03fca14_01e66ccfe515e42e55cf21194344523b

#include <config.h>

#define USING_DUNE_PYTHON 1


    #include <dune/fem/misc/boundaryidprovider.hh>

    template <class GridView, class Intersection>
    int boundaryId( const GridView& gv, const Intersection& i )
    {
      return Dune::Fem::boundaryId( gv, i );
    }

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( algorithm_18e0af3f7b117e06a5fbb3b9f03fca14_01e66ccfe515e42e55cf21194344523b, module )
{
  module.def( "run", [] ( Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView & arg0, Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView::Intersection & arg1 ) {
      return boundaryId( arg0, arg1 );
    } );
}
#endif
