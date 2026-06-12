#ifndef Guard_algorithm_ffe434ec87b4fb8d8045c4eb008446ba_30e6c9582f74f843ffa8db632cc01675
#define Guard_algorithm_ffe434ec87b4fb8d8045c4eb008446ba_30e6c9582f74f843ffa8db632cc01675

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>


    #include <dune/fem/misc/boundaryidprovider.hh>

    template <class GridView, class Intersection>
    int boundaryId( const GridView& gv, const Intersection& i )
    {
      return Dune::Fem::boundaryId( gv, i );
    }


PYBIND11_MODULE( algorithm_ffe434ec87b4fb8d8045c4eb008446ba_30e6c9582f74f843ffa8db632cc01675, module )
{
  module.def( "run", [] ( Dune::ALUGrid< 2, 2, Dune::cube >::LeafGridView & arg0, Dune::ALUGrid< 2, 2, Dune::cube >::LeafGridView::Intersection & arg1 ) {
      return boundaryId( arg0, arg1 );
    } );
}
#endif
