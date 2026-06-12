#ifndef Guard_algorithm_1f9546f3ea5f15618b40e32978119e8a_8f1172ed0b8eea17d2a3b06f30f1cdd6
#define Guard_algorithm_1f9546f3ea5f15618b40e32978119e8a_8f1172ed0b8eea17d2a3b06f30f1cdd6

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/filter/simple.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>


#include <dune/fem/misc/boundaryidprovider.hh>

template <class DF>
void inspectBoundary( DF& df )
{
  Dune::Fem::projectBoundaryIds( df );
}


PYBIND11_MODULE( algorithm_1f9546f3ea5f15618b40e32978119e8a_8f1172ed0b8eea17d2a3b06f30f1cdd6, module )
{
  module.def( "run", [] ( Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::FilteredGridPart< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::Fem::SimpleFilter< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > >, true > >, 0 , Dune::Fem::SimpleStorage> > & arg0 ) {
      return inspectBoundary( arg0 );
    } );
}
#endif
