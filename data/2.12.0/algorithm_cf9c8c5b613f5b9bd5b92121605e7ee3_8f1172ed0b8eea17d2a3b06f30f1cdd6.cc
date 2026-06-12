#ifndef Guard_algorithm_cf9c8c5b613f5b9bd5b92121605e7ee3_8f1172ed0b8eea17d2a3b06f30f1cdd6
#define Guard_algorithm_cf9c8c5b613f5b9bd5b92121605e7ee3_8f1172ed0b8eea17d2a3b06f30f1cdd6

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>


#include <dune/fem/misc/boundaryidprovider.hh>

template <class DF>
void inspectBoundary( DF& df )
{
  Dune::Fem::projectBoundaryIds( df );
}


PYBIND11_MODULE( algorithm_cf9c8c5b613f5b9bd5b92121605e7ee3_8f1172ed0b8eea17d2a3b06f30f1cdd6, module )
{
  module.def( "run", [] ( Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> > & arg0 ) {
      return inspectBoundary( arg0 );
    } );
}
#endif
