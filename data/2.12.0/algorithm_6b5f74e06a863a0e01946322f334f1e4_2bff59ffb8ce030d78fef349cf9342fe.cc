#ifndef Guard_algorithm_6b5f74e06a863a0e01946322f334f1e4_2bff59ffb8ce030d78fef349cf9342fe
#define Guard_algorithm_6b5f74e06a863a0e01946322f334f1e4_2bff59ffb8ce030d78fef349cf9342fe

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/filter/simple.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>


#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/common/bindguard.hh>
template <class GF, class DF>
void filteredInterpolate( const GF &gf, DF& df )
{
  Dune::Fem::LocalInterpolation< typename DF::DiscreteFunctionSpaceType > interpolation( df.space() );
  Dune::Fem::ConstLocalFunction< GF > gfLocal( gf );
  Dune::Fem::SetSelectedLocalContribution< DF > dfLocal( df );
  for( const auto &entity : df.space() ) // iterate over filter coming from df
  {
    auto gfGuard = Dune::Fem::bindGuard( gfLocal, entity ); // note: entities on filtergp are not wrapped
    auto dfGuard = Dune::Fem::bindGuard( dfLocal, entity );
    auto iGuard  = Dune::Fem::bindGuard( interpolation, entity );
    interpolation(gfLocal,dfLocal);
  }
}


PYBIND11_MODULE( algorithm_6b5f74e06a863a0e01946322f334f1e4_2bff59ffb8ce030d78fef349cf9342fe, module )
{
  module.def( "run", [] ( Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > & arg0, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::FilteredGridPart< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::Fem::SimpleFilter< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > >, true > >, Dune::Fem::CodegenStorage> > & arg1 ) {
      return filteredInterpolate( arg0, arg1 );
    } );
}
#endif
