#ifndef Guard_algorithm_1d94df8d475f997b0a12511ed7a73975_4d8db901b63fef6236c55e0085cd3e73
#define Guard_algorithm_1d94df8d475f997b0a12511ed7a73975_4d8db901b63fef6236c55e0085cd3e73

#include <config.h>

#define USING_DUNE_PYTHON 1


#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/common/bindguard.hh>
template <class GF, class DF>
void filteredInterpolate( const GF &gf, DF& df )
{
  Dune::Fem::LocalInterpolation< typename DF::DiscreteFunctionSpaceType > interpolation( df.space() );
  Dune::Fem::ConstLocalFunction< GF > gfLocal( gf );
  Dune::Fem::SetLocalContribution< DF > dfLocal( df );
  for( const auto &entity : gf.space() ) // iterate over filter coming from gf
  {
    auto gfGuard = Dune::Fem::bindGuard( gfLocal, entity );
    auto dfGuard = Dune::Fem::bindGuard( dfLocal, entity ); // note: entities on filtergp are not wrapped
    auto iGuard  = Dune::Fem::bindGuard( interpolation, entity );
    interpolation(gfLocal,dfLocal);
  }
}

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/filter/simple.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( algorithm_1d94df8d475f997b0a12511ed7a73975_4d8db901b63fef6236c55e0085cd3e73, module )
{
  module.def( "run", [] ( Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::FilteredGridPart< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::Fem::SimpleFilter< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > >, true > >, Dune::Fem::CodegenStorage> > & arg0, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > & arg1 ) {
      return filteredInterpolate( arg0, arg1 );
    } );
}
#endif
