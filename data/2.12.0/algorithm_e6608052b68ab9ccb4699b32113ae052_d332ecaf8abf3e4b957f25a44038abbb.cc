#ifndef Guard_algorithm_e6608052b68ab9ccb4699b32113ae052_d332ecaf8abf3e4b957f25a44038abbb
#define Guard_algorithm_e6608052b68ab9ccb4699b32113ae052_d332ecaf8abf3e4b957f25a44038abbb

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/geometrygridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>


#include <dune/geometry/quadraturerules.hh>
template< class Surface >
double calcRadius( const Surface &surface ) {
  double R = 0, vol = 0.;
  for( const auto &entity : elements( surface ) ) {
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(entity.type(), 4);
    for ( const auto &p : rule ) {
      const auto geo = entity.geometry();
      const double weight = geo.volume() * p.weight();
      R   += geo.global(p.position()).two_norm() * weight;
      vol += weight;
    }
  }
  return R/vol;
}


PYBIND11_MODULE( algorithm_e6608052b68ab9ccb4699b32113ae052_d332ecaf8abf3e4b957f25a44038abbb, module )
{
  module.def( "run", [] ( Dune::Fem::GeometryGridPart< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > > & arg0 ) {
      return calcRadius( arg0 );
    } );
}
#endif
