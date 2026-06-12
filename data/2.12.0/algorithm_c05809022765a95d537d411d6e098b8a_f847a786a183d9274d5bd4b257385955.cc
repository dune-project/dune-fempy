#ifndef Guard_algorithm_c05809022765a95d537d411d6e098b8a_f847a786a183d9274d5bd4b257385955
#define Guard_algorithm_c05809022765a95d537d411d6e098b8a_f847a786a183d9274d5bd4b257385955

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/common/fvector.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <python/dune/generated/localfunction_6fc6a1a15c1fbc7b790bd9cb62a416b2_af122c1df944c95cd395ec0f91d0f970.cc>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>

#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/common/bindguard.hh>

// Note: the coefficients of the functional will be ADDED to the storage of the functional
template <class Point, class Functional, class Error>
double pointFunctional(const Point &point, Functional &functional, Error &error)
{
  // first find the entity containing `point`
  typedef typename Functional::DiscreteFunctionSpaceType::GridPartType GridPartType;
  Dune::Fem::EntitySearch<GridPartType> search(functional.space().gridPart());
  const auto &entity = search(point);
  const auto localPoint = entity.geometry().local(point);

  // add the contributions from the basis functions to `functional`
  Dune::Fem::AddLocalContribution< Functional > wLocal( functional );
  {
    auto guard = Dune::Fem::bindGuard( wLocal, entity );
    typename Functional::RangeType one( 1 );
    wLocal.axpy(localPoint, one);
  }

  // in addition we also want to compute the value of `error` at `point`
  auto localError = constLocalFunction(error);
  localError.bind(entity);
  return localError.evaluate(localPoint)[0];
}


PYBIND11_MODULE( algorithm_c05809022765a95d537d411d6e098b8a_f847a786a183d9274d5bd4b257385955, module )
{
  module.def( "run", [] ( Dune::FieldVector< double, 2 > & arg0, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > >, Dune::Fem::CodegenStorage> > & arg1, UFLLocalFunctions_6fc6a1a15c1fbc7b790bd9cb62a416b2::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > > & arg2 ) {
      return pointFunctional( arg0, arg1, arg2 );
    } );
}
#endif
