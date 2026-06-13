#ifndef Guard_femoperator_766cdba53035482d2041b76ad5f503da
#define Guard_femoperator_766cdba53035482d2041b76ad5f503da

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem-dg/solver/rungekuttasolver.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/operator.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femoperator_766cdba53035482d2041b76ad5f503da, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::SimpleRungeKuttaSolver< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >>;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Operator",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::SimpleRungeKuttaSolver< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >>"), Dune::Python::IncludeFiles{"dune/fem-dg/solver/rungekuttasolver.hh","dune/fem/function/adaptivefunction.hh","dune/fem/space/finitevolume.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/operator.hh","dune/fempy/py/space.hh","dune/grid/spgrid.hh","dune/grid/spgrid/dgfparser.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerOperator( cls0, cls );
    cls.def( pybind11::init( [] ( Dune::Fem::SpaceOperatorInterface< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >> &op, Dune::Fem::SpaceOperatorInterface< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >> &explOp, Dune::Fem::SpaceOperatorInterface< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >> &implOp, const int imexId, const pybind11::dict &parameters ) {return new Dune::Fem::SimpleRungeKuttaSolver< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >>(op, explOp, implOp, imexId, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ));
      } ), "op"_a, "explOp"_a, "implOp"_a, "imexId"_a, "parameters"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >(), pybind11::keep_alive< 1, 6>() );
    cls.def( "solve", []( DuneType &self, typename DuneType::DestinationType &u) { self.solve(u); } );
    cls.def( "setTimeStepSize", &DuneType::setTimeStepSize );
    cls.def( "deltaT", &DuneType::deltaT );
    cls.def( "stages", &DuneType::stages );
  }
}
#endif