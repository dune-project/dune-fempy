#ifndef Guard_femoperator_279225db56c27a4f65fd4a69af85c38c
#define Guard_femoperator_279225db56c27a4f65fd4a69af85c38c

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem-dg/operator/dg/dgpyoperator.hh>
#include <dune/fem-dg/operator/limiter/indicatorbase.hh>
#include <dune/fem-dg/python/operator.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/schemes/conservationlawmodel.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/operator.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>
#include <dune/python/grid/hierarchical.hh>
#include <python/dune/generated/modelFemDgAdv_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a.cc>
#include <python/dune/generated/modelFemDgDiff_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a.cc>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Additional_2b612bf20bda652dcb1a22d48081ea34
// -------------------------------------------

template< class FunctionSpace >
struct Additional_2b612bf20bda652dcb1a22d48081ea34
{
  typedef typename FunctionSpace::DomainType DomainType;
  typedef typename FunctionSpace::RangeType RangeType;
  typedef typename FunctionSpace::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpace::HessianRangeType HessianRangeType;
  static const int limitedDimRange = FunctionSpace :: dimRange;
  static const bool hasAdvection = true;
  static const bool hasDiffusion = false;
  static const bool hasStiffSource = false;
  static const bool hasNonStiffSource = false;
  static const bool hasFlux = true;
  static const bool threading = false;
  static const Dune::Fem::Solver::Enum solverId = Dune::Fem::Solver::Enum::fem;
  static const Dune::Fem::Formulation::Enum formId = Dune::Fem::Formulation::Enum::primal;
  static const Dune::Fem::AdvectionLimiter::Enum limiterId = Dune::Fem::AdvectionLimiter::Enum::limited;
  static const Dune::Fem::AdvectionLimiterFunction::Enum limiterFunctionId = Dune::Fem::AdvectionLimiterFunction::Enum::minmod;
  static const Dune::Fem::AdvectionFlux::Enum advFluxId = Dune::Fem::AdvectionFlux::Enum::upwind;
  static const Dune::Fem::DiffusionFlux::Enum diffFluxId = Dune::Fem::DiffusionFlux::Enum::none;
  static const bool defaultQuadrature = true;
};

PYBIND11_MODULE( femoperator_279225db56c27a4f65fd4a69af85c38c, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::DGOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >, ModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > >, ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > >, Additional_2b612bf20bda652dcb1a22d48081ea34< typename Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage>::FunctionSpaceType > >;
    Dune::Python::insertClass< Dune::Fem::SpaceOperatorInterface< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >>>( module, "cls0", Dune::Python::GenerateTypeName("Dune::Fem::SpaceOperatorInterface< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >>"), Dune::Python::IncludeFiles{"dune/fem-dg/operator/dg/dgpyoperator.hh","dune/fem-dg/operator/limiter/indicatorbase.hh","dune/fem-dg/python/operator.hh","dune/fem/function/adaptivefunction.hh","dune/fem/schemes/conservationlawmodel.hh","dune/fem/space/finitevolume.hh","dune/fempy/parameter.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/operator.hh","dune/fempy/py/space.hh","dune/grid/spgrid.hh","dune/grid/spgrid/dgfparser.hh","dune/python/grid/hierarchical.hh","python/dune/generated/modelFemDgAdv_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a.cc","python/dune/generated/modelFemDgDiff_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a.cc"});
    auto cls = Dune::Python::insertClass< DuneType, Dune::Fem::SpaceOperatorInterface< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >> >( cls0, "Operator",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::DGOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >, ModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > >, ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > >, Additional_2b612bf20bda652dcb1a22d48081ea34< typename Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage>::FunctionSpaceType > >"), Dune::Python::IncludeFiles{"dune/fem-dg/operator/dg/dgpyoperator.hh","dune/fem-dg/operator/limiter/indicatorbase.hh","dune/fem-dg/python/operator.hh","dune/fem/function/adaptivefunction.hh","dune/fem/schemes/conservationlawmodel.hh","dune/fem/space/finitevolume.hh","dune/fempy/parameter.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/operator.hh","dune/fempy/py/space.hh","dune/grid/spgrid.hh","dune/grid/spgrid/dgfparser.hh","dune/python/grid/hierarchical.hh","python/dune/generated/modelFemDgAdv_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a.cc","python/dune/generated/modelFemDgDiff_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a.cc"}).first;
    Dune::FemPy::registerOperator( cls0, cls );
    cls.def( "estimateMark", []( DuneType &self, const typename DuneType::DestinationType &u, const double dt) { self.estimateMark(u, dt); } );
    cls.def("setTroubledCellIndicator", [] ( DuneType &self, typename DuneType::TroubledCellIndicatorType indicator ) {self.setTroubledCellIndicator(indicator);
      } , pybind11::keep_alive<0,1>() );
    cls.def( "_indicator", &DuneType::indicator );
    cls.def( "info", &DuneType::counter );
  }
}
#endif
