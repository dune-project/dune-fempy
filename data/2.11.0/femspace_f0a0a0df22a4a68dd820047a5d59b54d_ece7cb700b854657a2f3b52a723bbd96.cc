#ifndef Guard_femspace_f0a0a0df22a4a68dd820047a5d59b54d_ece7cb700b854657a2f3b52a723bbd96
#define Guard_femspace_f0a0a0df22a4a68dd820047a5d59b54d_ece7cb700b854657a2f3b52a723bbd96

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "as_numpy"
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/polygongrid/dgf.hh>
#include <dune/polygongrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( femspace_f0a0a0df22a4a68dd820047a5d59b54d_ece7cb700b854657a2f3b52a723bbd96, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 4 >, Dune::FemPy::GridPart< Dune::PolygonGrid< double >::LeafGridView >, 1, Dune::Fem::CodegenStorage >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 4 >, Dune::FemPy::GridPart< Dune::PolygonGrid< double >::LeafGridView >, 1, Dune::Fem::CodegenStorage >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/discontinuousgalerkin.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/polygongrid/dgf.hh","dune/polygongrid/grid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 4 >, Dune::FemPy::GridPart< Dune::PolygonGrid< double >::LeafGridView >, 1, Dune::Fem::CodegenStorage > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 4 >, Dune::FemPy::GridPart< Dune::PolygonGrid< double >::LeafGridView >, 1, Dune::Fem::CodegenStorage > >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/discontinuousgalerkin.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/polygongrid/dgf.hh","dune/polygongrid/grid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 4 >, Dune::FemPy::GridPart< Dune::PolygonGrid< double >::LeafGridView >, 1, Dune::Fem::CodegenStorage >&space, pybind11::array_t<double> dofVector ) {double *dof = static_cast< double* >( dofVector.request(false).ptr );
        return new DuneType(name,space,dof);
      } ), "name"_a, "space"_a, "dofVector"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif