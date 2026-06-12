#ifndef Guard_femspace_08eac4afbccc577e8e400c93df95836b_82ccbff1620d47da982577e64177ece9
#define Guard_femspace_08eac4afbccc577e8e400c93df95836b_82ccbff1620d47da982577e64177ece9

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "as_numpy"
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/hpdg/legendre.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femspace_08eac4afbccc577e8e400c93df95836b_82ccbff1620d47da982577e64177ece9, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView >, 1, Dune::Fem::CodegenStorage>;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView >, 1, Dune::Fem::CodegenStorage>"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/hpdg/legendre.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView >, 1, Dune::Fem::CodegenStorage> >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView >, 1, Dune::Fem::CodegenStorage> >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/hpdg/legendre.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView >, 1, Dune::Fem::CodegenStorage>&space, pybind11::array_t<double> dofVector ) {double *dof = static_cast< double* >( dofVector.request(false).ptr );
        return new DuneType(name,space,dof);
      } ), "name"_a, "space"_a, "dofVector"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
