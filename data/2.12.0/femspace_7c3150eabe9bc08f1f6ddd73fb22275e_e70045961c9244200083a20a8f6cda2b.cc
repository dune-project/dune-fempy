#ifndef Guard_femspace_7c3150eabe9bc08f1f6ddd73fb22275e_e70045961c9244200083a20a8f6cda2b
#define Guard_femspace_7c3150eabe9bc08f1f6ddd73fb22275e_e70045961c9244200083a20a8f6cda2b

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "as_numpy"
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femspace_7c3150eabe9bc08f1f6ddd73fb22275e_e70045961c9244200083a20a8f6cda2b, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1, Dune::Fem::CodegenStorage >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1, Dune::Fem::CodegenStorage >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/discontinuousgalerkin.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/misc/grid.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1, Dune::Fem::CodegenStorage > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1, Dune::Fem::CodegenStorage > >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/discontinuousgalerkin.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/misc/grid.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Fem::DiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1, Dune::Fem::CodegenStorage >&space, pybind11::array_t<double> dofVector ) {double *dof = static_cast< double* >( dofVector.request(false).ptr );
        return new DuneType(name,space,dof);
      } ), "name"_a, "space"_a, "dofVector"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
