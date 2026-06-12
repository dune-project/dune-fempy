#ifndef Guard_femspace_90f0a9524a8cb701e8ee5027b7658a0e_ea2b37cdbbff85ca503c01ecc2747bb2
#define Guard_femspace_90f0a9524a8cb701e8ee5027b7658a0e_ea2b37cdbbff85ca503c01ecc2747bb2

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "as_istl"
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femspace_90f0a9524a8cb701e8ee5027b7658a0e_ea2b37cdbbff85ca503c01ecc2747bb2, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>"), Dune::Python::IncludeFiles{"dune/fem/function/blockvectorfunction.hh","dune/fem/space/lagrange.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >"), Dune::Python::IncludeFiles{"dune/fem/function/blockvectorfunction.hh","dune/fem/space/lagrange.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>&space, typename DuneType::DofContainerType &dofVector ) {return new DuneType(name,space,dofVector);
      } ), "name"_a, "space"_a, "dofVector"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
