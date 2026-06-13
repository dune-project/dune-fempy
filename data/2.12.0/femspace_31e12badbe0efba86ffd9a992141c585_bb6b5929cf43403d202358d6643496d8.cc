#ifndef Guard_femspace_31e12badbe0efba86ffd9a992141c585_bb6b5929cf43403d202358d6643496d8
#define Guard_femspace_31e12badbe0efba86ffd9a992141c585_bb6b5929cf43403d202358d6643496d8

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "as_numpy"
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femspace_31e12badbe0efba86ffd9a992141c585_bb6b5929cf43403d202358d6643496d8, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage>;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage>"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/finitevolume.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/spgrid.hh","dune/grid/spgrid/dgfparser.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage> >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/finitevolume.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/spgrid.hh","dune/grid/spgrid/dgfparser.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Fem::FiniteVolumeSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, 0 , Dune::Fem::SimpleStorage>&space, pybind11::array_t<double> dofVector ) {double *dof = static_cast< double* >( dofVector.request(false).ptr );
        return new DuneType(name,space,dof);
      } ), "name"_a, "space"_a, "dofVector"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif