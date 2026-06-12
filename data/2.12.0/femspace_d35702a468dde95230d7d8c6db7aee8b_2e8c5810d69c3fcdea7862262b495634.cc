#ifndef Guard_femspace_d35702a468dde95230d7d8c6db7aee8b_2e8c5810d69c3fcdea7862262b495634
#define Guard_femspace_d35702a468dde95230d7d8c6db7aee8b_2e8c5810d69c3fcdea7862262b495634

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "as_petsc"
#include </home/dedner/dune-core/lib/python3.12/site-packages/petsc4py/include/petsc4py/petsc4py.h>
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femspace_d35702a468dde95230d7d8c6db7aee8b_2e8c5810d69c3fcdea7862262b495634, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>"), Dune::Python::IncludeFiles{"/home/dedner/dune-core/lib/python3.12/site-packages/petsc4py/include/petsc4py/petsc4py.h","dune/fem/function/petscdiscretefunction.hh","dune/fem/space/lagrange.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::PetscDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::PetscDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >"), Dune::Python::IncludeFiles{"/home/dedner/dune-core/lib/python3.12/site-packages/petsc4py/include/petsc4py/petsc4py.h","dune/fem/function/petscdiscretefunction.hh","dune/fem/space/lagrange.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>&space, pybind11::handle dofVector ) {if (import_petsc4py() != 0) {
          throw std::runtime_error("Error during import of petsc4py");
        }
        Vec petscVec = PyPetscVec_Get(dofVector.ptr());
        typename DuneType::DofVectorType *dofStorage = new typename DuneType::DofVectorType(space,petscVec);
        pybind11::cpp_function remove_petscStorage( [ dofStorage, dofVector, petscVec] ( pybind11::handle weakref ) {
          delete dofStorage;
          weakref.dec_ref();
        } );
        pybind11::weakref weakref( dofVector, remove_petscStorage );
        weakref.release();
        return new DuneType( name, space, *dofStorage );
      } ), "name"_a, "space"_a, "dofVector"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
