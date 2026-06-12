#ifndef Guard_femspace_5f752739f11fd9d7796acc936bd0f32a_3ceaa7153c8b876e2ec9f9e52c46e983
#define Guard_femspace_5f752739f11fd9d7796acc936bd0f32a_3ceaa7153c8b876e2ec9f9e52c46e983

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "as_numpy"
#include <dune/common/quadmath.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>
#include <dune/vem/space/hk.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( femspace_5f752739f11fd9d7796acc936bd0f32a_3ceaa7153c8b876e2ec9f9e52c46e983, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double >"), Dune::Python::IncludeFiles{"dune/common/quadmath.hh","dune/fem/function/adaptivefunction.hh","dune/fem/gridpart/common/gridpart.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/misc/grid.hh","dune/vem/space/hk.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
    cls.def( pybind11::init( [] ( pybind11::handle gridView, const pybind11::object agglomerate, unsigned int order, const std::array<unsigned int,3> orderTuple, const std::array<std::vector<int>,3> &testSpaces, int basisChoice, bool edgeInterpolation, bool rotatedBB ) {typedef Dune::Vem::Agglomeration<Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >> AggloType;
        AggloType *agglo = nullptr;
        // check storage
        pybind11::function retrieve = pybind11::module::import("dune.vem.vem").attr("retrieveAgglo_");
        pybind11::function insert   = pybind11::module::import("dune.vem.vem").attr("insertAgglo_");
        pybind11::function remove   = pybind11::module::import("dune.vem.vem").attr("removeAgglo_");
        pybind11::handle aggloExists = retrieve(gridView);
        if (! pybind11::isinstance<pybind11::none>(aggloExists) ) {
          // std::cout << "retrieve agglo\n";
          agglo = static_cast<AggloType*>(aggloExists.cast<void*>());
        }
        else {
          agglo = new AggloType(
                   Dune::FemPy::gridPart<Dune::Vem::Grid<2,2>::LeafGridView>(gridView),rotatedBB,
              [agglomerate](const auto& idx)
              { return agglomerate(idx).template cast<unsigned int>(); } );
          insert(gridView,(void*)agglo);
          pybind11::cpp_function aggloCleanup(
                [agglo,remove](pybind11::handle weakref) {
                // std::cout << "remove agglo\n";
                remove((void*)agglo);
                delete agglo;
                weakref.dec_ref();
            });
          (void) pybind11::weakref(gridView, aggloCleanup).release();
          // std::cout << "new agglo\n";
        }
        auto obj = new DuneType(*agglo, order, orderTuple, testSpaces, basisChoice, edgeInterpolation );
        return obj;
      } ), "gridView"_a, "agglomerate"_a, "order"_a, "orderTuple"_a, "testSpaces"_a, "basisChoice"_a, "edgeInterpolation"_a, "rotatedBB"_a, pybind11::keep_alive< 1, 2 >() );
    cls.def( "diameters", []( DuneType &self ) { return self.blockMapper().indexSet().diameters(); } );
    cls.def( "update", []( DuneType &self ) { self.update(); } );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double > >"), Dune::Python::IncludeFiles{"dune/common/quadmath.hh","dune/fem/function/adaptivefunction.hh","dune/fem/gridpart/common/gridpart.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/misc/grid.hh","dune/vem/space/hk.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Vem::AgglomerationVEMSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, false , false , double, long double >&space, pybind11::array_t<double> dofVector ) {double *dof = static_cast< double* >( dofVector.request(false).ptr );
        return new DuneType(name,space,dof);
      } ), "name"_a, "space"_a, "dofVector"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
