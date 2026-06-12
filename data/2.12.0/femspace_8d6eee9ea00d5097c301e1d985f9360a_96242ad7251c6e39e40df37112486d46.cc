#ifndef Guard_femspace_8d6eee9ea00d5097c301e1d985f9360a_96242ad7251c6e39e40df37112486d46
#define Guard_femspace_8d6eee9ea00d5097c301e1d985f9360a_96242ad7251c6e39e40df37112486d46

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "as_numpy"
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/agglomeration/dgspace.hh>
#include <dune/vem/misc/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femspace_8d6eee9ea00d5097c301e1d985f9360a_96242ad7251c6e39e40df37112486d46, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/agglomeration/dgspace.hh","dune/vem/misc/grid.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
    cls.def( pybind11::init( [] ( pybind11::handle gridView, const pybind11::object agglomerate, bool rotatedBB ) {typedef Dune::Vem::Agglomeration<Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >> AggloType;
        AggloType *agglo = nullptr;
        // check storage
        pybind11::function retrieve = pybind11::module::import("dune.vem.vem").attr("retrieveAgglo_");
        pybind11::function insert   = pybind11::module::import("dune.vem.vem").attr("insertAgglo_");
        pybind11::function remove   = pybind11::module::import("dune.vem.vem").attr("removeAgglo_");
        pybind11::handle aggloExists = retrieve(gridView);
        if (! pybind11::isinstance<pybind11::none>(aggloExists) ) {
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
                remove((void*)agglo);
                delete agglo;
                weakref.dec_ref();
            });
          (void) pybind11::weakref(gridView, aggloCleanup).release();
        }
        auto obj = new DuneType( *agglo );
        return obj;
      } ), "gridView"_a, "agglomerate"_a, "rotatedBB"_a, pybind11::keep_alive< 1, 2 >() );
    cls.def( "update", []( DuneType &self ) { self.update(); } );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 > >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/agglomeration/dgspace.hh","dune/vem/misc/grid.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >&space, pybind11::array_t<double> dofVector ) {double *dof = static_cast< double* >( dofVector.request(false).ptr );
        return new DuneType(name,space,dof);
      } ), "name"_a, "space"_a, "dofVector"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
