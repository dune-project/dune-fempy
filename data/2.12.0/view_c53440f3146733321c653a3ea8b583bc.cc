#ifndef Guard_view_c53440f3146733321c653a3ea8b583bc
#define Guard_view_c53440f3146733321c653a3ea8b583bc

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem/gridpart/geometrygridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/grid/gridview.hh>
#include <python/dune/generated/localfunction_b10ff5627bd2f22e898386281281f7db_e14d829faa8914e3ff98a63559ccacdd.cc>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( view_c53440f3146733321c653a3ea8b583bc, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::GeometryGridPart< UFLLocalFunctions_b10ff5627bd2f22e898386281281f7db::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >::LeafGridView > > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "GridView",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::GeometryGridPart< UFLLocalFunctions_b10ff5627bd2f22e898386281281f7db::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >::LeafGridView > > >"), Dune::Python::IncludeFiles{"dune/fem/gridpart/geometrygridpart.hh","dune/fempy/py/gridview.hh","dune/python/grid/gridview.hh","python/dune/generated/localfunction_b10ff5627bd2f22e898386281281f7db_e14d829faa8914e3ff98a63559ccacdd.cc"}).first;
    Dune::Python::registerGridView( cls0, cls );

    Dune::FemPy::registerGridView ( cls );
    cls.def( pybind11::init( [] ( UFLLocalFunctions_b10ff5627bd2f22e898386281281f7db::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >::LeafGridView > > &coordFunction ) {return Dune::Fem::GeometryGridPart< UFLLocalFunctions_b10ff5627bd2f22e898386281281f7db::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >::LeafGridView > > >( coordFunction );
      } ), pybind11::keep_alive< 1, 2 >() );
    cls.def( pybind11::pickle(
      [](const pybind11::object &self) {

                auto& gv = self.cast<DuneType&>();
                /* Return a tuple that fully encodes the state of the object */
                pybind11::dict d;
                if (pybind11::hasattr(self, "__dict__")) {
                  d = self.attr("__dict__");
                }
                return pybind11::make_tuple(gv.gridFunction(),d);
          }
    , [](pybind11::tuple t) {

                if (t.size() != 2)
                    throw std::runtime_error("Invalid state in GeoGV with "+std::to_string(t.size())+"arguments!");
                pybind11::handle pyGF = t[0];
                const auto& gf = pyGF.cast<const typename DuneType::GridFunctionType&>();
                /* Create a new C++ instance */
                auto py_state = t[1].cast<pybind11::dict>();
                return std::make_pair(std::make_unique<DuneType>(gf), py_state);
          }
    ), pybind11::keep_alive<1,2>(), pybind11::prepend() );
  }
}
#endif
