#ifndef Guard_view_b6e87c20f8b7f70aa8b9315ea4255bb6
#define Guard_view_b6e87c20f8b7f70aa8b9315ea4255bb6

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem/gridpart/geometrygridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/grid/gridview.hh>
#include <python/dune/generated/localfunction_7bbf374aa1348eddf7b1e4fbb1008ee2_19659fe2f2f74bcfb5d532f1a0d38be9.cc>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( view_b6e87c20f8b7f70aa8b9315ea4255bb6, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::GeometryGridPart< UFLLocalFunctions_7bbf374aa1348eddf7b1e4fbb1008ee2::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "GridView",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::GeometryGridPart< UFLLocalFunctions_7bbf374aa1348eddf7b1e4fbb1008ee2::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > >"), Dune::Python::IncludeFiles{"dune/fem/gridpart/geometrygridpart.hh","dune/fempy/py/gridview.hh","dune/python/grid/gridview.hh","python/dune/generated/localfunction_7bbf374aa1348eddf7b1e4fbb1008ee2_19659fe2f2f74bcfb5d532f1a0d38be9.cc"}).first;
    Dune::Python::registerGridView( cls0, cls );

    Dune::FemPy::registerGridView ( cls );
    cls.def( pybind11::init( [] ( UFLLocalFunctions_7bbf374aa1348eddf7b1e4fbb1008ee2::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > &coordFunction ) {return Dune::Fem::GeometryGridPart< UFLLocalFunctions_7bbf374aa1348eddf7b1e4fbb1008ee2::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > >( coordFunction );
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
