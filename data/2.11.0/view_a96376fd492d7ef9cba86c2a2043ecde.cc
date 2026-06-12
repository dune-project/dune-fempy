#ifndef Guard_view_a96376fd492d7ef9cba86c2a2043ecde
#define Guard_view_a96376fd492d7ef9cba86c2a2043ecde

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/gridpart/geometrygridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( view_a96376fd492d7ef9cba86c2a2043ecde, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::GeometryGridPart< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "GridView",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::GeometryGridPart< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> > >"), Dune::Python::IncludeFiles{"dune/fem/function/blockvectorfunction.hh","dune/fem/gridpart/geometrygridpart.hh","dune/fem/space/lagrange.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/gridview.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/gridview.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::Python::registerGridView( cls0, cls );

    Dune::FemPy::registerGridView ( cls );
    cls.def( pybind11::init( [] ( Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> > &coordFunction ) {return Dune::Fem::GeometryGridPart< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> > >( coordFunction );
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
