#ifndef Guard_femspace_6727914baecf4df5f1fad3dd151bd363_e2bb3dac42f9e6d7ba50a2f381192e4e
#define Guard_femspace_6727914baecf4df5f1fad3dd151bd363_e2bb3dac42f9e6d7ba50a2f381192e4e

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "as_numpy"
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/combinedspace/tuplespace.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femspace_6727914baecf4df5f1fad3dd151bd363_e2bb3dac42f9e6d7ba50a2f381192e4e, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/combinedspace/tuplespace.hh","dune/fem/space/discontinuousgalerkin.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
    cls.def( pybind11::init( [] ( typename DuneType::DiscreteFunctionSpaceTupleType spaceTuple ) {return new DuneType( spaceTuple);
      } ), "spaceTuple"_a, pybind11::keep_alive<1,2>() );
    cls.def( pybind11::pickle(
      [](const pybind11::object &self) {

                auto& tsp = self.cast<DuneType&>();
                /* Return a tuple that fully encodes the state of the object */
                pybind11::dict d;
                if (pybind11::hasattr(self, "__dict__")) {
                  d = self.attr("__dict__");
                }
                return pybind11::make_tuple(tsp.spaceTuple(),d);
          }
    , [](pybind11::tuple t) {

                if (t.size() != 2)
                    throw std::runtime_error("Invalid state in AdaptGV with "+std::to_string(t.size())+"arguments!");
                pybind11::handle pyspaceTuple = t[0];
                auto spaceTuple = pyspaceTuple.cast<typename DuneType::DiscreteFunctionSpaceTupleType>();
                /* Create a new C++ instance */
                DuneType* tsp = new DuneType(spaceTuple);
                auto py_state = t[1].cast<pybind11::dict>();
                return std::make_pair(tsp, py_state);
          }
    ), pybind11::keep_alive<1,2>(), pybind11::prepend() );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> > >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/space/combinedspace/tuplespace.hh","dune/fem/space/discontinuousgalerkin.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 3 ,Dune::Fem::CodegenStorage>, Dune::Fem::LagrangeDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 2 ,Dune::Fem::CodegenStorage> >&space, pybind11::array_t<double> dofVector ) {double *dof = static_cast< double* >( dofVector.request(false).ptr );
        return new DuneType(name,space,dof);
      } ), "name"_a, "space"_a, "dofVector"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
