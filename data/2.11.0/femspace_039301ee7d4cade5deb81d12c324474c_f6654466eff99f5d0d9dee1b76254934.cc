#ifndef Guard_femspace_039301ee7d4cade5deb81d12c324474c_f6654466eff99f5d0d9dee1b76254934
#define Guard_femspace_039301ee7d4cade5deb81d12c324474c_f6654466eff99f5d0d9dee1b76254934

#include <config.h>

#define USING_DUNE_PYTHON 1

#define BACKENDNAME "None"
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/combinedspace/tuplespace.hh>
#include <dune/fem/space/hpdg/legendre.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( femspace_039301ee7d4cade5deb81d12c324474c_f6654466eff99f5d0d9dee1b76254934, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage>, Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >;
    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "Space",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage>, Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/fem/function/adaptivefunction.hh","dune/fem/function/tuplediscretefunction.hh","dune/fem/gridpart/adaptiveleafgridpart.hh","dune/fem/space/combinedspace/tuplespace.hh","dune/fem/space/hpdg/legendre.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/gridview.hh","dune/fempy/py/space.hh","dune/python/grid/gridview.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerSpace( cls0, cls );
    cls.def( pybind11::init( [] ( typename DuneType::DiscreteFunctionSpaceTupleType spaceTuple ) {return new DuneType( spaceTuple);
      } ), "spaceTuple"_a, pybind11::keep_alive<1,2>() );
  }
  using pybind11::operator""_a;
  {
    using DuneType = Dune::Fem::TupleDiscreteFunction< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "DiscreteFunction",pybind11::buffer_protocol(),pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::TupleDiscreteFunction< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >, Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> > >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/fem/function/adaptivefunction.hh","dune/fem/function/tuplediscretefunction.hh","dune/fem/gridpart/adaptiveleafgridpart.hh","dune/fem/space/combinedspace/tuplespace.hh","dune/fem/space/hpdg/legendre.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/gridview.hh","dune/fempy/py/space.hh","dune/python/grid/gridview.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerDiscreteFunction( cls0, cls );
    cls.def( pybind11::init( [] ( const std::string &name, const Dune::Fem::TupleDiscreteFunctionSpace< Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage>, Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, 3, Dune::Fem::CodegenStorage> >&space, pybind11::handle vec ) {std::cerr <<"Can not use constructor with dof vector argument for this type of discrete function storage!\n";
        throw std::runtime_error("Can not use constructor with dof vector argument for this type of discrete function storage!");
        return new DuneType(name,space);
      } ), "name"_a, "space"_a, "vec"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
