#ifndef Guard_femscheme_6c4cb833cbf91db26b1c470fce859cda
#define Guard_femscheme_6c4cb833cbf91db26b1c470fce859cda

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/schemes/dirichletwrapper.hh>
#include <dune/fem/schemes/linearized.hh>
#include <dune/fem/schemes/molgalerkin.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/scheme.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( femscheme_6c4cb833cbf91db26b1c470fce859cda, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::LinearizedScheme< Dune::Fem::GalerkinScheme< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >, Dune::Fem::ISTLLinearOperator< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >,Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >>, Dune::Fem::ISTLInverseOperator< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> > >, true > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Scheme",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::LinearizedScheme< Dune::Fem::GalerkinScheme< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage>::GridPartType, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > >, std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > >, Dune::Fem::ISTLLinearOperator< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >,Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> >>, Dune::Fem::ISTLInverseOperator< Dune::Fem::ISTLBlockVectorDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Fem::CodegenStorage> > >, true > >"), Dune::Python::IncludeFiles{"dune/fem/function/blockvectorfunction.hh","dune/fem/operator/linear/istloperator.hh","dune/fem/schemes/dirichletwrapper.hh","dune/fem/schemes/linearized.hh","dune/fem/schemes/molgalerkin.hh","dune/fem/solver/istlinverseoperators.hh","dune/fem/space/lagrange.hh","dune/fempy/parameter.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/scheme.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerScheme( cls0, cls );
    cls.def( pybind11::init( [] ( typename DuneType::SchemeType &scheme, typename DuneType::DiscreteFunctionType &ubar, const pybind11::dict &parameters ) {return new DuneType( scheme, ubar, Dune::FemPy::pyParameter( "fem.solver.", parameters, std::make_shared< std::string >() ) );
      } ), "scheme"_a, "ubar"_a, "parameters"_a, pybind11::keep_alive< 1, 2 >() );
    cls.def( pybind11::init( [] ( typename DuneType::SchemeType &scheme, const pybind11::dict &parameters ) {return new DuneType( scheme, Dune::FemPy::pyParameter( "fem.solver.", parameters, std::make_shared< std::string >() ) );
      } ), "scheme"_a, "parameters"_a, pybind11::keep_alive< 1, 2 >() );
    cls.def( "setup", static_cast<void (DuneType::*)(const typename DuneType::DiscreteFunctionType &)>(&DuneType::setup) );
    cls.def( "setup", static_cast<void (DuneType::*)()>(&DuneType::setup) );
    cls.def( "linearScheme", []( const DuneType &self) { return pybind11::cast(self.linearScheme(), pybind11::return_value_policy::reference ); } );
  }
}
#endif