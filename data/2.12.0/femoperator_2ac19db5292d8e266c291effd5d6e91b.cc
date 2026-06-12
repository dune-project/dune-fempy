#ifndef Guard_femoperator_2ac19db5292d8e266c291effd5d6e91b
#define Guard_femoperator_2ac19db5292d8e266c291effd5d6e91b

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/space/brezzidouglasmarini.hh>
#include <dune/fem/space/hpdg/orthogonal.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/operator.hh>
#include <dune/fempy/py/space.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femoperator_2ac19db5292d8e266c291effd5d6e91b, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 4 ,Dune::Fem::CodegenStorage > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicBrezziDouglasMariniSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >,Dune::Fem::CodegenStorage> >,Dune::Fem::SparseRowMatrix<double,int>>;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "LinearOperator",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 4 ,Dune::Fem::CodegenStorage > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicBrezziDouglasMariniSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >,Dune::Fem::CodegenStorage> >,Dune::Fem::SparseRowMatrix<double,int>>"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/spoperator.hh","dune/fem/space/brezzidouglasmarini.hh","dune/fem/space/hpdg/orthogonal.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/grid/gridpart.hh","dune/fempy/py/operator.hh","dune/fempy/py/space.hh","dune/grid/io/file/dgfparser/dgfyasp.hh","dune/grid/yaspgrid.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerLinearOperator( cls0, cls );
    cls.def( pybind11::init( [] ( const Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 4 ,Dune::Fem::CodegenStorage >& dSpace, const Dune::Fem::DynamicBrezziDouglasMariniSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >,Dune::Fem::CodegenStorage> &rSpace, const pybind11::dict &parameters ) {return new Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, 4 ,Dune::Fem::CodegenStorage > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicBrezziDouglasMariniSpace< Dune::Fem::FunctionSpace< double, double, 2, 2 >, Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >,Dune::Fem::CodegenStorage> >,Dune::Fem::SparseRowMatrix<double,int>>( "tmp", dSpace, rSpace, Dune::FemPy::pyParameter( "fem.solver.", parameters, std::make_shared< std::string >() ));
      } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
