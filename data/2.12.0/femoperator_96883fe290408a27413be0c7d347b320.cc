#ifndef Guard_femoperator_96883fe290408a27413be0c7d347b320
#define Guard_femoperator_96883fe290408a27413be0c7d347b320

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/space/hpdg/orthogonal.hh>
#include <dune/fem/space/raviartthomas.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/operator.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femoperator_96883fe290408a27413be0c7d347b320, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicRaviartThomasSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >,Dune::Fem::CodegenStorage> >,Dune::Fem::SparseRowMatrix<double,int>>;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "LinearOperator",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicRaviartThomasSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >,Dune::Fem::CodegenStorage> >,Dune::Fem::SparseRowMatrix<double,int>>"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/spoperator.hh","dune/fem/space/hpdg/orthogonal.hh","dune/fem/space/raviartthomas.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/grid/gridpart.hh","dune/fempy/py/operator.hh","dune/fempy/py/space.hh","dune/python/alugrid/hierarchical.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerLinearOperator( cls0, cls );
    cls.def( pybind11::init( [] ( const Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage >& dSpace, const Dune::Fem::DynamicRaviartThomasSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >,Dune::Fem::CodegenStorage> &rSpace, const pybind11::dict &parameters ) {return new Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicRaviartThomasSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >,Dune::Fem::CodegenStorage> >,Dune::Fem::SparseRowMatrix<double,int>>( "tmp", dSpace, rSpace, Dune::FemPy::pyParameter( "fem.solver.", parameters, std::make_shared< std::string >() ));
      } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
