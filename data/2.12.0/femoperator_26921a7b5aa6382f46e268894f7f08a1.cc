#ifndef Guard_femoperator_26921a7b5aa6382f46e268894f7f08a1
#define Guard_femoperator_26921a7b5aa6382f46e268894f7f08a1

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/schemes/galerkin.hh>
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

PYBIND11_MODULE( femoperator_26921a7b5aa6382f46e268894f7f08a1, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::DifferentiableGalerkinOperator< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage >::GridPartType, std::tuple< Dune::FieldMatrix< double, 3, 3 > >, std::tuple< Dune::FieldVector< double, 1 > > >, Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicRaviartThomasSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >,Dune::Fem::CodegenStorage> >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage > >,Dune::Fem::SparseRowMatrix<double,int>> >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Operator",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::DifferentiableGalerkinOperator< Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage >::GridPartType, std::tuple< Dune::FieldMatrix< double, 3, 3 > >, std::tuple< Dune::FieldVector< double, 1 > > >, Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicRaviartThomasSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >,Dune::Fem::CodegenStorage> >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage > >,Dune::Fem::SparseRowMatrix<double,int>> >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/spoperator.hh","dune/fem/schemes/galerkin.hh","dune/fem/space/hpdg/orthogonal.hh","dune/fem/space/raviartthomas.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/grid/gridpart.hh","dune/fempy/py/operator.hh","dune/fempy/py/space.hh","dune/python/alugrid/hierarchical.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::FemPy::registerOperator( cls0, cls );
    cls.def( "setCommunicate", []( DuneType &self, const bool communicate ) { self.setCommunicate( communicate ); } );
    cls.def( pybind11::init( [] ( const Dune::Fem::DynamicRaviartThomasSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >,Dune::Fem::CodegenStorage>& dSpace, const Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage > &rSpace, Dune::Fem::VirtualizedIntegrands< typename Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, 1 ,Dune::Fem::CodegenStorage >::GridPartType, std::tuple< Dune::FieldMatrix< double, 3, 3 > >, std::tuple< Dune::FieldVector< double, 1 > > > &integrands ) {return new DuneType( dSpace, rSpace, integrands );
      } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
  }
}
#endif
