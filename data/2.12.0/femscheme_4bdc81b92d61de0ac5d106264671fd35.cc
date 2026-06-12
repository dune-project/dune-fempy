#ifndef Guard_femscheme_4bdc81b92d61de0ac5d106264671fd35
#define Guard_femscheme_4bdc81b92d61de0ac5d106264671fd35

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/schemes/conservationlawmodel.hh>
#include <dune/fem/schemes/dgelliptic.hh>
#include <dune/fem/schemes/femscheme.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/scheme.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/agglomeration/dgspace.hh>
#include <dune/vem/misc/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( femscheme_4bdc81b92d61de0ac5d106264671fd35, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::FemScheme< DifferentiableDGEllipticOperator< Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 > >,Dune::Fem::SparseRowMatrix<double,int>>,DGConservationLawModel< typename Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >::GridPartType, Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >::dimRange, Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >::dimRange, typename Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >::RangeFieldType >,Dune::Vem::BBDGPenalty<Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >>>, Dune::Fem::KrylovInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 > > > >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "Scheme",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::FemScheme< DifferentiableDGEllipticOperator< Dune::Fem::SparseRowLinearOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 > >,Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 > >,Dune::Fem::SparseRowMatrix<double,int>>,DGConservationLawModel< typename Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >::GridPartType, Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >::dimRange, Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >::dimRange, typename Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >::RangeFieldType >,Dune::Vem::BBDGPenalty<Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 >>>, Dune::Fem::KrylovInverseOperator< Dune::Fem::AdaptiveDiscreteFunction< Dune::Vem::AgglomerationDGSpace< Dune::Fem::FunctionSpace< double, double, 2, 1 >, Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, 1 > > > >"), Dune::Python::IncludeFiles{"dune/fem/function/adaptivefunction.hh","dune/fem/operator/linear/spoperator.hh","dune/fem/schemes/conservationlawmodel.hh","dune/fem/schemes/dgelliptic.hh","dune/fem/schemes/femscheme.hh","dune/fem/solver/krylovinverseoperators.hh","dune/fempy/parameter.hh","dune/fempy/py/discretefunction.hh","dune/fempy/py/scheme.hh","dune/fempy/py/space.hh","dune/python/grid/hierarchical.hh","dune/vem/agglomeration/dgspace.hh","dune/vem/misc/grid.hh"}).first;
    Dune::FemPy::registerScheme( cls0, cls );
  }
}
#endif
