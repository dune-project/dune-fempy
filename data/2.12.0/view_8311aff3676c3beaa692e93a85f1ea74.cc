#ifndef Guard_view_8311aff3676c3beaa692e93a85f1ea74
#define Guard_view_8311aff3676c3beaa692e93a85f1ea74

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/filter/simple.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE( view_8311aff3676c3beaa692e93a85f1ea74, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::Fem::FilteredGridPart< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::Fem::SimpleFilter< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > >, false >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "GridView",pybind11::dynamic_attr(), Dune::Python::GenerateTypeName("Dune::Fem::FilteredGridPart< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::Fem::SimpleFilter< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > >, false >"), Dune::Python::IncludeFiles{"dune/alugrid/dgf.hh","dune/alugrid/grid.hh","dune/fem/gridpart/filter/simple.hh","dune/fem/gridpart/filteredgridpart.hh","dune/fempy/py/gridview.hh","dune/python/alugrid/hierarchical.hh","dune/python/grid/gridview.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::Python::registerGridView( cls0, cls );

    Dune::FemPy::registerGridView ( cls );
    cls.def( pybind11::init( [] ( pybind11::handle hostGridView, pybind11::function contains, int domainId ) {auto containsCpp = [ contains ] ( const Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >::Codim< 0 >::EntityType &e ) {
            return contains( e ).template cast< int >();
          };
        Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > &hostGridPart = Dune::FemPy::gridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >( hostGridView );
        return Dune::Fem::FilteredGridPart< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::Fem::SimpleFilter< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > >, false > ( hostGridPart, Dune::Fem::SimpleFilter< Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > >( hostGridPart, containsCpp, domainId ) );
      } ), pybind11::keep_alive< 1, 2 >() );
  }
}
#endif
