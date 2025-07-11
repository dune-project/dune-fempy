#ifndef GuardIntegrands_7d578911193cd1d10bd8f77fd79f4a8bv1_3
#define GuardIntegrands_7d578911193cd1d10bd8f77fd79f4a8bv1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_7d578911193cd1d10bd8f77fd79f4a8bv1_3
{




  // Integrands
  // ----------

  template< class GridPart >
  struct Integrands
  {
    typedef GridPart GridPartType;
    typedef typename GridPartType::GridViewType GridView;
    typedef typename GridView::ctype ctype;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename GridPartType::IntersectionType IntersectionType;
    typedef typename EntityType::Geometry Geometry;
    typedef typename Geometry::GlobalCoordinate GlobalCoordinateType;
    typedef Dune::Fem::IntersectionSide Side;
    typedef std::tuple<> ConstantTupleType;
    typedef std::tuple<> CoefficientTupleType;
    static constexpr bool gridPartValid = true;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {}

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      geometry_.emplace( this->entity().geometry() );
      return true;
    }

    void unbind ()
    {}

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      return (intersection.boundary() && init( intersection.inside() ));
    }
    typedef std::tuple< Dune::FieldVector< double, 1 > > DomainValueType;
    typedef std::tuple< Dune::FieldMatrix< double, 2, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      const auto tmp0 = std::get< 0 >( u );
      return RangeValueType{ { { tmp0[ 0 ], 0.0 }, { 0.0, tmp0[ 0 ] } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { { (std::get< 0 >( phi ))[ 0 ], 0 }, { 0, (std::get< 0 >( phi ))[ 0 ] } } };
        };
    }
    typedef Dune::FieldVector< double, 2 >  RRangeType;
    typedef Dune::FieldMatrix< double, 2, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,2> DirichletComponentType;

    bool hasDirichletBoundary () const
    {
      return true;
    }

    bool isDirichletIntersection ( const IntersectionType &intersection, DirichletComponentType &dirichletComponent ) const
    {
      const int bndId = BoundaryIdProviderType::boundaryId( intersection );
      {
        std::fill( dirichletComponent.begin(), dirichletComponent.end(), bndId );
      }
      switch( bndId )
      {
      default:
        {
          {
            std::fill( dirichletComponent.begin(), dirichletComponent.end(), 1 );
            return true;
          }
          return false;
        }
      }
    }

    template< class Point >
    void dirichlet ( int bndId, const Point &x, RRangeType &result ) const
    {
      switch( bndId )
      {
      default:
        {
          using std::cos;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
          const auto tmp2 = std::cos( tmp1 );
          const auto tmp3 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp4 = std::sin( tmp3 );
          const auto tmp5 = -1 * tmp4;
          const auto tmp6 = 3.141592653589793 * tmp5;
          const auto tmp7 = tmp6 * tmp2;
          result[ 0 ] = tmp7;
        }
      }
    }

    template< class Point >
    void dDirichlet ( int bndId, const Point &x, RJacobianRangeType &result ) const
    {
      switch( bndId )
      {
      default:
        {
          using std::cos;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
          const auto tmp2 = std::cos( tmp1 );
          const auto tmp3 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp4 = std::cos( tmp3 );
          const auto tmp5 = 3.141592653589793 * tmp4;
          const auto tmp6 = -1 * tmp5;
          const auto tmp7 = 3.141592653589793 * tmp6;
          const auto tmp8 = tmp7 * tmp2;
          const auto tmp9 = std::sin( tmp1 );
          const auto tmp10 = -1 * tmp9;
          const auto tmp11 = 3.141592653589793 * tmp10;
          const auto tmp12 = std::sin( tmp3 );
          const auto tmp13 = -1 * tmp12;
          const auto tmp14 = 3.141592653589793 * tmp13;
          const auto tmp15 = tmp14 * tmp11;
          (result[ 0 ])[ 0 ] = tmp8;
          (result[ 0 ])[ 1 ] = tmp15;
        }
      }
    }

    template< std::size_t i >
    const ConstantType< i > &constant () const
    {
      return *std::get< i >( constants_ );
    }

    template< std::size_t i >
    ConstantType< i > &constant ()
    {
      return *std::get< i >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_;
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }
    EntityType entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    ConstantTupleType constants_;
    std::tuple<  > coefficients_;
  };

} // namespace Integrands_7d578911193cd1d10bd8f77fd79f4a8bv1_3

PYBIND11_MODULE( integrands_7d578911193cd1d10bd8f77fd79f4a8bv1_3_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > GridPart;
  typedef Integrands_7d578911193cd1d10bd8f77fd79f4a8bv1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_7d578911193cd1d10bd8f77fd79f4a8bv1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_7d578911193cd1d10bd8f77fd79f4a8bv1_3_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_7d578911193cd1d10bd8f77fd79f4a8bv1_3
