#ifndef GuardIntegrands_12d76ab9bd466d5cfc504103195a7925v1_3
#define GuardIntegrands_12d76ab9bd466d5cfc504103195a7925v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_12d76ab9bd466d5cfc504103195a7925v1_3
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
    typedef std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      const auto tmp0 = std::get< 1 >( u );
      const auto tmp1 = (tmp0[ 0 ])[ 0 ] + (tmp0[ 1 ])[ 1 ];
      const auto tmp2 = std::get< 0 >( u );
      const auto tmp3 = -1 * tmp2[ 2 ];
      const auto tmp4 = (tmp0[ 0 ])[ 0 ] + tmp3;
      const auto tmp5 = (tmp0[ 1 ])[ 1 ] + tmp3;
      return RangeValueType{ { 0.0, 0.0, tmp1 }, { { tmp4, (tmp0[ 0 ])[ 1 ] }, { (tmp0[ 1 ])[ 0 ], tmp5 }, { 0.0, 0.0 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0, 0, ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + ((std::get< 1 >( phi ))[ 1 ])[ 1 ] }, { { -1 * (std::get< 0 >( phi ))[ 2 ] + ((std::get< 1 >( phi ))[ 0 ])[ 0 ], ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { ((std::get< 1 >( phi ))[ 1 ])[ 0 ], -1 * (std::get< 0 >( phi ))[ 2 ] + ((std::get< 1 >( phi ))[ 1 ])[ 1 ] }, { 0, 0 } } };
        };
    }
    typedef Dune::FieldVector< double, 3 >  RRangeType;
    typedef Dune::FieldMatrix< double, 3, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,3> DirichletComponentType;

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
            dirichletComponent[2] = 0;
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
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = -1 * tmp0[ 1 ];
          const auto tmp2 = 1 + tmp1;
          const auto tmp3 = tmp0[ 1 ] * tmp2;
          result[ 0 ] = tmp3;
          result[ 1 ] = 0.0;
          result[ 2 ] = 0.0;
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
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = -1 * tmp0[ 1 ];
          const auto tmp2 = 1 + tmp1;
          const auto tmp3 = tmp2 + tmp1;
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = tmp3;
          (result[ 1 ])[ 0 ] = 0.0;
          (result[ 1 ])[ 1 ] = 0.0;
          (result[ 2 ])[ 0 ] = 0.0;
          (result[ 2 ])[ 1 ] = 0.0;
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

} // namespace Integrands_12d76ab9bd466d5cfc504103195a7925v1_3

PYBIND11_MODULE( integrands_12d76ab9bd466d5cfc504103195a7925v1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_12d76ab9bd466d5cfc504103195a7925v1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_12d76ab9bd466d5cfc504103195a7925v1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_12d76ab9bd466d5cfc504103195a7925v1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_12d76ab9bd466d5cfc504103195a7925v1_3
