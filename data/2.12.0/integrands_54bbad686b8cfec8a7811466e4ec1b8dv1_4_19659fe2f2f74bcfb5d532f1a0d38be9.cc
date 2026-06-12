#ifndef GuardIntegrands_54bbad686b8cfec8a7811466e4ec1b8dv1_4
#define GuardIntegrands_54bbad686b8cfec8a7811466e4ec1b8dv1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <cmath>
#include <dune/fem/misc/boundaryidprovider.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_54bbad686b8cfec8a7811466e4ec1b8dv1_4
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

    void decouple ()
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
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      const auto tmp0 = std::get< 1 >( u );
        assert(tmp0 == tmp0);
      return RangeValueType{ { 0 }, { { (tmp0[ 0 ])[ 0 ], (tmp0[ 0 ])[ 1 ] } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { ((std::get< 1 >( phi ))[ 0 ])[ 0 ], ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      const auto tmp0 = Dune::Fem::boundaryId<GridPartType>( (intersection_) );
        assert(tmp0 == tmp0);
      const auto tmp1 = -2 * (tmp0 == 3 ? 1.0 : 0.0);
        assert(tmp1 == tmp1);
      GlobalCoordinateType tmp2 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp2 == tmp2);
      const auto tmp3 = -0.5 + tmp2[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = std::abs( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = 10 * (tmp4 < 0.1 ? 1 : 0.0);
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = -1 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 * (tmp0 == 2 ? 1.0 : 0.0);
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp7 * (tmp0 == 1 ? 1.0 : 0.0);
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 + tmp8;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 + tmp1;
        assert(tmp11 == tmp11);
      return RangeValueType{ { tmp11 }, { { 0, 0 } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { 0, 0 } } };
        };
    }
    typedef Dune::FieldVector< double, 1 >  RRangeType;
    typedef Dune::FieldMatrix< double, 1, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,1> DirichletComponentType;

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
          int domainId;
          {
              using std::abs;
              using std::logical_and;
              auto tmp0 = intersection.geometry().center( );
                assert(tmp0 == tmp0);
              const auto tmp1 = -0.5 + tmp0[ 0 ];
                assert(tmp1 == tmp1);
              const auto tmp2 = std::abs( tmp1 );
                assert(tmp2 == tmp2);
              const auto tmp3 = Dune::Fem::boundaryId<GridPartType>( (intersection_) );
                assert(tmp3 == tmp3);
              const auto tmp4 = std::logical_and()( tmp3 == 4, tmp2 < 0.2 );
                assert(tmp4 == tmp4);
              domainId = tmp4;
              if (domainId)
              {
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 2 );
                return true;
              }
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
      case 2:
        {
          result[ 0 ] = 1;
        }
        break;
      default:
        {
          result = RRangeType( 0 );
        }
      }
    }

    template< class Point >
    void dDirichlet ( int bndId, const Point &x, RJacobianRangeType &result ) const
    {
      switch( bndId )
      {
      case 2:
        {
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = 0.0;
        }
        break;
      default:
        {
          result = RJacobianRangeType( 0 );
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

} // namespace Integrands_54bbad686b8cfec8a7811466e4ec1b8dv1_4

PYBIND11_MODULE( integrands_54bbad686b8cfec8a7811466e4ec1b8dv1_4_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_54bbad686b8cfec8a7811466e4ec1b8dv1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_54bbad686b8cfec8a7811466e4ec1b8dv1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_54bbad686b8cfec8a7811466e4ec1b8dv1_4_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_54bbad686b8cfec8a7811466e4ec1b8dv1_4
