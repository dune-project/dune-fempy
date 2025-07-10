#ifndef GuardIntegrands_ddf3ff4d527ccaee8da94c708bfe32f2v1_3
#define GuardIntegrands_ddf3ff4d527ccaee8da94c708bfe32f2v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>
#include <cmath>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/vem/py/integrands.hh>

namespace Integrands_ddf3ff4d527ccaee8da94c708bfe32f2v1_3
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
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > RangeValueType;
    static constexpr bool _nonlinear = true;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
      const auto tmp2 = -1 * tmp1;
      const auto tmp3 = tmp0[ 1 ] + tmp2;
      const auto tmp4 = tmp0[ 0 ] + tmp0[ 0 ];
      const auto tmp5 = -1 * tmp4;
      const auto tmp6 = 1 + tmp5;
      const auto tmp7 = tmp6 * tmp3;
      const auto tmp8 = 2 * tmp7;
      const auto tmp9 = tmp0[ 0 ] * tmp0[ 0 ];
      const auto tmp10 = -1 * tmp9;
      const auto tmp11 = tmp0[ 0 ] + tmp10;
      const auto tmp12 = tmp11 * tmp3;
      const auto tmp13 = tmp12 * tmp8;
      const auto tmp14 = tmp7 * tmp13;
      const auto tmp15 = -2 * tmp3;
      const auto tmp16 = std::pow( tmp12, 2 );
      const auto tmp17 = 1.0 + tmp16;
      const auto tmp18 = tmp17 * tmp15;
      const auto tmp19 = tmp18 + tmp14;
      const auto tmp20 = tmp0[ 1 ] + tmp0[ 1 ];
      const auto tmp21 = -1 * tmp20;
      const auto tmp22 = 1 + tmp21;
      const auto tmp23 = tmp22 * tmp11;
      const auto tmp24 = 2 * tmp23;
      const auto tmp25 = tmp12 * tmp24;
      const auto tmp26 = tmp23 * tmp25;
      const auto tmp27 = -2 * tmp11;
      const auto tmp28 = tmp17 * tmp27;
      const auto tmp29 = tmp28 + tmp26;
      const auto tmp30 = tmp29 + tmp19;
      const auto tmp31 = -1 * tmp30;
      const auto tmp32 = -1 * tmp31;
      const auto tmp33 = std::get< 0 >( u );
      const auto tmp34 = std::pow( tmp33[ 0 ], 2 );
      const auto tmp35 = 1.0 + tmp34;
      const auto tmp36 = std::get< 1 >( u );
      const auto tmp37 = (tmp36[ 0 ])[ 0 ] * tmp35;
      const auto tmp38 = (tmp36[ 0 ])[ 1 ] * tmp35;
      return RangeValueType{ { tmp32 }, { { tmp37, tmp38 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      const auto tmp0 = std::get< 0 >( u );
      const auto tmp1 = 2 * tmp0[ 0 ];
      const auto tmp2 = std::get< 1 >( u );
      const auto tmp3 = (tmp2[ 0 ])[ 0 ] * tmp1;
      const auto tmp4 = (tmp2[ 0 ])[ 1 ] * tmp1;
      const auto tmp5 = std::pow( tmp0[ 0 ], 2 );
      const auto tmp6 = 1.0 + tmp5;
      return [ tmp3, tmp4, tmp6 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { tmp3 * (std::get< 0 >( phi ))[ 0 ] + tmp6 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], tmp4 * (std::get< 0 >( phi ))[ 0 ] + tmp6 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
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
      case 1:
        {
          return true;
        }
        break;
      case 2:
        {
          return true;
        }
        break;
      case 3:
        {
          return true;
        }
        break;
      case 4:
        {
          return true;
        }
        break;
      default:
        {
          return false;
        }
      }
    }

    template< class Point >
    void dirichlet ( int bndId, const Point &x, RRangeType &result ) const
    {
      switch( bndId )
      {
      case 1:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp2 = -1 * tmp1;
          const auto tmp3 = tmp0[ 1 ] + tmp2;
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp5 = -1 * tmp4;
          const auto tmp6 = tmp0[ 0 ] + tmp5;
          const auto tmp7 = tmp6 * tmp3;
          result[ 0 ] = tmp7;
        }
        break;
      case 2:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp2 = -1 * tmp1;
          const auto tmp3 = tmp0[ 1 ] + tmp2;
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp5 = -1 * tmp4;
          const auto tmp6 = tmp0[ 0 ] + tmp5;
          const auto tmp7 = tmp6 * tmp3;
          result[ 0 ] = tmp7;
        }
        break;
      case 3:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp2 = -1 * tmp1;
          const auto tmp3 = tmp0[ 1 ] + tmp2;
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp5 = -1 * tmp4;
          const auto tmp6 = tmp0[ 0 ] + tmp5;
          const auto tmp7 = tmp6 * tmp3;
          result[ 0 ] = tmp7;
        }
        break;
      case 4:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp2 = -1 * tmp1;
          const auto tmp3 = tmp0[ 1 ] + tmp2;
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp5 = -1 * tmp4;
          const auto tmp6 = tmp0[ 0 ] + tmp5;
          const auto tmp7 = tmp6 * tmp3;
          result[ 0 ] = tmp7;
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
      case 1:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp2 = -1 * tmp1;
          const auto tmp3 = tmp0[ 1 ] + tmp2;
          const auto tmp4 = tmp0[ 0 ] + tmp0[ 0 ];
          const auto tmp5 = -1 * tmp4;
          const auto tmp6 = 1 + tmp5;
          const auto tmp7 = tmp6 * tmp3;
          const auto tmp8 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp9 = -1 * tmp8;
          const auto tmp10 = tmp0[ 0 ] + tmp9;
          const auto tmp11 = tmp0[ 1 ] + tmp0[ 1 ];
          const auto tmp12 = -1 * tmp11;
          const auto tmp13 = 1 + tmp12;
          const auto tmp14 = tmp13 * tmp10;
          (result[ 0 ])[ 0 ] = tmp7;
          (result[ 0 ])[ 1 ] = tmp14;
        }
        break;
      case 2:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp2 = -1 * tmp1;
          const auto tmp3 = tmp0[ 1 ] + tmp2;
          const auto tmp4 = tmp0[ 0 ] + tmp0[ 0 ];
          const auto tmp5 = -1 * tmp4;
          const auto tmp6 = 1 + tmp5;
          const auto tmp7 = tmp6 * tmp3;
          const auto tmp8 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp9 = -1 * tmp8;
          const auto tmp10 = tmp0[ 0 ] + tmp9;
          const auto tmp11 = tmp0[ 1 ] + tmp0[ 1 ];
          const auto tmp12 = -1 * tmp11;
          const auto tmp13 = 1 + tmp12;
          const auto tmp14 = tmp13 * tmp10;
          (result[ 0 ])[ 0 ] = tmp7;
          (result[ 0 ])[ 1 ] = tmp14;
        }
        break;
      case 3:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp2 = -1 * tmp1;
          const auto tmp3 = tmp0[ 1 ] + tmp2;
          const auto tmp4 = tmp0[ 0 ] + tmp0[ 0 ];
          const auto tmp5 = -1 * tmp4;
          const auto tmp6 = 1 + tmp5;
          const auto tmp7 = tmp6 * tmp3;
          const auto tmp8 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp9 = -1 * tmp8;
          const auto tmp10 = tmp0[ 0 ] + tmp9;
          const auto tmp11 = tmp0[ 1 ] + tmp0[ 1 ];
          const auto tmp12 = -1 * tmp11;
          const auto tmp13 = 1 + tmp12;
          const auto tmp14 = tmp13 * tmp10;
          (result[ 0 ])[ 0 ] = tmp7;
          (result[ 0 ])[ 1 ] = tmp14;
        }
        break;
      case 4:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp2 = -1 * tmp1;
          const auto tmp3 = tmp0[ 1 ] + tmp2;
          const auto tmp4 = tmp0[ 0 ] + tmp0[ 0 ];
          const auto tmp5 = -1 * tmp4;
          const auto tmp6 = 1 + tmp5;
          const auto tmp7 = tmp6 * tmp3;
          const auto tmp8 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp9 = -1 * tmp8;
          const auto tmp10 = tmp0[ 0 ] + tmp9;
          const auto tmp11 = tmp0[ 1 ] + tmp0[ 1 ];
          const auto tmp12 = -1 * tmp11;
          const auto tmp13 = 1 + tmp12;
          const auto tmp14 = tmp13 * tmp10;
          (result[ 0 ])[ 0 ] = tmp7;
          (result[ 0 ])[ 1 ] = tmp14;
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

  public:
    template< class Point, class DRangeType >
    RRangeType hessStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linHessStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType gradStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      using std::pow;
      const auto tmp0 = std::pow( u, 2 );
      const auto tmp1 = 1.0 + tmp0;
      result[ 0 ] = tmp1;
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linGradStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result[ 0 ] = 0.0;
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType massStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linMassStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }
  };

} // namespace Integrands_ddf3ff4d527ccaee8da94c708bfe32f2v1_3

PYBIND11_MODULE( vemintegrands_ddf3ff4d527ccaee8da94c708bfe32f2v1_3_f084d41ce9ff164524e5c3ea6f8082cc, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView > GridPart;
  typedef Integrands_ddf3ff4d527ccaee8da94c708bfe32f2v1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_ddf3ff4d527ccaee8da94c708bfe32f2v1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/vemintegrands_ddf3ff4d527ccaee8da94c708bfe32f2v1_3_f084d41ce9ff164524e5c3ea6f8082cc.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_ddf3ff4d527ccaee8da94c708bfe32f2v1_3
