#ifndef GuardIntegrands_b982bf4271de9597e756a64ec0d6c815v1_3
#define GuardIntegrands_b982bf4271de9597e756a64ec0d6c815v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_b982bf4271de9597e756a64ec0d6c815v1_3
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
    typedef std::tuple< Dune::FieldMatrix< double, 1, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldMatrix< double, 1, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      const auto tmp0 = std::get< 0 >( u );
      return RangeValueType{ { { (tmp0[ 0 ])[ 0 ], (tmp0[ 0 ])[ 1 ] } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { { ((std::get< 0 >( phi ))[ 0 ])[ 0 ], ((std::get< 0 >( phi ))[ 0 ])[ 1 ] } } };
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
          using std::atan2;
          using std::pow;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = std::atan2( tmp0[ 1 ], tmp0[ 0 ] );
          const auto tmp2 = (tmp0[ 1 ] < 0.0 ? 6.283185307179586 : 0.0) + tmp1;
          const auto tmp3 = 0.5625 * tmp2;
          const auto tmp4 = std::sin( tmp3 );
          const auto tmp5 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp6 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp7 = tmp6 + tmp5;
          const auto tmp8 = std::pow( tmp7, 0.28125 );
          const auto tmp9 = tmp8 * tmp4;
          result[ 0 ] = tmp9;
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
          using std::atan2;
          using std::cos;
          using std::pow;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = std::atan2( tmp0[ 1 ], tmp0[ 0 ] );
          const auto tmp2 = (tmp0[ 1 ] < 0.0 ? 6.283185307179586 : 0.0) + tmp1;
          const auto tmp3 = 0.5625 * tmp2;
          const auto tmp4 = std::sin( tmp3 );
          const auto tmp5 = tmp0[ 1 ] * tmp0[ 1 ];
          const auto tmp6 = tmp0[ 0 ] * tmp0[ 0 ];
          const auto tmp7 = tmp6 + tmp5;
          const auto tmp8 = std::pow( tmp7, -0.71875 );
          const auto tmp9 = tmp0[ 0 ] + tmp0[ 0 ];
          const auto tmp10 = 0.28125 * tmp9;
          const auto tmp11 = tmp10 * tmp8;
          const auto tmp12 = tmp11 * tmp4;
          const auto tmp13 = std::pow( tmp7, 0.28125 );
          const auto tmp14 = std::cos( tmp3 );
          const auto tmp15 = std::pow( tmp0[ 1 ], 2 );
          const auto tmp16 = std::pow( tmp0[ 0 ], 2 );
          const auto tmp17 = tmp16 + tmp15;
          const auto tmp18 = -1 * tmp0[ 1 ];
          const auto tmp19 = tmp18 / tmp17;
          const auto tmp20 = 0.5625 * tmp19;
          const auto tmp21 = tmp20 * tmp14;
          const auto tmp22 = tmp21 * tmp13;
          const auto tmp23 = tmp22 + tmp12;
          const auto tmp24 = tmp0[ 1 ] + tmp0[ 1 ];
          const auto tmp25 = 0.28125 * tmp24;
          const auto tmp26 = tmp25 * tmp8;
          const auto tmp27 = tmp26 * tmp4;
          const auto tmp28 = tmp0[ 0 ] / tmp17;
          const auto tmp29 = 0.5625 * tmp28;
          const auto tmp30 = tmp29 * tmp14;
          const auto tmp31 = tmp30 * tmp13;
          const auto tmp32 = tmp31 + tmp27;
          (result[ 0 ])[ 0 ] = tmp23;
          (result[ 0 ])[ 1 ] = tmp32;
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

} // namespace Integrands_b982bf4271de9597e756a64ec0d6c815v1_3

PYBIND11_MODULE( integrands_b982bf4271de9597e756a64ec0d6c815v1_3_af122c1df944c95cd395ec0f91d0f970, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > > GridPart;
  typedef Integrands_b982bf4271de9597e756a64ec0d6c815v1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_b982bf4271de9597e756a64ec0d6c815v1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_b982bf4271de9597e756a64ec0d6c815v1_3_af122c1df944c95cd395ec0f91d0f970.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_b982bf4271de9597e756a64ec0d6c815v1_3
