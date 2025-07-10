#ifndef GuardIntegrands_6442882efe640a3313d957a8885ec062v1_3
#define GuardIntegrands_6442882efe640a3313d957a8885ec062v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <cmath>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_6442882efe640a3313d957a8885ec062v1_3
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
      const auto tmp0 = std::get< 0 >( u );
      const auto tmp1 = 1 + tmp0[ 0 ];
      const auto tmp2 = std::pow( tmp1, 2 );
      const auto tmp3 = tmp0[ 0 ] * tmp2;
      GlobalCoordinateType tmp4 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp5 = std::pow( tmp4[ 0 ], 3 );
      const auto tmp6 = std::pow( tmp4[ 1 ], 3 );
      const auto tmp7 = -1 * tmp6;
      const auto tmp8 = tmp7 + tmp5;
      const auto tmp9 = 0.3333333333333333 * tmp8;
      const auto tmp10 = -1 * tmp9;
      const auto tmp11 = std::pow( tmp4[ 1 ], 2 );
      const auto tmp12 = std::pow( tmp4[ 0 ], 2 );
      const auto tmp13 = tmp12 + tmp11;
      const auto tmp14 = 0.5 * tmp13;
      const auto tmp15 = tmp14 + tmp10;
      const auto tmp16 = 1 + tmp15;
      const auto tmp17 = 1 + tmp16;
      const auto tmp18 = std::pow( tmp17, 2 );
      const auto tmp19 = tmp16 * tmp18;
      const auto tmp20 = 2 * tmp4[ 1 ];
      const auto tmp21 = 3 * tmp20;
      const auto tmp22 = -1 * tmp21;
      const auto tmp23 = 0.3333333333333333 * tmp22;
      const auto tmp24 = -1 * tmp23;
      const auto tmp25 = 1.0 + tmp24;
      const auto tmp26 = 2 * tmp4[ 0 ];
      const auto tmp27 = 3 * tmp26;
      const auto tmp28 = 0.3333333333333333 * tmp27;
      const auto tmp29 = -1 * tmp28;
      const auto tmp30 = 1.0 + tmp29;
      const auto tmp31 = tmp30 + tmp25;
      const auto tmp32 = -1 * tmp31;
      const auto tmp33 = tmp32 + tmp19;
      const auto tmp34 = -1 * tmp33;
      const auto tmp35 = tmp34 + tmp3;
      const auto tmp36 = std::get< 1 >( u );
      return RangeValueType{ { tmp35 }, { { (tmp36[ 0 ])[ 0 ], (tmp36[ 0 ])[ 1 ] } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      const auto tmp0 = std::get< 0 >( u );
      const auto tmp1 = 1 + tmp0[ 0 ];
      const auto tmp2 = std::pow( tmp1, 2 );
      const auto tmp3 = 2 * tmp1;
      const auto tmp4 = tmp0[ 0 ] * tmp3;
      const auto tmp5 = tmp4 + tmp2;
      return [ tmp5 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp5 * (std::get< 0 >( phi ))[ 0 ] }, { { ((std::get< 1 >( phi ))[ 0 ])[ 0 ], ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = std::pow( tmp0[ 0 ], 2 );
      const auto tmp2 = 3 * tmp1;
      const auto tmp3 = 0.3333333333333333 * tmp2;
      const auto tmp4 = -1 * tmp3;
      const auto tmp5 = 2 * tmp0[ 0 ];
      const auto tmp6 = 0.5 * tmp5;
      const auto tmp7 = tmp6 + tmp4;
      GlobalCoordinateType tmp8 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp9 = tmp8[ 0 ] * tmp7;
      const auto tmp10 = std::pow( tmp0[ 1 ], 2 );
      const auto tmp11 = 3 * tmp10;
      const auto tmp12 = -1 * tmp11;
      const auto tmp13 = 0.3333333333333333 * tmp12;
      const auto tmp14 = -1 * tmp13;
      const auto tmp15 = 2 * tmp0[ 1 ];
      const auto tmp16 = 0.5 * tmp15;
      const auto tmp17 = tmp16 + tmp14;
      const auto tmp18 = tmp8[ 1 ] * tmp17;
      const auto tmp19 = tmp18 + tmp9;
      const auto tmp20 = -1 * tmp19;
      return RangeValueType{ { tmp20 }, { { 0, 0 } } };
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
      return false;
    }

    bool isDirichletIntersection ( const IntersectionType &intersection, DirichletComponentType &dirichletComponent ) const
    {
      return false;
    }

    template< class Point >
    void dirichlet ( int bndId, const Point &x, RRangeType &result ) const
    {}

    template< class Point >
    void dDirichlet ( int bndId, const Point &x, RJacobianRangeType &result ) const
    {}

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

} // namespace Integrands_6442882efe640a3313d957a8885ec062v1_3

PYBIND11_MODULE( integrands_6442882efe640a3313d957a8885ec062v1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_6442882efe640a3313d957a8885ec062v1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_6442882efe640a3313d957a8885ec062v1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_6442882efe640a3313d957a8885ec062v1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_6442882efe640a3313d957a8885ec062v1_3
