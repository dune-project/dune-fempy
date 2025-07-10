#ifndef GuardIntegrands_a76beb02c9ee134ee7b76e6fe9418447v1_3
#define GuardIntegrands_a76beb02c9ee134ee7b76e6fe9418447v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <cmath>
#include <tuple>
#include <utility>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_a76beb02c9ee134ee7b76e6fe9418447v1_3
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
    typedef double Coneps;
    typedef std::tuple< std::shared_ptr< Coneps > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple<> CoefficientTupleType;
    static constexpr bool gridPartValid = true;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Coneps >( (Coneps(0)) );
    }

    bool init ( const EntityType &entity )
    {
      entity_[ (static_cast< std::size_t >( Side::in )) ] = entity;
      geometry_.emplace( this->entity().geometry() );
      cellVolume_[ (static_cast< std::size_t >( Side::in )) ] = (geometry().volume());
      return true;
    }

    void unbind ()
    {}

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      const auto& intersectionGeom = intersection_.geometry();
      facetArea_ = (intersectionGeom.volume());
      entity_[ (static_cast< std::size_t >( Side::in )) ] = (intersection.inside());
      geometry_.emplace( this->entity().geometry() );
      cellVolume_[ (static_cast< std::size_t >( Side::in )) ] = (geometry().volume());
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
        const Geometry& outsideGeom = entity_[ static_cast< std::size_t >( Side::out ) ].geometry();
        cellVolume_[ static_cast< std::size_t >( Side::out ) ] = outsideGeom.volume();
      }
      return true;
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
      double tmp0 = constant< 0 >();
      const auto tmp1 = std::get< 1 >( u );
      const auto tmp2 = (tmp1[ 0 ])[ 0 ] * tmp0;
      const auto tmp3 = std::get< 0 >( u );
      const auto tmp4 = -1 * tmp3[ 0 ];
      const auto tmp5 = tmp4 + tmp2;
      const auto tmp6 = (tmp1[ 0 ])[ 1 ] * tmp0;
      return RangeValueType{ { 0 }, { { tmp5, tmp6 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 0 >();
      return [ tmp0 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { -1 * (std::get< 0 >( phi ))[ 0 ] + tmp0 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], tmp0 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::atan;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = -1 * tmp0[ 0 ];
      const auto tmp2 = 1 + tmp1;
      const auto tmp3 = 1 + tmp0[ 0 ];
      const auto tmp4 = tmp3 * tmp2;
      const auto tmp5 = 10 * tmp0[ 1 ];
      const auto tmp6 = std::atan( tmp5 );
      GlobalCoordinateType tmp7 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp8 = std::abs( tmp7[ 0 ] );
      const auto tmp9 = tmp7[ 0 ] + tmp8;
      const auto tmp10 = tmp9 / 2.0;
      const auto tmp11 = -1 * tmp10;
      const auto tmp12 = tmp7[ 0 ] + tmp11;
      const auto tmp13 = tmp12 * (tmp0[ 0 ] < 0.0 ? tmp6 : 0.0);
      const auto tmp14 = std::get< 0 >( u );
      const auto tmp15 = tmp14[ 0 ] * tmp10;
      const auto tmp16 = tmp15 + tmp13;
      const auto tmp17 = tmp16 * (tmp4 < 1e-10 ? 1 : 0.0);
      auto tmp18 = facetArea_;
      auto tmp19 = cellVolume();
      const auto tmp20 = tmp19 / tmp18;
      double tmp21 = constant< 0 >();
      const auto tmp22 = 40 * tmp21;
      const auto tmp23 = tmp22 / tmp20;
      const auto tmp24 = -1 * (tmp0[ 0 ] < 0.0 ? tmp6 : 0.0);
      const auto tmp25 = tmp14[ 0 ] + tmp24;
      const auto tmp26 = tmp25 * tmp23;
      const auto tmp27 = tmp26 * (tmp4 < 1e-10 ? 1 : 0.0);
      const auto tmp28 = std::get< 1 >( u );
      const auto tmp29 = (tmp28[ 0 ])[ 1 ] * tmp7[ 1 ];
      const auto tmp30 = tmp7[ 0 ] * (tmp28[ 0 ])[ 0 ];
      const auto tmp31 = tmp30 + tmp29;
      const auto tmp32 = tmp21 * tmp31;
      const auto tmp33 = tmp32 * (tmp4 < 1e-10 ? 1 : 0.0);
      const auto tmp34 = -1 * tmp33;
      const auto tmp35 = tmp34 + tmp27;
      const auto tmp36 = tmp35 + tmp17;
      return RangeValueType{ { tmp36 }, { { 0, 0 } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = -1 * tmp0[ 0 ];
      const auto tmp2 = 1 + tmp1;
      const auto tmp3 = 1 + tmp0[ 0 ];
      const auto tmp4 = tmp3 * tmp2;
      auto tmp5 = facetArea_;
      auto tmp6 = cellVolume();
      const auto tmp7 = tmp6 / tmp5;
      double tmp8 = constant< 0 >();
      const auto tmp9 = 40 * tmp8;
      const auto tmp10 = tmp9 / tmp7;
      const auto tmp11 = tmp10 * (tmp4 < 1e-10 ? 1 : 0.0);
      GlobalCoordinateType tmp12 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp13 = std::abs( tmp12[ 0 ] );
      const auto tmp14 = tmp12[ 0 ] + tmp13;
      const auto tmp15 = tmp14 / 2.0;
      const auto tmp16 = tmp15 * (tmp4 < 1e-10 ? 1 : 0.0);
      const auto tmp17 = tmp16 + tmp11;
      const auto tmp18 = tmp12[ 0 ] * tmp8;
      const auto tmp19 = tmp18 * (tmp4 < 1e-10 ? 1 : 0.0);
      const auto tmp20 = -1 * tmp19;
      const auto tmp21 = tmp12[ 1 ] * tmp8;
      const auto tmp22 = tmp21 * (tmp4 < 1e-10 ? 1 : 0.0);
      const auto tmp23 = -1 * tmp22;
      return [ tmp17, tmp20, tmp23 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (tmp17 * (std::get< 0 >( phi ))[ 0 ] + tmp20 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp23 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { { 0, 0 } } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      auto tmp0 = facetArea_;
      auto tmp1 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
      const auto tmp3 = tmp2 + tmp1;
      const auto tmp4 = 0.5 * tmp3;
      const auto tmp5 = tmp4 / tmp0;
      double tmp6 = constant< 0 >();
      const auto tmp7 = 40 * tmp6;
      const auto tmp8 = tmp7 / tmp5;
      const auto tmp9 = std::get< 0 >( uOut );
      const auto tmp10 = -1 * tmp9[ 0 ];
      const auto tmp11 = std::get< 0 >( uIn );
      const auto tmp12 = tmp11[ 0 ] + tmp10;
      const auto tmp13 = tmp12 * tmp8;
      const auto tmp14 = std::get< 1 >( uOut );
      const auto tmp15 = std::get< 1 >( uIn );
      const auto tmp16 = (tmp15[ 0 ])[ 1 ] + (tmp14[ 0 ])[ 1 ];
      const auto tmp17 = 0.5 * tmp16;
      GlobalCoordinateType tmp18 = intersection_.unitOuterNormal( xIn.localPosition() );
      const auto tmp19 = tmp18[ 1 ] * tmp17;
      const auto tmp20 = (tmp15[ 0 ])[ 0 ] + (tmp14[ 0 ])[ 0 ];
      const auto tmp21 = 0.5 * tmp20;
      const auto tmp22 = tmp18[ 0 ] * tmp21;
      const auto tmp23 = tmp22 + tmp19;
      const auto tmp24 = tmp6 * tmp23;
      const auto tmp25 = -1 * tmp24;
      const auto tmp26 = tmp25 + tmp13;
      const auto tmp27 = std::abs( tmp18[ 0 ] );
      const auto tmp28 = tmp18[ 0 ] + tmp27;
      const auto tmp29 = tmp28 / 2.0;
      const auto tmp30 = tmp11[ 0 ] * tmp29;
      const auto tmp31 = -1 * tmp18[ 0 ];
      const auto tmp32 = std::abs( tmp31 );
      const auto tmp33 = tmp31 + tmp32;
      const auto tmp34 = tmp33 / 2.0;
      const auto tmp35 = tmp9[ 0 ] * tmp34;
      const auto tmp36 = -1 * tmp35;
      const auto tmp37 = tmp36 + tmp30;
      const auto tmp38 = tmp37 + tmp26;
      const auto tmp39 = tmp6 * tmp12;
      const auto tmp40 = 0.5 * tmp18[ 0 ];
      const auto tmp41 = tmp40 * tmp39;
      const auto tmp42 = -1 * tmp41;
      const auto tmp43 = 0.5 * tmp18[ 1 ];
      const auto tmp44 = tmp43 * tmp39;
      const auto tmp45 = -1 * tmp44;
      const auto tmp46 = -1 * tmp37;
      const auto tmp47 = -1 * tmp13;
      const auto tmp48 = -1 * tmp25;
      const auto tmp49 = tmp48 + tmp47;
      const auto tmp50 = tmp49 + tmp46;
      return std::make_pair( RangeValueType{ { tmp38 }, { { tmp42, tmp45 } } }, RangeValueType{ { tmp50 }, { { tmp42, tmp45 } } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      auto tmp0 = facetArea_;
      auto tmp1 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
      const auto tmp3 = tmp2 + tmp1;
      const auto tmp4 = 0.5 * tmp3;
      const auto tmp5 = tmp4 / tmp0;
      double tmp6 = constant< 0 >();
      const auto tmp7 = 40 * tmp6;
      const auto tmp8 = tmp7 / tmp5;
      GlobalCoordinateType tmp9 = intersection_.unitOuterNormal( xIn.localPosition() );
      const auto tmp10 = std::abs( tmp9[ 0 ] );
      const auto tmp11 = tmp9[ 0 ] + tmp10;
      const auto tmp12 = tmp11 / 2.0;
      const auto tmp13 = tmp12 + tmp8;
      const auto tmp14 = 0.5 * tmp9[ 0 ];
      const auto tmp15 = tmp6 * tmp14;
      const auto tmp16 = -1 * tmp15;
      const auto tmp17 = 0.5 * tmp9[ 1 ];
      const auto tmp18 = tmp6 * tmp17;
      const auto tmp19 = -1 * tmp18;
      const auto tmp20 = -1 * tmp8;
      const auto tmp21 = -1 * tmp12;
      const auto tmp22 = tmp21 + tmp20;
      const auto tmp23 = -1 * tmp16;
      const auto tmp24 = -1 * tmp19;
      const auto tmp25 = -1 * tmp9[ 0 ];
      const auto tmp26 = std::abs( tmp25 );
      const auto tmp27 = tmp25 + tmp26;
      const auto tmp28 = tmp27 / 2.0;
      const auto tmp29 = -1 * tmp28;
      const auto tmp30 = tmp29 + tmp20;
      const auto tmp31 = -1 * tmp6;
      const auto tmp32 = tmp14 * tmp31;
      const auto tmp33 = -1 * tmp32;
      const auto tmp34 = tmp17 * tmp31;
      const auto tmp35 = -1 * tmp34;
      const auto tmp36 = -1 * tmp20;
      const auto tmp37 = -1 * tmp29;
      const auto tmp38 = tmp37 + tmp36;
      return std::make_pair( [ tmp13, tmp16, tmp19, tmp22, tmp23, tmp24 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { (tmp13 * (std::get< 0 >( phiIn ))[ 0 ] + tmp16 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp19 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ] }, { { tmp16 * (std::get< 0 >( phiIn ))[ 0 ], tmp19 * (std::get< 0 >( phiIn ))[ 0 ] } } }, RangeValueType{ { (tmp22 * (std::get< 0 >( phiIn ))[ 0 ] + tmp23 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp24 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ] }, { { tmp16 * (std::get< 0 >( phiIn ))[ 0 ], tmp19 * (std::get< 0 >( phiIn ))[ 0 ] } } } );
        }, [ tmp16, tmp19, tmp23, tmp24, tmp30, tmp33, tmp35, tmp38 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { (tmp30 * (std::get< 0 >( phiOut ))[ 0 ] + tmp16 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp19 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ] }, { { tmp33 * (std::get< 0 >( phiOut ))[ 0 ], tmp35 * (std::get< 0 >( phiOut ))[ 0 ] } } }, RangeValueType{ { (tmp38 * (std::get< 0 >( phiOut ))[ 0 ] + tmp23 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp24 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ] }, { { tmp33 * (std::get< 0 >( phiOut ))[ 0 ], tmp35 * (std::get< 0 >( phiOut ))[ 0 ] } } } );
        } );
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

    const Coneps &coneps () const
    {
      return *std::get< 0 >( constants_ );
    }

    Coneps &coneps ()
    {
      return *std::get< 0 >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_[ (static_cast< std::size_t >( Side::in )) ];
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }

    const ctype cellVolume () const
    {
      return cellVolume_[ (static_cast< std::size_t >( Side::in )) ];
    }
    std::array< EntityType, 2 > entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    std::array< ctype, 2 > cellVolume_;
    ctype facetArea_;
    ConstantTupleType constants_;
    std::array< std::tuple<  >, 2 > coefficients_;
  };

} // namespace Integrands_a76beb02c9ee134ee7b76e6fe9418447v1_3

PYBIND11_MODULE( integrands_a76beb02c9ee134ee7b76e6fe9418447v1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_a76beb02c9ee134ee7b76e6fe9418447v1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_a76beb02c9ee134ee7b76e6fe9418447v1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_a76beb02c9ee134ee7b76e6fe9418447v1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property( "eps", [] ( Integrands &self ) -> Integrands::Coneps { return self.coneps(); }, [] ( Integrands &self, const Integrands::Coneps &v ) { self.coneps() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_a76beb02c9ee134ee7b76e6fe9418447v1_3
