#ifndef GuardIntegrands_f4d7a331d8f4abd1acfa3dd142f69f59v1_4
#define GuardIntegrands_f4d7a331d8f4abd1acfa3dd142f69f59v1_4
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
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_f4d7a331d8f4abd1acfa3dd142f69f59v1_4
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

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< Coneps >( (constant<0>()) );
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
        assert(tmp0 == tmp0);
      const auto tmp1 = std::get< 1 >( u );
        assert(tmp1 == tmp1);
      const auto tmp2 = (tmp1[ 0 ])[ 0 ] * tmp0;
        assert(tmp2 == tmp2);
      const auto tmp3 = std::get< 0 >( u );
        assert(tmp3 == tmp3);
      const auto tmp4 = -1 * tmp3[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 + tmp2;
        assert(tmp5 == tmp5);
      const auto tmp6 = (tmp1[ 0 ])[ 1 ] * tmp0;
        assert(tmp6 == tmp6);
      return RangeValueType{ { 0 }, { { tmp5, tmp6 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
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
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = 1 + tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = 1 + tmp0[ 0 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 * tmp2;
        assert(tmp4 == tmp4);
      const auto tmp5 = 10 * tmp0[ 1 ];
        assert(tmp5 == tmp5);
      const auto tmp6 = std::atan( tmp5 );
        assert(tmp6 == tmp6);
      GlobalCoordinateType tmp7 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp7 == tmp7);
      const auto tmp8 = std::abs( tmp7[ 0 ] );
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp7[ 0 ] + tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 / 2.0;
        assert(tmp10 == tmp10);
      const auto tmp11 = -1 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp7[ 0 ] + tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 * (tmp0[ 0 ] < 0.0 ? tmp6 : 0.0);
        assert(tmp13 == tmp13);
      const auto tmp14 = std::get< 0 >( u );
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14[ 0 ] * tmp10;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 + tmp13;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 * (tmp4 < 1e-10 ? 1 : 0.0);
        assert(tmp17 == tmp17);
      auto tmp18 = facetArea_;
        assert(tmp18 == tmp18);
      auto tmp19 = cellVolume();
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 / tmp18;
        assert(tmp20 == tmp20);
      double tmp21 = constant< 0 >();
        assert(tmp21 == tmp21);
      const auto tmp22 = 40 * tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 / tmp20;
        assert(tmp23 == tmp23);
      const auto tmp24 = -1 * (tmp0[ 0 ] < 0.0 ? tmp6 : 0.0);
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp14[ 0 ] + tmp24;
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp25 * tmp23;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 * (tmp4 < 1e-10 ? 1 : 0.0);
        assert(tmp27 == tmp27);
      const auto tmp28 = std::get< 1 >( u );
        assert(tmp28 == tmp28);
      const auto tmp29 = (tmp28[ 0 ])[ 1 ] * tmp7[ 1 ];
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp7[ 0 ] * (tmp28[ 0 ])[ 0 ];
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp30 + tmp29;
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 * tmp21;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 * (tmp4 < 1e-10 ? 1 : 0.0);
        assert(tmp33 == tmp33);
      const auto tmp34 = -1 * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp34 + tmp27;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 + tmp17;
        assert(tmp36 == tmp36);
      return RangeValueType{ { tmp36 }, { { 0, 0 } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = 1 + tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = 1 + tmp0[ 0 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 * tmp2;
        assert(tmp4 == tmp4);
      auto tmp5 = facetArea_;
        assert(tmp5 == tmp5);
      auto tmp6 = cellVolume();
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6 / tmp5;
        assert(tmp7 == tmp7);
      double tmp8 = constant< 0 >();
        assert(tmp8 == tmp8);
      const auto tmp9 = 40 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 / tmp7;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * (tmp4 < 1e-10 ? 1 : 0.0);
        assert(tmp11 == tmp11);
      GlobalCoordinateType tmp12 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp12 == tmp12);
      const auto tmp13 = std::abs( tmp12[ 0 ] );
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp12[ 0 ] + tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 / 2.0;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 * (tmp4 < 1e-10 ? 1 : 0.0);
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 + tmp11;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp12[ 0 ] * tmp8;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp18 * (tmp4 < 1e-10 ? 1 : 0.0);
        assert(tmp19 == tmp19);
      const auto tmp20 = -1 * tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp12[ 1 ] * tmp8;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 * (tmp4 < 1e-10 ? 1 : 0.0);
        assert(tmp22 == tmp22);
      const auto tmp23 = -1 * tmp22;
        assert(tmp23 == tmp23);
      return [ tmp17, tmp20, tmp23 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (tmp17 * (std::get< 0 >( phi ))[ 0 ] + tmp20 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp23 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { { 0, 0 } } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      auto tmp0 = facetArea_;
        assert(tmp0 == tmp0);
      auto tmp1 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
        assert(tmp1 == tmp1);
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = 0.5 * tmp3;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 / tmp0;
        assert(tmp5 == tmp5);
      double tmp6 = constant< 0 >();
        assert(tmp6 == tmp6);
      const auto tmp7 = 40 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 / tmp5;
        assert(tmp8 == tmp8);
      const auto tmp9 = std::get< 0 >( uOut );
        assert(tmp9 == tmp9);
      const auto tmp10 = -1 * tmp9[ 0 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = std::get< 0 >( uIn );
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11[ 0 ] + tmp10;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 * tmp8;
        assert(tmp13 == tmp13);
      GlobalCoordinateType tmp14 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp14 == tmp14);
      const auto tmp15 = -1 * tmp14[ 1 ];
        assert(tmp15 == tmp15);
      const auto tmp16 = std::get< 1 >( uOut );
        assert(tmp16 == tmp16);
      const auto tmp17 = (tmp16[ 0 ])[ 1 ] * tmp15;
        assert(tmp17 == tmp17);
      const auto tmp18 = -1 * tmp14[ 0 ];
        assert(tmp18 == tmp18);
      const auto tmp19 = (tmp16[ 0 ])[ 0 ] * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 + tmp17;
        assert(tmp20 == tmp20);
      const auto tmp21 = std::get< 1 >( uIn );
        assert(tmp21 == tmp21);
      const auto tmp22 = (tmp21[ 0 ])[ 1 ] * tmp14[ 1 ];
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp14[ 0 ] * (tmp21[ 0 ])[ 0 ];
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 + tmp22;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 + tmp20;
        assert(tmp25 == tmp25);
      const auto tmp26 = 0.5 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 * tmp6;
        assert(tmp27 == tmp27);
      const auto tmp28 = -1 * tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 + tmp13;
        assert(tmp29 == tmp29);
      const auto tmp30 = std::abs( tmp14[ 0 ] );
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp14[ 0 ] + tmp30;
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 / 2.0;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp11[ 0 ] * tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = std::abs( tmp18 );
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp18 + tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 / 2.0;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp9[ 0 ] * tmp36;
        assert(tmp37 == tmp37);
      const auto tmp38 = -1 * tmp37;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 + tmp33;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 + tmp29;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp12 * tmp6;
        assert(tmp41 == tmp41);
      const auto tmp42 = 0.5 * tmp14[ 0 ];
        assert(tmp42 == tmp42);
      const auto tmp43 = tmp42 * tmp41;
        assert(tmp43 == tmp43);
      const auto tmp44 = -1 * tmp43;
        assert(tmp44 == tmp44);
      const auto tmp45 = 0.5 * tmp14[ 1 ];
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 * tmp41;
        assert(tmp46 == tmp46);
      const auto tmp47 = -1 * tmp46;
        assert(tmp47 == tmp47);
      const auto tmp48 = -1 * tmp39;
        assert(tmp48 == tmp48);
      const auto tmp49 = -1 * tmp13;
        assert(tmp49 == tmp49);
      const auto tmp50 = -1 * tmp28;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp50 + tmp49;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 + tmp48;
        assert(tmp52 == tmp52);
      const auto tmp53 = 0.5 * tmp18;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 * tmp41;
        assert(tmp54 == tmp54);
      const auto tmp55 = -1 * tmp54;
        assert(tmp55 == tmp55);
      const auto tmp56 = 0.5 * tmp15;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp56 * tmp41;
        assert(tmp57 == tmp57);
      const auto tmp58 = -1 * tmp57;
        assert(tmp58 == tmp58);
      return std::make_pair( RangeValueType{ { tmp40 }, { { tmp44, tmp47 } } }, RangeValueType{ { tmp52 }, { { tmp55, tmp58 } } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      auto tmp0 = facetArea_;
        assert(tmp0 == tmp0);
      auto tmp1 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
        assert(tmp1 == tmp1);
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = 0.5 * tmp3;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 / tmp0;
        assert(tmp5 == tmp5);
      double tmp6 = constant< 0 >();
        assert(tmp6 == tmp6);
      const auto tmp7 = 40 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 / tmp5;
        assert(tmp8 == tmp8);
      GlobalCoordinateType tmp9 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp9 == tmp9);
      const auto tmp10 = std::abs( tmp9[ 0 ] );
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp9[ 0 ] + tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 / 2.0;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 + tmp8;
        assert(tmp13 == tmp13);
      const auto tmp14 = 0.5 * tmp9[ 0 ];
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp6;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = 0.5 * tmp9[ 1 ];
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp6;
        assert(tmp18 == tmp18);
      const auto tmp19 = -1 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = -1 * tmp8;
        assert(tmp20 == tmp20);
      const auto tmp21 = -1 * tmp12;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 + tmp20;
        assert(tmp22 == tmp22);
      const auto tmp23 = -1 * tmp16;
        assert(tmp23 == tmp23);
      const auto tmp24 = -1 * tmp19;
        assert(tmp24 == tmp24);
      const auto tmp25 = -1 * tmp9[ 0 ];
        assert(tmp25 == tmp25);
      const auto tmp26 = 0.5 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 * tmp6;
        assert(tmp27 == tmp27);
      const auto tmp28 = -1 * tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = -1 * tmp9[ 1 ];
        assert(tmp29 == tmp29);
      const auto tmp30 = 0.5 * tmp29;
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp30 * tmp6;
        assert(tmp31 == tmp31);
      const auto tmp32 = -1 * tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = std::abs( tmp25 );
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp25 + tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp34 / 2.0;
        assert(tmp35 == tmp35);
      const auto tmp36 = -1 * tmp35;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 + tmp20;
        assert(tmp37 == tmp37);
      const auto tmp38 = -1 * tmp6;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp14 * tmp38;
        assert(tmp39 == tmp39);
      const auto tmp40 = -1 * tmp39;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp17 * tmp38;
        assert(tmp41 == tmp41);
      const auto tmp42 = -1 * tmp41;
        assert(tmp42 == tmp42);
      const auto tmp43 = -1 * tmp20;
        assert(tmp43 == tmp43);
      const auto tmp44 = -1 * tmp36;
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp44 + tmp43;
        assert(tmp45 == tmp45);
      const auto tmp46 = -1 * tmp28;
        assert(tmp46 == tmp46);
      const auto tmp47 = -1 * tmp32;
        assert(tmp47 == tmp47);
      const auto tmp48 = tmp26 * tmp38;
        assert(tmp48 == tmp48);
      const auto tmp49 = -1 * tmp48;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp30 * tmp38;
        assert(tmp50 == tmp50);
      const auto tmp51 = -1 * tmp50;
        assert(tmp51 == tmp51);
      return std::make_pair( [ tmp13, tmp16, tmp19, tmp22, tmp23, tmp24, tmp28, tmp32 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { (tmp13 * (std::get< 0 >( phiIn ))[ 0 ] + tmp16 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp19 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ] }, { { tmp16 * (std::get< 0 >( phiIn ))[ 0 ], tmp19 * (std::get< 0 >( phiIn ))[ 0 ] } } }, RangeValueType{ { (tmp22 * (std::get< 0 >( phiIn ))[ 0 ] + tmp23 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp24 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ] }, { { tmp28 * (std::get< 0 >( phiIn ))[ 0 ], tmp32 * (std::get< 0 >( phiIn ))[ 0 ] } } } );
        }, [ tmp28, tmp32, tmp37, tmp40, tmp42, tmp45, tmp46, tmp47, tmp49, tmp51 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { (tmp37 * (std::get< 0 >( phiOut ))[ 0 ] + tmp28 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp32 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ] }, { { tmp40 * (std::get< 0 >( phiOut ))[ 0 ], tmp42 * (std::get< 0 >( phiOut ))[ 0 ] } } }, RangeValueType{ { (tmp45 * (std::get< 0 >( phiOut ))[ 0 ] + tmp46 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp47 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ] }, { { tmp49 * (std::get< 0 >( phiOut ))[ 0 ], tmp51 * (std::get< 0 >( phiOut ))[ 0 ] } } } );
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

} // namespace Integrands_f4d7a331d8f4abd1acfa3dd142f69f59v1_4

PYBIND11_MODULE( integrands_f4d7a331d8f4abd1acfa3dd142f69f59v1_4_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_f4d7a331d8f4abd1acfa3dd142f69f59v1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_f4d7a331d8f4abd1acfa3dd142f69f59v1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_f4d7a331d8f4abd1acfa3dd142f69f59v1_4_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property( "eps", [] ( Integrands &self ) -> Integrands::Coneps { return self.coneps(); }, [] ( Integrands &self, const Integrands::Coneps &v ) { self.coneps() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_f4d7a331d8f4abd1acfa3dd142f69f59v1_4
