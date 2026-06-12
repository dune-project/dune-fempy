#ifndef GuardIntegrands_1896a171b906827df87904471d1a8fa6v1_4
#define GuardIntegrands_1896a171b906827df87904471d1a8fa6v1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/alugrid/hierarchical.hh>
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

namespace Integrands_1896a171b906827df87904471d1a8fa6v1_4
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
      entity_[ (static_cast< std::size_t >( Side::in )) ] = entity;
      geometry_.emplace( this->entity().geometry() );
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
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
      }
      return true;
    }
    typedef std::tuple< Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::cos;
      using std::pow;
      using std::sin;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::sin( tmp1 );
        assert(tmp2 == tmp2);
      const auto tmp3 = std::pow( tmp2, 2 );
        assert(tmp3 == tmp3);
      const auto tmp4 = 3.141592653589793 * tmp0[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = std::sin( tmp4 );
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 3.141592653589793 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = 3.141592653589793 * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = 3.141592653589793 * tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = 3.141592653589793 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = 2 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 * tmp5;
        assert(tmp13 == tmp13);
      const auto tmp14 = std::cos( tmp4 );
        assert(tmp14 == tmp14);
      const auto tmp15 = 3.141592653589793 * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = 3.141592653589793 * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = 3.141592653589793 * tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = 2 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 * tmp15;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 + tmp13;
        assert(tmp21 == tmp21);
      const auto tmp22 = 2 * tmp8;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 * tmp8;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 + tmp20;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 + tmp21;
        assert(tmp25 == tmp25);
      const auto tmp26 = 2 * tmp15;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp18 * tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp23 + tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 + tmp24;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 + tmp25;
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp30 * tmp3;
        assert(tmp31 == tmp31);
      const auto tmp32 = -1 * tmp2;
        assert(tmp32 == tmp32);
      const auto tmp33 = 3.141592653589793 * tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = 3.141592653589793 * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = 2 * tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 * tmp2;
        assert(tmp36 == tmp36);
      const auto tmp37 = std::cos( tmp1 );
        assert(tmp37 == tmp37);
      const auto tmp38 = 3.141592653589793 * tmp37;
        assert(tmp38 == tmp38);
      const auto tmp39 = 2 * tmp38;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 * tmp38;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 + tmp36;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp22 * tmp5;
        assert(tmp42 == tmp42);
      const auto tmp43 = tmp26 * tmp15;
        assert(tmp43 == tmp43);
      const auto tmp44 = tmp43 + tmp42;
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp44 * tmp41;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 + tmp31;
        assert(tmp46 == tmp46);
      const auto tmp47 = std::pow( tmp5, 2 );
        assert(tmp47 == tmp47);
      const auto tmp48 = -1 * tmp34;
        assert(tmp48 == tmp48);
      const auto tmp49 = 3.141592653589793 * tmp48;
        assert(tmp49 == tmp49);
      const auto tmp50 = 3.141592653589793 * tmp49;
        assert(tmp50 == tmp50);
      const auto tmp51 = 2 * tmp50;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 * tmp2;
        assert(tmp52 == tmp52);
      const auto tmp53 = -1 * tmp38;
        assert(tmp53 == tmp53);
      const auto tmp54 = 3.141592653589793 * tmp53;
        assert(tmp54 == tmp54);
      const auto tmp55 = 3.141592653589793 * tmp54;
        assert(tmp55 == tmp55);
      const auto tmp56 = 2 * tmp55;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp56 * tmp38;
        assert(tmp57 == tmp57);
      const auto tmp58 = tmp57 + tmp52;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp35 * tmp34;
        assert(tmp59 == tmp59);
      const auto tmp60 = tmp59 + tmp57;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp60 + tmp58;
        assert(tmp61 == tmp61);
      const auto tmp62 = tmp55 * tmp39;
        assert(tmp62 == tmp62);
      const auto tmp63 = tmp59 + tmp62;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp63 + tmp60;
        assert(tmp64 == tmp64);
      const auto tmp65 = tmp64 + tmp61;
        assert(tmp65 == tmp65);
      const auto tmp66 = tmp65 * tmp47;
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp45 + tmp66;
        assert(tmp67 == tmp67);
      const auto tmp68 = tmp67 + tmp46;
        assert(tmp68 == tmp68);
      const auto tmp69 = -1 * tmp68;
        assert(tmp69 == tmp69);
      const auto tmp70 = std::get< 1 >( u );
        assert(tmp70 == tmp70);
      const auto tmp71 = ((tmp70[ 0 ])[ 0 ])[ 0 ] + ((tmp70[ 0 ])[ 1 ])[ 1 ];
        assert(tmp71 == tmp71);
      return RangeValueType{ { tmp69 }, { { 0, 0 } }, { { { tmp71, 0.0 }, { 0.0, tmp71 } } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { 0, 0 } }, { { { (((std::get< 1 >( phi ))[ 0 ])[ 0 ])[ 0 ] + (((std::get< 1 >( phi ))[ 0 ])[ 1 ])[ 1 ], 0 }, { 0, (((std::get< 1 >( phi ))[ 0 ])[ 0 ])[ 0 ] + (((std::get< 1 >( phi ))[ 0 ])[ 1 ])[ 1 ] } } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::cos;
      using std::pow;
      using std::sin;
      auto tmp0 = facetArea_;
        assert(tmp0 == tmp0);
      const auto tmp1 = 40 / tmp0;
        assert(tmp1 == tmp1);
      GlobalCoordinateType tmp2 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp2 == tmp2);
      const auto tmp3 = 3.141592653589793 * tmp2[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = std::sin( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = std::pow( tmp4, 2 );
        assert(tmp5 == tmp5);
      const auto tmp6 = 3.141592653589793 * tmp2[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = std::sin( tmp6 );
        assert(tmp7 == tmp7);
      const auto tmp8 = std::cos( tmp6 );
        assert(tmp8 == tmp8);
      const auto tmp9 = 3.141592653589793 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = 2 * tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * tmp7;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 * tmp5;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = std::get< 0 >( u );
        assert(tmp14 == tmp14);
      const auto tmp15 = (tmp14[ 0 ])[ 0 ] + tmp13;
        assert(tmp15 == tmp15);
      GlobalCoordinateType tmp16 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16[ 0 ] * tmp15;
        assert(tmp17 == tmp17);
      const auto tmp18 = std::pow( tmp7, 2 );
        assert(tmp18 == tmp18);
      const auto tmp19 = std::cos( tmp3 );
        assert(tmp19 == tmp19);
      const auto tmp20 = 3.141592653589793 * tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = 2 * tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 * tmp4;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 * tmp18;
        assert(tmp23 == tmp23);
      const auto tmp24 = -1 * tmp23;
        assert(tmp24 == tmp24);
      const auto tmp25 = (tmp14[ 0 ])[ 1 ] + tmp24;
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp16[ 1 ] * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 + tmp17;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 * tmp1;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp16[ 0 ] * tmp28;
        assert(tmp29 == tmp29);
      const auto tmp30 = std::get< 1 >( u );
        assert(tmp30 == tmp30);
      const auto tmp31 = ((tmp30[ 0 ])[ 0 ])[ 0 ] + ((tmp30[ 0 ])[ 1 ])[ 1 ];
        assert(tmp31 == tmp31);
      const auto tmp32 = -1 * tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp16[ 0 ] * tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp33 + tmp29;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp16[ 1 ] * tmp28;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp16[ 1 ] * tmp32;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 + tmp35;
        assert(tmp37 == tmp37);
      const auto tmp38 = (tmp14[ 0 ])[ 1 ] * tmp16[ 1 ];
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp16[ 0 ] * (tmp14[ 0 ])[ 0 ];
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 + tmp38;
        assert(tmp40 == tmp40);
      const auto tmp41 = -1 * tmp40;
        assert(tmp41 == tmp41);
      return RangeValueType{ { 0 }, { { tmp34, tmp37 } }, { { { tmp41, 0.0 }, { 0.0, tmp41 } } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      auto tmp0 = facetArea_;
        assert(tmp0 == tmp0);
      const auto tmp1 = 40 / tmp0;
        assert(tmp1 == tmp1);
      GlobalCoordinateType tmp2 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2[ 0 ] * tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp2[ 0 ] * tmp3;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp2[ 1 ] * tmp1;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp2[ 0 ] * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp2[ 1 ] * tmp3;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp2[ 1 ] * tmp5;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp2[ 0 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = -1 * tmp2[ 1 ];
        assert(tmp10 == tmp10);
      return [ tmp10, tmp4, tmp6, tmp7, tmp8, tmp9 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { ((tmp4 * ((std::get< 0 >( phi ))[ 0 ])[ 0 ] + tmp6 * ((std::get< 0 >( phi ))[ 0 ])[ 1 ]) + tmp9 * (((std::get< 1 >( phi ))[ 0 ])[ 0 ])[ 0 ]) + tmp9 * (((std::get< 1 >( phi ))[ 0 ])[ 1 ])[ 1 ], ((tmp7 * ((std::get< 0 >( phi ))[ 0 ])[ 0 ] + tmp8 * ((std::get< 0 >( phi ))[ 0 ])[ 1 ]) + tmp10 * (((std::get< 1 >( phi ))[ 0 ])[ 0 ])[ 0 ]) + tmp10 * (((std::get< 1 >( phi ))[ 0 ])[ 1 ])[ 1 ] } }, { { { tmp9 * ((std::get< 0 >( phi ))[ 0 ])[ 0 ] + tmp10 * ((std::get< 0 >( phi ))[ 0 ])[ 1 ], 0 }, { 0, tmp9 * ((std::get< 0 >( phi ))[ 0 ])[ 0 ] + tmp10 * ((std::get< 0 >( phi ))[ 0 ])[ 1 ] } } } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      auto tmp0 = facetArea_;
        assert(tmp0 == tmp0);
      const auto tmp1 = 40 / tmp0;
        assert(tmp1 == tmp1);
      GlobalCoordinateType tmp2 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp2 == tmp2);
      const auto tmp3 = -1 * tmp2[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = std::get< 0 >( uOut );
        assert(tmp4 == tmp4);
      const auto tmp5 = (tmp4[ 0 ])[ 1 ] * tmp3;
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp2[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = (tmp4[ 0 ])[ 0 ] * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 + tmp5;
        assert(tmp8 == tmp8);
      const auto tmp9 = std::get< 0 >( uIn );
        assert(tmp9 == tmp9);
      const auto tmp10 = (tmp9[ 0 ])[ 1 ] * tmp2[ 1 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp2[ 0 ] * (tmp9[ 0 ])[ 0 ];
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 + tmp10;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 + tmp8;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp2[ 0 ] * tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp1;
        assert(tmp15 == tmp15);
      const auto tmp16 = std::get< 1 >( uOut );
        assert(tmp16 == tmp16);
      const auto tmp17 = ((tmp16[ 0 ])[ 0 ])[ 0 ] + ((tmp16[ 0 ])[ 1 ])[ 1 ];
        assert(tmp17 == tmp17);
      const auto tmp18 = std::get< 1 >( uIn );
        assert(tmp18 == tmp18);
      const auto tmp19 = ((tmp18[ 0 ])[ 0 ])[ 0 ] + ((tmp18[ 0 ])[ 1 ])[ 1 ];
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 + tmp17;
        assert(tmp20 == tmp20);
      const auto tmp21 = 0.5 * tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = -1 * tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp2[ 0 ] * tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 + tmp15;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp2[ 1 ] * tmp13;
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp25 * tmp1;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp2[ 1 ] * tmp22;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 + tmp26;
        assert(tmp28 == tmp28);
      const auto tmp29 = -0.5 * tmp13;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp13 * tmp6;
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp30 * tmp1;
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp6 * tmp22;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 + tmp31;
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp13 * tmp3;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp34 * tmp1;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp3 * tmp22;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 + tmp35;
        assert(tmp37 == tmp37);
      return std::make_pair( RangeValueType{ { 0 }, { { tmp24, tmp28 } }, { { { tmp29, 0.0 }, { 0.0, tmp29 } } } }, RangeValueType{ { 0 }, { { tmp33, tmp37 } }, { { { tmp29, 0.0 }, { 0.0, tmp29 } } } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      auto tmp0 = facetArea_;
        assert(tmp0 == tmp0);
      const auto tmp1 = 40 / tmp0;
        assert(tmp1 == tmp1);
      GlobalCoordinateType tmp2 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2[ 0 ] * tmp2[ 0 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 * tmp1;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp2[ 0 ] * tmp2[ 1 ];
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 * tmp1;
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp2[ 1 ] * tmp2[ 1 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 * tmp1;
        assert(tmp8 == tmp8);
      const auto tmp9 = -0.5 * tmp2[ 0 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = -0.5 * tmp2[ 1 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = -1 * tmp2[ 0 ];
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp2[ 0 ] * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 * tmp1;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp2[ 1 ] * tmp11;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp1;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp2[ 1 ];
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp2[ 0 ] * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp1;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp2[ 1 ] * tmp16;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 * tmp1;
        assert(tmp20 == tmp20);
      const auto tmp21 = -0.5 * tmp11;
        assert(tmp21 == tmp21);
      const auto tmp22 = -0.5 * tmp16;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp11 * tmp11;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 * tmp1;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp11 * tmp16;
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp25 * tmp1;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp16 * tmp16;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 * tmp1;
        assert(tmp28 == tmp28);
      return std::make_pair( [ tmp10, tmp13, tmp15, tmp18, tmp20, tmp21, tmp22, tmp4, tmp6, tmp8, tmp9 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { 0 }, { { ((tmp4 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp6 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ]) + tmp9 * (((std::get< 1 >( phiIn ))[ 0 ])[ 0 ])[ 0 ]) + tmp9 * (((std::get< 1 >( phiIn ))[ 0 ])[ 1 ])[ 1 ], ((tmp6 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp8 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ]) + tmp10 * (((std::get< 1 >( phiIn ))[ 0 ])[ 0 ])[ 0 ]) + tmp10 * (((std::get< 1 >( phiIn ))[ 0 ])[ 1 ])[ 1 ] } }, { { { tmp9 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp10 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ], 0 }, { 0, tmp9 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp10 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ] } } } }, RangeValueType{ { 0 }, { { ((tmp13 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp15 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ]) + tmp21 * (((std::get< 1 >( phiIn ))[ 0 ])[ 0 ])[ 0 ]) + tmp21 * (((std::get< 1 >( phiIn ))[ 0 ])[ 1 ])[ 1 ], ((tmp18 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp20 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ]) + tmp22 * (((std::get< 1 >( phiIn ))[ 0 ])[ 0 ])[ 0 ]) + tmp22 * (((std::get< 1 >( phiIn ))[ 0 ])[ 1 ])[ 1 ] } }, { { { tmp9 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp10 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ], 0 }, { 0, tmp9 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp10 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ] } } } } );
        }, [ tmp10, tmp13, tmp15, tmp18, tmp20, tmp21, tmp22, tmp24, tmp26, tmp28, tmp9 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { 0 }, { { ((tmp13 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp18 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ]) + tmp9 * (((std::get< 1 >( phiOut ))[ 0 ])[ 0 ])[ 0 ]) + tmp9 * (((std::get< 1 >( phiOut ))[ 0 ])[ 1 ])[ 1 ], ((tmp15 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp20 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ]) + tmp10 * (((std::get< 1 >( phiOut ))[ 0 ])[ 0 ])[ 0 ]) + tmp10 * (((std::get< 1 >( phiOut ))[ 0 ])[ 1 ])[ 1 ] } }, { { { tmp21 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp22 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ], 0 }, { 0, tmp21 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp22 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ] } } } }, RangeValueType{ { 0 }, { { ((tmp24 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp26 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ]) + tmp21 * (((std::get< 1 >( phiOut ))[ 0 ])[ 0 ])[ 0 ]) + tmp21 * (((std::get< 1 >( phiOut ))[ 0 ])[ 1 ])[ 1 ], ((tmp26 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp28 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ]) + tmp22 * (((std::get< 1 >( phiOut ))[ 0 ])[ 0 ])[ 0 ]) + tmp22 * (((std::get< 1 >( phiOut ))[ 0 ])[ 1 ])[ 1 ] } }, { { { tmp21 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp22 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ], 0 }, { 0, tmp21 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp22 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ] } } } } );
        } );
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
          using std::pow;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
            assert(tmp1 == tmp1);
          const auto tmp2 = std::sin( tmp1 );
            assert(tmp2 == tmp2);
          const auto tmp3 = std::pow( tmp2, 2 );
            assert(tmp3 == tmp3);
          const auto tmp4 = 3.141592653589793 * tmp0[ 0 ];
            assert(tmp4 == tmp4);
          const auto tmp5 = std::sin( tmp4 );
            assert(tmp5 == tmp5);
          const auto tmp6 = std::pow( tmp5, 2 );
            assert(tmp6 == tmp6);
          const auto tmp7 = tmp6 * tmp3;
            assert(tmp7 == tmp7);
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
          using std::pow;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
            assert(tmp1 == tmp1);
          const auto tmp2 = std::sin( tmp1 );
            assert(tmp2 == tmp2);
          const auto tmp3 = std::pow( tmp2, 2 );
            assert(tmp3 == tmp3);
          const auto tmp4 = 3.141592653589793 * tmp0[ 0 ];
            assert(tmp4 == tmp4);
          const auto tmp5 = std::sin( tmp4 );
            assert(tmp5 == tmp5);
          const auto tmp6 = std::cos( tmp4 );
            assert(tmp6 == tmp6);
          const auto tmp7 = 3.141592653589793 * tmp6;
            assert(tmp7 == tmp7);
          const auto tmp8 = 2 * tmp7;
            assert(tmp8 == tmp8);
          const auto tmp9 = tmp8 * tmp5;
            assert(tmp9 == tmp9);
          const auto tmp10 = tmp9 * tmp3;
            assert(tmp10 == tmp10);
          const auto tmp11 = std::pow( tmp5, 2 );
            assert(tmp11 == tmp11);
          const auto tmp12 = std::cos( tmp1 );
            assert(tmp12 == tmp12);
          const auto tmp13 = 3.141592653589793 * tmp12;
            assert(tmp13 == tmp13);
          const auto tmp14 = 2 * tmp13;
            assert(tmp14 == tmp14);
          const auto tmp15 = tmp14 * tmp2;
            assert(tmp15 == tmp15);
          const auto tmp16 = tmp15 * tmp11;
            assert(tmp16 == tmp16);
          (result[ 0 ])[ 0 ] = tmp10;
          (result[ 0 ])[ 1 ] = tmp16;
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
      return entity_[ (static_cast< std::size_t >( Side::in )) ];
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }
    std::array< EntityType, 2 > entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    ctype facetArea_;
    ConstantTupleType constants_;
    std::array< std::tuple<  >, 2 > coefficients_;
  };

} // namespace Integrands_1896a171b906827df87904471d1a8fa6v1_4

PYBIND11_MODULE( integrands_1896a171b906827df87904471d1a8fa6v1_4_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > GridPart;
  typedef Integrands_1896a171b906827df87904471d1a8fa6v1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_1896a171b906827df87904471d1a8fa6v1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_1896a171b906827df87904471d1a8fa6v1_4_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_1896a171b906827df87904471d1a8fa6v1_4
