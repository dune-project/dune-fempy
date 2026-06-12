#ifndef GuardIntegrands_6df37dd8a3c7afefb581db5ae9568645v1_4
#define GuardIntegrands_6df37dd8a3c7afefb581db5ae9568645v1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <utility>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_6df37dd8a3c7afefb581db5ae9568645v1_4
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffsolution, class Coeffsolution_old, class Coeffiterate >
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
    typedef double ContimeStep;
    typedef double Conpenalty;
    typedef std::tuple< std::shared_ptr< ContimeStep >, std::shared_ptr< Conpenalty > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffsolution, Coeffsolution_old, Coeffiterate > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsolution>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsolution_old>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffiterate>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffsolution &coeffsolution, const Coeffsolution_old &coeffsolution_old, const Coeffiterate &coeffiterate, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_{{ std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution >, Dune::Fem::ConstLocalFunction< Coeffsolution_old >, Dune::Fem::ConstLocalFunction< Coeffiterate > >( Dune::Fem::ConstLocalFunction< Coeffsolution >( coeffsolution ), Dune::Fem::ConstLocalFunction< Coeffsolution_old >( coeffsolution_old ), Dune::Fem::ConstLocalFunction< Coeffiterate >( coeffiterate ) ), std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution >, Dune::Fem::ConstLocalFunction< Coeffsolution_old >, Dune::Fem::ConstLocalFunction< Coeffiterate > >( Dune::Fem::ConstLocalFunction< Coeffsolution >( coeffsolution ), Dune::Fem::ConstLocalFunction< Coeffsolution_old >( coeffsolution_old ), Dune::Fem::ConstLocalFunction< Coeffiterate >( coeffiterate ) ) }}
    {
      std::get< 0 >( constants_ ) = std::make_shared< ContimeStep >( (ContimeStep(0)) );
      std::get< 1 >( constants_ ) = std::make_shared< Conpenalty >( (Conpenalty(0)) );
    }

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< ContimeStep >( (constant<0>()) );
      std::get< 1 >( constants_ ) = std::make_shared< Conpenalty >( (constant<1>()) );
    }

    bool init ( const EntityType &entity )
    {
      entity_[ (static_cast< std::size_t >( Side::in )) ] = entity;
      geometry_.emplace( this->entity().geometry() );
      cellVolume_[ (static_cast< std::size_t >( Side::in )) ] = (geometry().volume());
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
    }

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      const auto& intersectionGeom = intersection_.geometry();
      facetArea_ = (intersectionGeom.volume());
      entity_[ (static_cast< std::size_t >( Side::in )) ] = (intersection.inside());
      geometry_.emplace( this->entity().geometry() );
      cellVolume_[ (static_cast< std::size_t >( Side::in )) ] = (geometry().volume());
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
        const Geometry& outsideGeom = entity_[ static_cast< std::size_t >( Side::out ) ].geometry();
        cellVolume_[ static_cast< std::size_t >( Side::out ) ] = outsideGeom.volume();
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
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
      using std::abs;
      using std::pow;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = -0.45 + tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::abs( tmp1 );
        assert(tmp2 == tmp2);
      const auto tmp3 = -0.49 + tmp0[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = std::abs( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = (tmp4 < 0.03 ? 1.0 : 0.0) * (tmp2 < 0.11 ? 1.0 : 0.0);
        assert(tmp5 == tmp5);
      const auto tmp6 = -0.010000000000000009 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 0.4 + tmp6;
        assert(tmp7 == tmp7);
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp8 = evaluateCoefficient< 1 >( x );
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8[ 1 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = std::get< 0 >( u );
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10[ 0 ] + tmp9;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 * tmp7;
        assert(tmp12 == tmp12);
      double tmp13 = constant< 0 >();
        assert(tmp13 == tmp13);
      const auto tmp14 = -0.01999999999999999 * tmp5;
        assert(tmp14 == tmp14);
      const auto tmp15 = 0.12 + tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = 1.0 + tmp16;
        assert(tmp17 == tmp17);
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp18 = evaluateCoefficient< 2 >( x );
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp18[ 1 ] / tmp17;
        assert(tmp19 == tmp19);
      const auto tmp20 = std::pow( tmp19, 2 );
        assert(tmp20 == tmp20);
      const auto tmp21 = -0.7000000000000002 * tmp5;
        assert(tmp21 == tmp21);
      const auto tmp22 = 2.7 + tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = 2.0 + tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 / tmp22;
        assert(tmp24 == tmp24);
      const auto tmp25 = -1 * tmp18[ 1 ];
        assert(tmp25 == tmp25);
      const auto tmp26 = 1 + tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 + tmp16;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 / tmp17;
        assert(tmp28 == tmp28);
      const auto tmp29 = std::pow( tmp28, tmp24 );
        assert(tmp29 == tmp29);
      const auto tmp30 = -1 * tmp29;
        assert(tmp30 == tmp30);
      const auto tmp31 = 1.0 + tmp30;
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 * tmp20;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 / 0.0009;
        assert(tmp33 == tmp33);
      const auto tmp34 = 5e-11 * tmp5;
        assert(tmp34 == tmp34);
      const auto tmp35 = -5e-11 + tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 * tmp33;
        assert(tmp36 == tmp36);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp37 = jacobianCoefficient< 0 >( x );
        assert(tmp37 == tmp37);
      const auto tmp38 = 14322.6 + (tmp37[ 0 ])[ 1 ];
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 * tmp36;
        assert(tmp39 == tmp39);
      const auto tmp40 = -9.993360000000001e-11 * tmp5;
        assert(tmp40 == tmp40);
      const auto tmp41 = 1e-10 + tmp40;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp41 * tmp33;
        assert(tmp42 == tmp42);
      const auto tmp43 = (tmp37[ 0 ])[ 0 ] * tmp42;
        assert(tmp43 == tmp43);
      const auto tmp44 = tmp43 + tmp39;
        assert(tmp44 == tmp44);
      const auto tmp45 = -1.0 / tmp17;
        assert(tmp45 == tmp45);
      const auto tmp46 = -1.0 / tmp22;
        assert(tmp46 == tmp46);
      const auto tmp47 = -1.0 + tmp46;
        assert(tmp47 == tmp47);
      const auto tmp48 = std::pow( tmp28, tmp47 );
        assert(tmp48 == tmp48);
      const auto tmp49 = 4245.0 * tmp5;
        assert(tmp49 == tmp49);
      const auto tmp50 = 755.0 + tmp49;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp50 * tmp46;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 * tmp48;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp52 * tmp45;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp36 * tmp53;
        assert(tmp54 == tmp54);
      const auto tmp55 = std::get< 1 >( u );
        assert(tmp55 == tmp55);
      const auto tmp56 = (tmp55[ 0 ])[ 1 ] * tmp54;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp42 * tmp53;
        assert(tmp57 == tmp57);
      const auto tmp58 = (tmp55[ 0 ])[ 0 ] * tmp57;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp58 + tmp56;
        assert(tmp59 == tmp59);
      const auto tmp60 = tmp59 + tmp44;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp60 * tmp13;
        assert(tmp61 == tmp61);
      const auto tmp62 = (tmp37[ 0 ])[ 0 ] * tmp36;
        assert(tmp62 == tmp62);
      const auto tmp63 = tmp38 * tmp42;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp63 + tmp62;
        assert(tmp64 == tmp64);
      const auto tmp65 = (tmp55[ 0 ])[ 0 ] * tmp54;
        assert(tmp65 == tmp65);
      const auto tmp66 = (tmp55[ 0 ])[ 1 ] * tmp57;
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp66 + tmp65;
        assert(tmp67 == tmp67);
      const auto tmp68 = tmp67 + tmp64;
        assert(tmp68 == tmp68);
      const auto tmp69 = tmp68 * tmp13;
        assert(tmp69 == tmp69);
      return RangeValueType{ { tmp12 }, { { tmp61, tmp69 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = -0.45 + tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::abs( tmp1 );
        assert(tmp2 == tmp2);
      const auto tmp3 = -0.49 + tmp0[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = std::abs( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = (tmp4 < 0.03 ? 1.0 : 0.0) * (tmp2 < 0.11 ? 1.0 : 0.0);
        assert(tmp5 == tmp5);
      const auto tmp6 = -0.010000000000000009 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 0.4 + tmp6;
        assert(tmp7 == tmp7);
      double tmp8 = constant< 0 >();
        assert(tmp8 == tmp8);
      const auto tmp9 = -0.01999999999999999 * tmp5;
        assert(tmp9 == tmp9);
      const auto tmp10 = 0.12 + tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = -1 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = 1.0 + tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1.0 / tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = -0.7000000000000002 * tmp5;
        assert(tmp14 == tmp14);
      const auto tmp15 = 2.7 + tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1.0 / tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = -1.0 + tmp16;
        assert(tmp17 == tmp17);
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp18 = evaluateCoefficient< 2 >( x );
        assert(tmp18 == tmp18);
      const auto tmp19 = -1 * tmp18[ 1 ];
        assert(tmp19 == tmp19);
      const auto tmp20 = 1 + tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 + tmp11;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 / tmp12;
        assert(tmp22 == tmp22);
      const auto tmp23 = std::pow( tmp22, tmp17 );
        assert(tmp23 == tmp23);
      const auto tmp24 = 4245.0 * tmp5;
        assert(tmp24 == tmp24);
      const auto tmp25 = 755.0 + tmp24;
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp25 * tmp16;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 * tmp23;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 * tmp13;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp18[ 1 ] / tmp12;
        assert(tmp29 == tmp29);
      const auto tmp30 = std::pow( tmp29, 2 );
        assert(tmp30 == tmp30);
      const auto tmp31 = 2.0 + tmp15;
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 / tmp15;
        assert(tmp32 == tmp32);
      const auto tmp33 = std::pow( tmp22, tmp32 );
        assert(tmp33 == tmp33);
      const auto tmp34 = -1 * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = 1.0 + tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 * tmp30;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 / 0.0009;
        assert(tmp37 == tmp37);
      const auto tmp38 = -9.993360000000001e-11 * tmp5;
        assert(tmp38 == tmp38);
      const auto tmp39 = 1e-10 + tmp38;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 * tmp37;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 * tmp28;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp41 * tmp8;
        assert(tmp42 == tmp42);
      const auto tmp43 = 5e-11 * tmp5;
        assert(tmp43 == tmp43);
      const auto tmp44 = -5e-11 + tmp43;
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp44 * tmp37;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 * tmp28;
        assert(tmp46 == tmp46);
      const auto tmp47 = tmp46 * tmp8;
        assert(tmp47 == tmp47);
      return [ tmp42, tmp47, tmp7 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp7 * (std::get< 0 >( phi ))[ 0 ] }, { { tmp42 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp47 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], tmp47 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp42 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      GlobalCoordinateType tmp1 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp1 == tmp1);
      const auto tmp2 = -0.9 + tmp1[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = std::abs( tmp2 );
        assert(tmp3 == tmp3);
      const auto tmp4 = std::abs( tmp1[ 0 ] );
        assert(tmp4 == tmp4);
      const auto tmp5 = (tmp4 < 1e-08 ? 1.0 : 0.0) + (tmp3 < 1e-08 ? 1.0 : 0.0);
        assert(tmp5 == tmp5);
      const auto tmp6 = -0.45 + tmp1[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = std::abs( tmp6 );
        assert(tmp7 == tmp7);
      const auto tmp8 = -0.49 + tmp1[ 1 ];
        assert(tmp8 == tmp8);
      const auto tmp9 = std::abs( tmp8 );
        assert(tmp9 == tmp9);
      const auto tmp10 = (tmp9 < 0.03 ? 1.0 : 0.0) * (tmp7 < 0.11 ? 1.0 : 0.0);
        assert(tmp10 == tmp10);
      const auto tmp11 = -0.01999999999999999 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = 0.12 + tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = 1.0 + tmp13;
        assert(tmp14 == tmp14);
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp15 = evaluateCoefficient< 2 >( x );
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15[ 1 ] / tmp14;
        assert(tmp16 == tmp16);
      const auto tmp17 = std::pow( tmp16, 2 );
        assert(tmp17 == tmp17);
      const auto tmp18 = -0.7000000000000002 * tmp10;
        assert(tmp18 == tmp18);
      const auto tmp19 = 2.7 + tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = 2.0 + tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 / tmp19;
        assert(tmp21 == tmp21);
      const auto tmp22 = -1 * tmp15[ 1 ];
        assert(tmp22 == tmp22);
      const auto tmp23 = 1 + tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 + tmp13;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 / tmp14;
        assert(tmp25 == tmp25);
      const auto tmp26 = std::pow( tmp25, tmp21 );
        assert(tmp26 == tmp26);
      const auto tmp27 = -1 * tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = 1.0 + tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 * tmp17;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 / 0.0009;
        assert(tmp30 == tmp30);
      const auto tmp31 = 5e-11 * tmp10;
        assert(tmp31 == tmp31);
      const auto tmp32 = -5e-11 + tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 * tmp30;
        assert(tmp33 == tmp33);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp34 = jacobianCoefficient< 0 >( x );
        assert(tmp34 == tmp34);
      const auto tmp35 = 14322.6 + (tmp34[ 0 ])[ 1 ];
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 * tmp33;
        assert(tmp36 == tmp36);
      const auto tmp37 = -9.993360000000001e-11 * tmp10;
        assert(tmp37 == tmp37);
      const auto tmp38 = 1e-10 + tmp37;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 * tmp30;
        assert(tmp39 == tmp39);
      const auto tmp40 = (tmp34[ 0 ])[ 0 ] * tmp39;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 + tmp36;
        assert(tmp41 == tmp41);
      const auto tmp42 = -1.0 / tmp14;
        assert(tmp42 == tmp42);
      const auto tmp43 = -1.0 / tmp19;
        assert(tmp43 == tmp43);
      const auto tmp44 = -1.0 + tmp43;
        assert(tmp44 == tmp44);
      const auto tmp45 = std::pow( tmp25, tmp44 );
        assert(tmp45 == tmp45);
      const auto tmp46 = 4245.0 * tmp10;
        assert(tmp46 == tmp46);
      const auto tmp47 = 755.0 + tmp46;
        assert(tmp47 == tmp47);
      const auto tmp48 = tmp47 * tmp43;
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp48 * tmp45;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp49 * tmp42;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp33 * tmp50;
        assert(tmp51 == tmp51);
      const auto tmp52 = std::get< 1 >( u );
        assert(tmp52 == tmp52);
      const auto tmp53 = (tmp52[ 0 ])[ 1 ] * tmp51;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp39 * tmp50;
        assert(tmp54 == tmp54);
      const auto tmp55 = (tmp52[ 0 ])[ 0 ] * tmp54;
        assert(tmp55 == tmp55);
      const auto tmp56 = tmp55 + tmp53;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp56 + tmp41;
        assert(tmp57 == tmp57);
      GlobalCoordinateType tmp58 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp58[ 0 ] * tmp57;
        assert(tmp59 == tmp59);
      const auto tmp60 = (tmp34[ 0 ])[ 0 ] * tmp33;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp35 * tmp39;
        assert(tmp61 == tmp61);
      const auto tmp62 = tmp61 + tmp60;
        assert(tmp62 == tmp62);
      const auto tmp63 = (tmp52[ 0 ])[ 0 ] * tmp51;
        assert(tmp63 == tmp63);
      const auto tmp64 = (tmp52[ 0 ])[ 1 ] * tmp54;
        assert(tmp64 == tmp64);
      const auto tmp65 = tmp64 + tmp63;
        assert(tmp65 == tmp65);
      const auto tmp66 = tmp65 + tmp62;
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp58[ 1 ] * tmp66;
        assert(tmp67 == tmp67);
      const auto tmp68 = tmp67 + tmp59;
        assert(tmp68 == tmp68);
      const auto tmp69 = tmp68 * tmp5;
        assert(tmp69 == tmp69);
      const auto tmp70 = -1 * tmp69;
        assert(tmp70 == tmp70);
      const auto tmp71 = tmp70 * tmp0;
        assert(tmp71 == tmp71);
      const auto tmp72 = -0.65 + tmp1[ 1 ];
        assert(tmp72 == tmp72);
      const auto tmp73 = std::abs( tmp72 );
        assert(tmp73 == tmp73);
      const auto tmp74 = (tmp7 < 0.06 ? 1.0 : 0.0) * (tmp73 < 1e-08 ? 1.0 : 0.0);
        assert(tmp74 == tmp74);
      const auto tmp75 = -5.137e-05 * tmp74;
        assert(tmp75 == tmp75);
      const auto tmp76 = tmp75 * tmp0;
        assert(tmp76 == tmp76);
      auto tmp77 = facetArea_;
        assert(tmp77 == tmp77);
      auto tmp78 = cellVolume();
        assert(tmp78 == tmp78);
      const auto tmp79 = tmp78 / tmp77;
        assert(tmp79 == tmp79);
      double tmp80 = constant< 1 >();
        assert(tmp80 == tmp80);
      const auto tmp81 = tmp58[ 1 ] * tmp32;
        assert(tmp81 == tmp81);
      const auto tmp82 = tmp58[ 0 ] * tmp38;
        assert(tmp82 == tmp82);
      const auto tmp83 = tmp82 + tmp81;
        assert(tmp83 == tmp83);
      const auto tmp84 = tmp58[ 0 ] * tmp83;
        assert(tmp84 == tmp84);
      const auto tmp85 = tmp58[ 0 ] * tmp32;
        assert(tmp85 == tmp85);
      const auto tmp86 = tmp58[ 1 ] * tmp38;
        assert(tmp86 == tmp86);
      const auto tmp87 = tmp86 + tmp85;
        assert(tmp87 == tmp87);
      const auto tmp88 = tmp58[ 1 ] * tmp87;
        assert(tmp88 == tmp88);
      const auto tmp89 = tmp88 + tmp84;
        assert(tmp89 == tmp89);
      const auto tmp90 = tmp89 * tmp80;
        assert(tmp90 == tmp90);
      const auto tmp91 = 0.5 / tmp14;
        assert(tmp91 == tmp91);
      const auto tmp92 = std::pow( tmp91, 2 );
        assert(tmp92 == tmp92);
      const auto tmp93 = 0.5 + tmp13;
        assert(tmp93 == tmp93);
      const auto tmp94 = tmp93 / tmp14;
        assert(tmp94 == tmp94);
      const auto tmp95 = std::pow( tmp94, tmp21 );
        assert(tmp95 == tmp95);
      const auto tmp96 = -1 * tmp95;
        assert(tmp96 == tmp96);
      const auto tmp97 = 1.0 + tmp96;
        assert(tmp97 == tmp97);
      const auto tmp98 = tmp97 * tmp92;
        assert(tmp98 == tmp98);
      const auto tmp99 = tmp98 / 0.0009;
        assert(tmp99 == tmp99);
      const auto tmp100 = std::pow( tmp94, tmp44 );
        assert(tmp100 == tmp100);
      const auto tmp101 = tmp48 * tmp100;
        assert(tmp101 == tmp101);
      const auto tmp102 = tmp101 * tmp42;
        assert(tmp102 == tmp102);
      const auto tmp103 = tmp102 * tmp99;
        assert(tmp103 == tmp103);
      const auto tmp104 = tmp103 * tmp90;
        assert(tmp104 == tmp104);
      const auto tmp105 = tmp104 / tmp79;
        assert(tmp105 == tmp105);
      const auto tmp106 = std::get< 0 >( u );
        assert(tmp106 == tmp106);
      const auto tmp107 = tmp106[ 0 ] * tmp105;
        assert(tmp107 == tmp107);
      const auto tmp108 = tmp5 * tmp107;
        assert(tmp108 == tmp108);
      const auto tmp109 = tmp108 * tmp0;
        assert(tmp109 == tmp109);
      const auto tmp110 = tmp109 + tmp76;
        assert(tmp110 == tmp110);
      const auto tmp111 = tmp110 + tmp71;
        assert(tmp111 == tmp111);
      return RangeValueType{ { tmp111 }, { { 0, 0 } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      auto tmp1 = facetArea_;
        assert(tmp1 == tmp1);
      auto tmp2 = cellVolume();
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 / tmp1;
        assert(tmp3 == tmp3);
      double tmp4 = constant< 1 >();
        assert(tmp4 == tmp4);
      GlobalCoordinateType tmp5 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp5 == tmp5);
      const auto tmp6 = -0.45 + tmp5[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = std::abs( tmp6 );
        assert(tmp7 == tmp7);
      const auto tmp8 = -0.49 + tmp5[ 1 ];
        assert(tmp8 == tmp8);
      const auto tmp9 = std::abs( tmp8 );
        assert(tmp9 == tmp9);
      const auto tmp10 = (tmp9 < 0.03 ? 1.0 : 0.0) * (tmp7 < 0.11 ? 1.0 : 0.0);
        assert(tmp10 == tmp10);
      const auto tmp11 = 5e-11 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = -5e-11 + tmp11;
        assert(tmp12 == tmp12);
      GlobalCoordinateType tmp13 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp13[ 1 ] * tmp12;
        assert(tmp14 == tmp14);
      const auto tmp15 = -9.993360000000001e-11 * tmp10;
        assert(tmp15 == tmp15);
      const auto tmp16 = 1e-10 + tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp13[ 0 ] * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 + tmp14;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp13[ 0 ] * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp13[ 0 ] * tmp12;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp13[ 1 ] * tmp16;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 + tmp20;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp13[ 1 ] * tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 + tmp19;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 * tmp4;
        assert(tmp25 == tmp25);
      const auto tmp26 = -0.01999999999999999 * tmp10;
        assert(tmp26 == tmp26);
      const auto tmp27 = 0.12 + tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = -1 * tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = 1.0 + tmp28;
        assert(tmp29 == tmp29);
      const auto tmp30 = 0.5 / tmp29;
        assert(tmp30 == tmp30);
      const auto tmp31 = std::pow( tmp30, 2 );
        assert(tmp31 == tmp31);
      const auto tmp32 = -0.7000000000000002 * tmp10;
        assert(tmp32 == tmp32);
      const auto tmp33 = 2.7 + tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = 2.0 + tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp34 / tmp33;
        assert(tmp35 == tmp35);
      const auto tmp36 = 0.5 + tmp28;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 / tmp29;
        assert(tmp37 == tmp37);
      const auto tmp38 = std::pow( tmp37, tmp35 );
        assert(tmp38 == tmp38);
      const auto tmp39 = -1 * tmp38;
        assert(tmp39 == tmp39);
      const auto tmp40 = 1.0 + tmp39;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 * tmp31;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp41 / 0.0009;
        assert(tmp42 == tmp42);
      const auto tmp43 = -1.0 / tmp29;
        assert(tmp43 == tmp43);
      const auto tmp44 = -1.0 / tmp33;
        assert(tmp44 == tmp44);
      const auto tmp45 = -1.0 + tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = std::pow( tmp37, tmp45 );
        assert(tmp46 == tmp46);
      const auto tmp47 = 4245.0 * tmp10;
        assert(tmp47 == tmp47);
      const auto tmp48 = 755.0 + tmp47;
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp48 * tmp44;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp49 * tmp46;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp50 * tmp43;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 * tmp42;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp52 * tmp25;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 / tmp3;
        assert(tmp54 == tmp54);
      const auto tmp55 = -0.9 + tmp5[ 0 ];
        assert(tmp55 == tmp55);
      const auto tmp56 = std::abs( tmp55 );
        assert(tmp56 == tmp56);
      const auto tmp57 = std::abs( tmp5[ 0 ] );
        assert(tmp57 == tmp57);
      const auto tmp58 = (tmp57 < 1e-08 ? 1.0 : 0.0) + (tmp56 < 1e-08 ? 1.0 : 0.0);
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp58 * tmp54;
        assert(tmp59 == tmp59);
      const auto tmp60 = tmp59 * tmp0;
        assert(tmp60 == tmp60);
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp61 = evaluateCoefficient< 2 >( x );
        assert(tmp61 == tmp61);
      const auto tmp62 = -1 * tmp61[ 1 ];
        assert(tmp62 == tmp62);
      const auto tmp63 = 1 + tmp62;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp63 + tmp28;
        assert(tmp64 == tmp64);
      const auto tmp65 = tmp64 / tmp29;
        assert(tmp65 == tmp65);
      const auto tmp66 = std::pow( tmp65, tmp45 );
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp49 * tmp66;
        assert(tmp67 == tmp67);
      const auto tmp68 = tmp67 * tmp43;
        assert(tmp68 == tmp68);
      const auto tmp69 = tmp61[ 1 ] / tmp29;
        assert(tmp69 == tmp69);
      const auto tmp70 = std::pow( tmp69, 2 );
        assert(tmp70 == tmp70);
      const auto tmp71 = std::pow( tmp65, tmp35 );
        assert(tmp71 == tmp71);
      const auto tmp72 = -1 * tmp71;
        assert(tmp72 == tmp72);
      const auto tmp73 = 1.0 + tmp72;
        assert(tmp73 == tmp73);
      const auto tmp74 = tmp73 * tmp70;
        assert(tmp74 == tmp74);
      const auto tmp75 = tmp74 / 0.0009;
        assert(tmp75 == tmp75);
      const auto tmp76 = tmp12 * tmp75;
        assert(tmp76 == tmp76);
      const auto tmp77 = tmp76 * tmp68;
        assert(tmp77 == tmp77);
      const auto tmp78 = tmp13[ 1 ] * tmp77;
        assert(tmp78 == tmp78);
      const auto tmp79 = tmp16 * tmp75;
        assert(tmp79 == tmp79);
      const auto tmp80 = tmp79 * tmp68;
        assert(tmp80 == tmp80);
      const auto tmp81 = tmp13[ 0 ] * tmp80;
        assert(tmp81 == tmp81);
      const auto tmp82 = tmp81 + tmp78;
        assert(tmp82 == tmp82);
      const auto tmp83 = tmp82 * tmp58;
        assert(tmp83 == tmp83);
      const auto tmp84 = -1 * tmp83;
        assert(tmp84 == tmp84);
      const auto tmp85 = tmp84 * tmp0;
        assert(tmp85 == tmp85);
      const auto tmp86 = tmp13[ 0 ] * tmp77;
        assert(tmp86 == tmp86);
      const auto tmp87 = tmp13[ 1 ] * tmp80;
        assert(tmp87 == tmp87);
      const auto tmp88 = tmp87 + tmp86;
        assert(tmp88 == tmp88);
      const auto tmp89 = tmp88 * tmp58;
        assert(tmp89 == tmp89);
      const auto tmp90 = -1 * tmp89;
        assert(tmp90 == tmp90);
      const auto tmp91 = tmp90 * tmp0;
        assert(tmp91 == tmp91);
      return [ tmp60, tmp85, tmp91 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (tmp60 * (std::get< 0 >( phi ))[ 0 ] + tmp85 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp91 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { { 0, 0 } } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      using std::max;
      using std::pow;
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      auto tmp1 = facetArea_;
        assert(tmp1 == tmp1);
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
        assert(tmp2 == tmp2);
      auto tmp3 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 + tmp2;
        assert(tmp4 == tmp4);
      const auto tmp5 = 0.5 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 / tmp1;
        assert(tmp6 == tmp6);
      GlobalCoordinateType tmp7 = geometry().global( Dune::Fem::coordinate( xIn ) );
        assert(tmp7 == tmp7);
      const auto tmp8 = -0.45 + tmp7[ 0 ];
        assert(tmp8 == tmp8);
      const auto tmp9 = std::abs( tmp8 );
        assert(tmp9 == tmp9);
      const auto tmp10 = -0.49 + tmp7[ 1 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = std::abs( tmp10 );
        assert(tmp11 == tmp11);
      const auto tmp12 = (tmp11 < 0.03 ? 1.0 : 0.0) * (tmp9 < 0.11 ? 1.0 : 0.0);
        assert(tmp12 == tmp12);
      const auto tmp13 = -0.01999999999999999 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = 0.12 + tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = -1 * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = 1.0 + tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = 0.5 / tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = std::pow( tmp17, 2 );
        assert(tmp18 == tmp18);
      const auto tmp19 = -0.7000000000000002 * tmp12;
        assert(tmp19 == tmp19);
      const auto tmp20 = 2.7 + tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = 2.0 + tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 / tmp20;
        assert(tmp22 == tmp22);
      const auto tmp23 = 0.5 + tmp15;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 / tmp16;
        assert(tmp24 == tmp24);
      const auto tmp25 = std::pow( tmp24, tmp22 );
        assert(tmp25 == tmp25);
      const auto tmp26 = -1 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = 1.0 + tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 * tmp18;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 / 0.0009;
        assert(tmp29 == tmp29);
      const auto tmp30 = -1.0 / tmp16;
        assert(tmp30 == tmp30);
      const auto tmp31 = -1.0 / tmp20;
        assert(tmp31 == tmp31);
      const auto tmp32 = -1.0 + tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = std::pow( tmp24, tmp32 );
        assert(tmp33 == tmp33);
      const auto tmp34 = 4245.0 * tmp12;
        assert(tmp34 == tmp34);
      const auto tmp35 = 755.0 + tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 * tmp31;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 * tmp33;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp37 * tmp30;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 * tmp29;
        assert(tmp39 == tmp39);
      const auto tmp40 = std::max( tmp39, tmp39 );
        assert(tmp40 == tmp40);
      double tmp41 = constant< 1 >();
        assert(tmp41 == tmp41);
      GlobalCoordinateType tmp42 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp42 == tmp42);
      const auto tmp43 = -1 * tmp42[ 1 ];
        assert(tmp43 == tmp43);
      const auto tmp44 = -9.993360000000001e-11 * tmp12;
        assert(tmp44 == tmp44);
      const auto tmp45 = 1e-10 + tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 * tmp43;
        assert(tmp46 == tmp46);
      const auto tmp47 = -1 * tmp42[ 0 ];
        assert(tmp47 == tmp47);
      const auto tmp48 = 5e-11 * tmp12;
        assert(tmp48 == tmp48);
      const auto tmp49 = -5e-11 + tmp48;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp49 * tmp47;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp50 + tmp46;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 * tmp43;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp49 * tmp43;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp45 * tmp47;
        assert(tmp54 == tmp54);
      const auto tmp55 = tmp54 + tmp53;
        assert(tmp55 == tmp55);
      const auto tmp56 = tmp55 * tmp47;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp56 + tmp52;
        assert(tmp57 == tmp57);
      const auto tmp58 = tmp42[ 1 ] * tmp49;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp42[ 0 ] * tmp45;
        assert(tmp59 == tmp59);
      const auto tmp60 = tmp59 + tmp58;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp42[ 0 ] * tmp60;
        assert(tmp61 == tmp61);
      const auto tmp62 = tmp42[ 0 ] * tmp49;
        assert(tmp62 == tmp62);
      const auto tmp63 = tmp42[ 1 ] * tmp45;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp63 + tmp62;
        assert(tmp64 == tmp64);
      const auto tmp65 = tmp42[ 1 ] * tmp64;
        assert(tmp65 == tmp65);
      const auto tmp66 = tmp65 + tmp61;
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp66 + tmp57;
        assert(tmp67 == tmp67);
      const auto tmp68 = 0.5 * tmp67;
        assert(tmp68 == tmp68);
      const auto tmp69 = tmp66 * tmp57;
        assert(tmp69 == tmp69);
      const auto tmp70 = tmp69 / tmp68;
        assert(tmp70 == tmp70);
      const auto tmp71 = tmp70 * tmp41;
        assert(tmp71 == tmp71);
      const auto tmp72 = tmp71 * tmp40;
        assert(tmp72 == tmp72);
      const auto tmp73 = tmp72 / tmp6;
        assert(tmp73 == tmp73);
      const auto tmp74 = std::get< 0 >( uOut );
        assert(tmp74 == tmp74);
      const auto tmp75 = -1 * tmp74[ 0 ];
        assert(tmp75 == tmp75);
      const auto tmp76 = std::get< 0 >( uIn );
        assert(tmp76 == tmp76);
      const auto tmp77 = tmp76[ 0 ] + tmp75;
        assert(tmp77 == tmp77);
      const auto tmp78 = tmp77 * tmp73;
        assert(tmp78 == tmp78);
      const auto tmp79 = tmp78 * tmp0;
        assert(tmp79 == tmp79);
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp80 = evaluateCoefficient< 2, Side::in >( xIn );
        assert(tmp80 == tmp80);
      const auto tmp81 = tmp80[ 1 ] / tmp16;
        assert(tmp81 == tmp81);
      const auto tmp82 = std::pow( tmp81, 2 );
        assert(tmp82 == tmp82);
      const auto tmp83 = -1 * tmp80[ 1 ];
        assert(tmp83 == tmp83);
      const auto tmp84 = 1 + tmp83;
        assert(tmp84 == tmp84);
      const auto tmp85 = tmp84 + tmp15;
        assert(tmp85 == tmp85);
      const auto tmp86 = tmp85 / tmp16;
        assert(tmp86 == tmp86);
      const auto tmp87 = std::pow( tmp86, tmp22 );
        assert(tmp87 == tmp87);
      const auto tmp88 = -1 * tmp87;
        assert(tmp88 == tmp88);
      const auto tmp89 = 1.0 + tmp88;
        assert(tmp89 == tmp89);
      const auto tmp90 = tmp89 * tmp82;
        assert(tmp90 == tmp90);
      const auto tmp91 = tmp90 / 0.0009;
        assert(tmp91 == tmp91);
      const auto tmp92 = tmp49 * tmp91;
        assert(tmp92 == tmp92);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp93 = jacobianCoefficient< 0, Side::in >( xIn );
        assert(tmp93 == tmp93);
      const auto tmp94 = 14322.6 + (tmp93[ 0 ])[ 1 ];
        assert(tmp94 == tmp94);
      const auto tmp95 = tmp94 * tmp92;
        assert(tmp95 == tmp95);
      const auto tmp96 = tmp45 * tmp91;
        assert(tmp96 == tmp96);
      const auto tmp97 = (tmp93[ 0 ])[ 0 ] * tmp96;
        assert(tmp97 == tmp97);
      const auto tmp98 = tmp97 + tmp95;
        assert(tmp98 == tmp98);
      const auto tmp99 = std::pow( tmp86, tmp32 );
        assert(tmp99 == tmp99);
      const auto tmp100 = tmp36 * tmp99;
        assert(tmp100 == tmp100);
      const auto tmp101 = tmp100 * tmp30;
        assert(tmp101 == tmp101);
      const auto tmp102 = tmp92 * tmp101;
        assert(tmp102 == tmp102);
      const auto tmp103 = std::get< 1 >( uIn );
        assert(tmp103 == tmp103);
      const auto tmp104 = (tmp103[ 0 ])[ 1 ] * tmp102;
        assert(tmp104 == tmp104);
      const auto tmp105 = tmp96 * tmp101;
        assert(tmp105 == tmp105);
      const auto tmp106 = (tmp103[ 0 ])[ 0 ] * tmp105;
        assert(tmp106 == tmp106);
      const auto tmp107 = tmp106 + tmp104;
        assert(tmp107 == tmp107);
      const auto tmp108 = tmp107 + tmp98;
        assert(tmp108 == tmp108);
      const auto tmp109 = tmp108 * tmp57;
        assert(tmp109 == tmp109);
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp110 = evaluateCoefficient< 2, Side::out >( xOut );
        assert(tmp110 == tmp110);
      const auto tmp111 = tmp110[ 1 ] / tmp16;
        assert(tmp111 == tmp111);
      const auto tmp112 = std::pow( tmp111, 2 );
        assert(tmp112 == tmp112);
      const auto tmp113 = -1 * tmp110[ 1 ];
        assert(tmp113 == tmp113);
      const auto tmp114 = 1 + tmp113;
        assert(tmp114 == tmp114);
      const auto tmp115 = tmp114 + tmp15;
        assert(tmp115 == tmp115);
      const auto tmp116 = tmp115 / tmp16;
        assert(tmp116 == tmp116);
      const auto tmp117 = std::pow( tmp116, tmp22 );
        assert(tmp117 == tmp117);
      const auto tmp118 = -1 * tmp117;
        assert(tmp118 == tmp118);
      const auto tmp119 = 1.0 + tmp118;
        assert(tmp119 == tmp119);
      const auto tmp120 = tmp119 * tmp112;
        assert(tmp120 == tmp120);
      const auto tmp121 = tmp120 / 0.0009;
        assert(tmp121 == tmp121);
      const auto tmp122 = tmp49 * tmp121;
        assert(tmp122 == tmp122);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp123 = jacobianCoefficient< 0, Side::out >( xOut );
        assert(tmp123 == tmp123);
      const auto tmp124 = 14322.6 + (tmp123[ 0 ])[ 1 ];
        assert(tmp124 == tmp124);
      const auto tmp125 = tmp124 * tmp122;
        assert(tmp125 == tmp125);
      const auto tmp126 = tmp45 * tmp121;
        assert(tmp126 == tmp126);
      const auto tmp127 = (tmp123[ 0 ])[ 0 ] * tmp126;
        assert(tmp127 == tmp127);
      const auto tmp128 = tmp127 + tmp125;
        assert(tmp128 == tmp128);
      const auto tmp129 = std::pow( tmp116, tmp32 );
        assert(tmp129 == tmp129);
      const auto tmp130 = tmp36 * tmp129;
        assert(tmp130 == tmp130);
      const auto tmp131 = tmp130 * tmp30;
        assert(tmp131 == tmp131);
      const auto tmp132 = tmp122 * tmp131;
        assert(tmp132 == tmp132);
      const auto tmp133 = std::get< 1 >( uOut );
        assert(tmp133 == tmp133);
      const auto tmp134 = (tmp133[ 0 ])[ 1 ] * tmp132;
        assert(tmp134 == tmp134);
      const auto tmp135 = tmp126 * tmp131;
        assert(tmp135 == tmp135);
      const auto tmp136 = (tmp133[ 0 ])[ 0 ] * tmp135;
        assert(tmp136 == tmp136);
      const auto tmp137 = tmp136 + tmp134;
        assert(tmp137 == tmp137);
      const auto tmp138 = tmp137 + tmp128;
        assert(tmp138 == tmp138);
      const auto tmp139 = tmp138 * tmp66;
        assert(tmp139 == tmp139);
      const auto tmp140 = tmp139 + tmp109;
        assert(tmp140 == tmp140);
      const auto tmp141 = tmp140 / tmp67;
        assert(tmp141 == tmp141);
      const auto tmp142 = tmp42[ 0 ] * tmp141;
        assert(tmp142 == tmp142);
      const auto tmp143 = (tmp93[ 0 ])[ 0 ] * tmp92;
        assert(tmp143 == tmp143);
      const auto tmp144 = tmp94 * tmp96;
        assert(tmp144 == tmp144);
      const auto tmp145 = tmp144 + tmp143;
        assert(tmp145 == tmp145);
      const auto tmp146 = (tmp103[ 0 ])[ 0 ] * tmp102;
        assert(tmp146 == tmp146);
      const auto tmp147 = (tmp103[ 0 ])[ 1 ] * tmp105;
        assert(tmp147 == tmp147);
      const auto tmp148 = tmp147 + tmp146;
        assert(tmp148 == tmp148);
      const auto tmp149 = tmp148 + tmp145;
        assert(tmp149 == tmp149);
      const auto tmp150 = tmp149 * tmp57;
        assert(tmp150 == tmp150);
      const auto tmp151 = (tmp123[ 0 ])[ 0 ] * tmp122;
        assert(tmp151 == tmp151);
      const auto tmp152 = tmp124 * tmp126;
        assert(tmp152 == tmp152);
      const auto tmp153 = tmp152 + tmp151;
        assert(tmp153 == tmp153);
      const auto tmp154 = (tmp133[ 0 ])[ 0 ] * tmp132;
        assert(tmp154 == tmp154);
      const auto tmp155 = (tmp133[ 0 ])[ 1 ] * tmp135;
        assert(tmp155 == tmp155);
      const auto tmp156 = tmp155 + tmp154;
        assert(tmp156 == tmp156);
      const auto tmp157 = tmp156 + tmp153;
        assert(tmp157 == tmp157);
      const auto tmp158 = tmp157 * tmp66;
        assert(tmp158 == tmp158);
      const auto tmp159 = tmp158 + tmp150;
        assert(tmp159 == tmp159);
      const auto tmp160 = tmp159 / tmp67;
        assert(tmp160 == tmp160);
      const auto tmp161 = tmp42[ 1 ] * tmp160;
        assert(tmp161 == tmp161);
      const auto tmp162 = tmp161 + tmp142;
        assert(tmp162 == tmp162);
      const auto tmp163 = -1 * tmp162;
        assert(tmp163 == tmp163);
      const auto tmp164 = tmp163 * tmp0;
        assert(tmp164 == tmp164);
      const auto tmp165 = tmp164 + tmp79;
        assert(tmp165 == tmp165);
      const auto tmp166 = -1 * tmp78;
        assert(tmp166 == tmp166);
      const auto tmp167 = tmp166 * tmp0;
        assert(tmp167 == tmp167);
      const auto tmp168 = -1 * tmp163;
        assert(tmp168 == tmp168);
      const auto tmp169 = tmp168 * tmp0;
        assert(tmp169 == tmp169);
      const auto tmp170 = tmp169 + tmp167;
        assert(tmp170 == tmp170);
      return std::make_pair( RangeValueType{ { tmp165 }, { { 0, 0 } } }, RangeValueType{ { tmp170 }, { { 0, 0 } } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      using std::max;
      using std::pow;
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      auto tmp1 = facetArea_;
        assert(tmp1 == tmp1);
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
        assert(tmp2 == tmp2);
      auto tmp3 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 + tmp2;
        assert(tmp4 == tmp4);
      const auto tmp5 = 0.5 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 / tmp1;
        assert(tmp6 == tmp6);
      GlobalCoordinateType tmp7 = geometry().global( Dune::Fem::coordinate( xIn ) );
        assert(tmp7 == tmp7);
      const auto tmp8 = -0.45 + tmp7[ 0 ];
        assert(tmp8 == tmp8);
      const auto tmp9 = std::abs( tmp8 );
        assert(tmp9 == tmp9);
      const auto tmp10 = -0.49 + tmp7[ 1 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = std::abs( tmp10 );
        assert(tmp11 == tmp11);
      const auto tmp12 = (tmp11 < 0.03 ? 1.0 : 0.0) * (tmp9 < 0.11 ? 1.0 : 0.0);
        assert(tmp12 == tmp12);
      const auto tmp13 = -0.01999999999999999 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = 0.12 + tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = -1 * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = 1.0 + tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = 0.5 / tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = std::pow( tmp17, 2 );
        assert(tmp18 == tmp18);
      const auto tmp19 = -0.7000000000000002 * tmp12;
        assert(tmp19 == tmp19);
      const auto tmp20 = 2.7 + tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = 2.0 + tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 / tmp20;
        assert(tmp22 == tmp22);
      const auto tmp23 = 0.5 + tmp15;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 / tmp16;
        assert(tmp24 == tmp24);
      const auto tmp25 = std::pow( tmp24, tmp22 );
        assert(tmp25 == tmp25);
      const auto tmp26 = -1 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = 1.0 + tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 * tmp18;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 / 0.0009;
        assert(tmp29 == tmp29);
      const auto tmp30 = -1.0 / tmp16;
        assert(tmp30 == tmp30);
      const auto tmp31 = -1.0 / tmp20;
        assert(tmp31 == tmp31);
      const auto tmp32 = -1.0 + tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = std::pow( tmp24, tmp32 );
        assert(tmp33 == tmp33);
      const auto tmp34 = 4245.0 * tmp12;
        assert(tmp34 == tmp34);
      const auto tmp35 = 755.0 + tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 * tmp31;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 * tmp33;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp37 * tmp30;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 * tmp29;
        assert(tmp39 == tmp39);
      const auto tmp40 = std::max( tmp39, tmp39 );
        assert(tmp40 == tmp40);
      double tmp41 = constant< 1 >();
        assert(tmp41 == tmp41);
      GlobalCoordinateType tmp42 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp42 == tmp42);
      const auto tmp43 = -1 * tmp42[ 1 ];
        assert(tmp43 == tmp43);
      const auto tmp44 = -9.993360000000001e-11 * tmp12;
        assert(tmp44 == tmp44);
      const auto tmp45 = 1e-10 + tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 * tmp43;
        assert(tmp46 == tmp46);
      const auto tmp47 = -1 * tmp42[ 0 ];
        assert(tmp47 == tmp47);
      const auto tmp48 = 5e-11 * tmp12;
        assert(tmp48 == tmp48);
      const auto tmp49 = -5e-11 + tmp48;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp49 * tmp47;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp50 + tmp46;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 * tmp43;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp49 * tmp43;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp45 * tmp47;
        assert(tmp54 == tmp54);
      const auto tmp55 = tmp54 + tmp53;
        assert(tmp55 == tmp55);
      const auto tmp56 = tmp55 * tmp47;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp56 + tmp52;
        assert(tmp57 == tmp57);
      const auto tmp58 = tmp42[ 1 ] * tmp49;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp42[ 0 ] * tmp45;
        assert(tmp59 == tmp59);
      const auto tmp60 = tmp59 + tmp58;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp42[ 0 ] * tmp60;
        assert(tmp61 == tmp61);
      const auto tmp62 = tmp42[ 0 ] * tmp49;
        assert(tmp62 == tmp62);
      const auto tmp63 = tmp42[ 1 ] * tmp45;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp63 + tmp62;
        assert(tmp64 == tmp64);
      const auto tmp65 = tmp42[ 1 ] * tmp64;
        assert(tmp65 == tmp65);
      const auto tmp66 = tmp65 + tmp61;
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp66 + tmp57;
        assert(tmp67 == tmp67);
      const auto tmp68 = 0.5 * tmp67;
        assert(tmp68 == tmp68);
      const auto tmp69 = tmp66 * tmp57;
        assert(tmp69 == tmp69);
      const auto tmp70 = tmp69 / tmp68;
        assert(tmp70 == tmp70);
      const auto tmp71 = tmp70 * tmp41;
        assert(tmp71 == tmp71);
      const auto tmp72 = tmp71 * tmp40;
        assert(tmp72 == tmp72);
      const auto tmp73 = tmp72 / tmp6;
        assert(tmp73 == tmp73);
      const auto tmp74 = tmp73 * tmp0;
        assert(tmp74 == tmp74);
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp75 = evaluateCoefficient< 2, Side::in >( xIn );
        assert(tmp75 == tmp75);
      const auto tmp76 = -1 * tmp75[ 1 ];
        assert(tmp76 == tmp76);
      const auto tmp77 = 1 + tmp76;
        assert(tmp77 == tmp77);
      const auto tmp78 = tmp77 + tmp15;
        assert(tmp78 == tmp78);
      const auto tmp79 = tmp78 / tmp16;
        assert(tmp79 == tmp79);
      const auto tmp80 = std::pow( tmp79, tmp32 );
        assert(tmp80 == tmp80);
      const auto tmp81 = tmp36 * tmp80;
        assert(tmp81 == tmp81);
      const auto tmp82 = tmp81 * tmp30;
        assert(tmp82 == tmp82);
      const auto tmp83 = tmp75[ 1 ] / tmp16;
        assert(tmp83 == tmp83);
      const auto tmp84 = std::pow( tmp83, 2 );
        assert(tmp84 == tmp84);
      const auto tmp85 = std::pow( tmp79, tmp22 );
        assert(tmp85 == tmp85);
      const auto tmp86 = -1 * tmp85;
        assert(tmp86 == tmp86);
      const auto tmp87 = 1.0 + tmp86;
        assert(tmp87 == tmp87);
      const auto tmp88 = tmp87 * tmp84;
        assert(tmp88 == tmp88);
      const auto tmp89 = tmp88 / 0.0009;
        assert(tmp89 == tmp89);
      const auto tmp90 = tmp49 * tmp89;
        assert(tmp90 == tmp90);
      const auto tmp91 = tmp90 * tmp82;
        assert(tmp91 == tmp91);
      const auto tmp92 = tmp57 * tmp91;
        assert(tmp92 == tmp92);
      const auto tmp93 = tmp92 / tmp67;
        assert(tmp93 == tmp93);
      const auto tmp94 = tmp42[ 1 ] * tmp93;
        assert(tmp94 == tmp94);
      const auto tmp95 = tmp45 * tmp89;
        assert(tmp95 == tmp95);
      const auto tmp96 = tmp95 * tmp82;
        assert(tmp96 == tmp96);
      const auto tmp97 = tmp57 * tmp96;
        assert(tmp97 == tmp97);
      const auto tmp98 = tmp97 / tmp67;
        assert(tmp98 == tmp98);
      const auto tmp99 = tmp42[ 0 ] * tmp98;
        assert(tmp99 == tmp99);
      const auto tmp100 = tmp99 + tmp94;
        assert(tmp100 == tmp100);
      const auto tmp101 = -1 * tmp100;
        assert(tmp101 == tmp101);
      const auto tmp102 = tmp101 * tmp0;
        assert(tmp102 == tmp102);
      const auto tmp103 = tmp42[ 0 ] * tmp93;
        assert(tmp103 == tmp103);
      const auto tmp104 = tmp42[ 1 ] * tmp98;
        assert(tmp104 == tmp104);
      const auto tmp105 = tmp104 + tmp103;
        assert(tmp105 == tmp105);
      const auto tmp106 = -1 * tmp105;
        assert(tmp106 == tmp106);
      const auto tmp107 = tmp106 * tmp0;
        assert(tmp107 == tmp107);
      const auto tmp108 = -1 * tmp73;
        assert(tmp108 == tmp108);
      const auto tmp109 = tmp108 * tmp0;
        assert(tmp109 == tmp109);
      const auto tmp110 = -1 * tmp101;
        assert(tmp110 == tmp110);
      const auto tmp111 = tmp110 * tmp0;
        assert(tmp111 == tmp111);
      const auto tmp112 = -1 * tmp106;
        assert(tmp112 == tmp112);
      const auto tmp113 = tmp112 * tmp0;
        assert(tmp113 == tmp113);
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp114 = evaluateCoefficient< 2, Side::out >( xOut );
        assert(tmp114 == tmp114);
      const auto tmp115 = -1 * tmp114[ 1 ];
        assert(tmp115 == tmp115);
      const auto tmp116 = 1 + tmp115;
        assert(tmp116 == tmp116);
      const auto tmp117 = tmp116 + tmp15;
        assert(tmp117 == tmp117);
      const auto tmp118 = tmp117 / tmp16;
        assert(tmp118 == tmp118);
      const auto tmp119 = std::pow( tmp118, tmp32 );
        assert(tmp119 == tmp119);
      const auto tmp120 = tmp36 * tmp119;
        assert(tmp120 == tmp120);
      const auto tmp121 = tmp120 * tmp30;
        assert(tmp121 == tmp121);
      const auto tmp122 = tmp114[ 1 ] / tmp16;
        assert(tmp122 == tmp122);
      const auto tmp123 = std::pow( tmp122, 2 );
        assert(tmp123 == tmp123);
      const auto tmp124 = std::pow( tmp118, tmp22 );
        assert(tmp124 == tmp124);
      const auto tmp125 = -1 * tmp124;
        assert(tmp125 == tmp125);
      const auto tmp126 = 1.0 + tmp125;
        assert(tmp126 == tmp126);
      const auto tmp127 = tmp126 * tmp123;
        assert(tmp127 == tmp127);
      const auto tmp128 = tmp127 / 0.0009;
        assert(tmp128 == tmp128);
      const auto tmp129 = tmp49 * tmp128;
        assert(tmp129 == tmp129);
      const auto tmp130 = tmp129 * tmp121;
        assert(tmp130 == tmp130);
      const auto tmp131 = tmp66 * tmp130;
        assert(tmp131 == tmp131);
      const auto tmp132 = tmp131 / tmp67;
        assert(tmp132 == tmp132);
      const auto tmp133 = tmp42[ 1 ] * tmp132;
        assert(tmp133 == tmp133);
      const auto tmp134 = tmp45 * tmp128;
        assert(tmp134 == tmp134);
      const auto tmp135 = tmp134 * tmp121;
        assert(tmp135 == tmp135);
      const auto tmp136 = tmp66 * tmp135;
        assert(tmp136 == tmp136);
      const auto tmp137 = tmp136 / tmp67;
        assert(tmp137 == tmp137);
      const auto tmp138 = tmp42[ 0 ] * tmp137;
        assert(tmp138 == tmp138);
      const auto tmp139 = tmp138 + tmp133;
        assert(tmp139 == tmp139);
      const auto tmp140 = -1 * tmp139;
        assert(tmp140 == tmp140);
      const auto tmp141 = tmp140 * tmp0;
        assert(tmp141 == tmp141);
      const auto tmp142 = tmp42[ 0 ] * tmp132;
        assert(tmp142 == tmp142);
      const auto tmp143 = tmp42[ 1 ] * tmp137;
        assert(tmp143 == tmp143);
      const auto tmp144 = tmp143 + tmp142;
        assert(tmp144 == tmp144);
      const auto tmp145 = -1 * tmp144;
        assert(tmp145 == tmp145);
      const auto tmp146 = tmp145 * tmp0;
        assert(tmp146 == tmp146);
      const auto tmp147 = -1 * tmp108;
        assert(tmp147 == tmp147);
      const auto tmp148 = tmp147 * tmp0;
        assert(tmp148 == tmp148);
      const auto tmp149 = -1 * tmp140;
        assert(tmp149 == tmp149);
      const auto tmp150 = tmp149 * tmp0;
        assert(tmp150 == tmp150);
      const auto tmp151 = -1 * tmp145;
        assert(tmp151 == tmp151);
      const auto tmp152 = tmp151 * tmp0;
        assert(tmp152 == tmp152);
      return std::make_pair( [ tmp102, tmp107, tmp109, tmp111, tmp113, tmp74 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { (tmp74 * (std::get< 0 >( phiIn ))[ 0 ] + tmp102 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp107 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ] }, { { 0, 0 } } }, RangeValueType{ { (tmp109 * (std::get< 0 >( phiIn ))[ 0 ] + tmp111 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp113 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ] }, { { 0, 0 } } } );
        }, [ tmp109, tmp141, tmp146, tmp148, tmp150, tmp152 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { (tmp109 * (std::get< 0 >( phiOut ))[ 0 ] + tmp141 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp146 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ] }, { { 0, 0 } } }, RangeValueType{ { (tmp148 * (std::get< 0 >( phiOut ))[ 0 ] + tmp150 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp152 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ] }, { { 0, 0 } } } );
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

    const ContimeStep &contimeStep () const
    {
      return *std::get< 0 >( constants_ );
    }

    ContimeStep &contimeStep ()
    {
      return *std::get< 0 >( constants_ );
    }

    const Conpenalty &conpenalty () const
    {
      return *std::get< 1 >( constants_ );
    }

    Conpenalty &conpenalty ()
    {
      return *std::get< 1 >( constants_ );
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

    template< std::size_t i, Side side, class Point >
    typename CoefficientFunctionSpaceType< i >::RangeType evaluateCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::RangeType result;
      std::get< i >( coefficients_[ static_cast< std::size_t >( side ) ] ).evaluate( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::RangeType evaluateCoefficient ( const Point &x ) const
    {
      return evaluateCoefficient< i, Side::in >( x );
    }

    template< std::size_t i, Side side, class Point >
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::JacobianRangeType result;
      std::get< i >( coefficients_[ static_cast< std::size_t >( side ) ] ).jacobian( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
    {
      return jacobianCoefficient< i, Side::in >( x );
    }

    template< std::size_t i, Side side, class Point >
    typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::HessianRangeType result;
      std::get< i >( coefficients_[ static_cast< std::size_t >( side ) ] ).hessian( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
    {
      return hessianCoefficient< i, Side::in >( x );
    }
    std::array< EntityType, 2 > entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    std::array< ctype, 2 > cellVolume_;
    ctype facetArea_;
    ConstantTupleType constants_;
    std::array< std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution >, Dune::Fem::ConstLocalFunction< Coeffsolution_old >, Dune::Fem::ConstLocalFunction< Coeffiterate > >, 2 > coefficients_;
  };

} // namespace Integrands_6df37dd8a3c7afefb581db5ae9568645v1_4

PYBIND11_MODULE( integrands_6df37dd8a3c7afefb581db5ae9568645v1_4_220094643888a7089a01c0741f938529, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > > GridPart;
  typedef Integrands_6df37dd8a3c7afefb581db5ae9568645v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_6df37dd8a3c7afefb581db5ae9568645v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_6df37dd8a3c7afefb581db5ae9568645v1_4_220094643888a7089a01c0741f938529.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffsolution, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffsolution_old, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffiterate ) { return new Integrands( coeffsolution, coeffsolution_old, coeffiterate ); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def_property( "timeStep", [] ( Integrands &self ) -> Integrands::ContimeStep { return self.contimeStep(); }, [] ( Integrands &self, const Integrands::ContimeStep &v ) { self.contimeStep() = v; } );
      cls.def_property( "penalty", [] ( Integrands &self ) -> Integrands::Conpenalty { return self.conpenalty(); }, [] ( Integrands &self, const Integrands::Conpenalty &v ) { self.conpenalty() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_6df37dd8a3c7afefb581db5ae9568645v1_4
