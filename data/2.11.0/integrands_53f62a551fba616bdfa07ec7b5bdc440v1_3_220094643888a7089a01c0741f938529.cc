#ifndef GuardIntegrands_53f62a551fba616bdfa07ec7b5bdc440v1_3
#define GuardIntegrands_53f62a551fba616bdfa07ec7b5bdc440v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
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
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_53f62a551fba616bdfa07ec7b5bdc440v1_3
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
      const auto tmp1 = -0.45 + tmp0[ 0 ];
      const auto tmp2 = std::abs( tmp1 );
      const auto tmp3 = -0.49 + tmp0[ 1 ];
      const auto tmp4 = std::abs( tmp3 );
      const auto tmp5 = (tmp4 < 0.03 ? 1.0 : 0.0) * (tmp2 < 0.11 ? 1.0 : 0.0);
      const auto tmp6 = -0.010000000000000009 * tmp5;
      const auto tmp7 = 0.4 + tmp6;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp8 = evaluateCoefficient< 1 >( x );
      const auto tmp9 = -1 * tmp8[ 1 ];
      const auto tmp10 = std::get< 0 >( u );
      const auto tmp11 = tmp10[ 0 ] + tmp9;
      const auto tmp12 = tmp11 * tmp7;
      const auto tmp13 = -0.01999999999999999 * tmp5;
      const auto tmp14 = 0.12 + tmp13;
      const auto tmp15 = -1 * tmp14;
      const auto tmp16 = 1.0 + tmp15;
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp17 = evaluateCoefficient< 2 >( x );
      const auto tmp18 = tmp17[ 1 ] / tmp16;
      const auto tmp19 = std::pow( tmp18, 2 );
      const auto tmp20 = -0.7000000000000002 * tmp5;
      const auto tmp21 = 2.7 + tmp20;
      const auto tmp22 = 2.0 + tmp21;
      const auto tmp23 = tmp22 / tmp21;
      const auto tmp24 = -1 * tmp17[ 1 ];
      const auto tmp25 = 1 + tmp24;
      const auto tmp26 = tmp25 + tmp15;
      const auto tmp27 = tmp26 / tmp16;
      const auto tmp28 = std::pow( tmp27, tmp23 );
      const auto tmp29 = -1 * tmp28;
      const auto tmp30 = 1.0 + tmp29;
      const auto tmp31 = tmp30 * tmp19;
      const auto tmp32 = tmp31 / 0.0009;
      const auto tmp33 = 5e-11 * tmp5;
      const auto tmp34 = -5e-11 + tmp33;
      const auto tmp35 = tmp34 * tmp32;
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp36 = jacobianCoefficient< 0 >( x );
      const auto tmp37 = 14322.6 + (tmp36[ 0 ])[ 1 ];
      const auto tmp38 = tmp37 * tmp35;
      const auto tmp39 = -9.993360000000001e-11 * tmp5;
      const auto tmp40 = 1e-10 + tmp39;
      const auto tmp41 = tmp40 * tmp32;
      const auto tmp42 = (tmp36[ 0 ])[ 0 ] * tmp41;
      const auto tmp43 = tmp42 + tmp38;
      const auto tmp44 = -1.0 / tmp16;
      const auto tmp45 = -1.0 / tmp21;
      const auto tmp46 = -1.0 + tmp45;
      const auto tmp47 = std::pow( tmp27, tmp46 );
      const auto tmp48 = 4245.0 * tmp5;
      const auto tmp49 = 755.0 + tmp48;
      const auto tmp50 = tmp49 * tmp45;
      const auto tmp51 = tmp50 * tmp47;
      const auto tmp52 = tmp51 * tmp44;
      const auto tmp53 = tmp35 * tmp52;
      const auto tmp54 = std::get< 1 >( u );
      const auto tmp55 = (tmp54[ 0 ])[ 1 ] * tmp53;
      const auto tmp56 = tmp41 * tmp52;
      const auto tmp57 = (tmp54[ 0 ])[ 0 ] * tmp56;
      const auto tmp58 = tmp57 + tmp55;
      const auto tmp59 = tmp58 + tmp43;
      double tmp60 = constant< 0 >();
      const auto tmp61 = tmp60 * tmp59;
      const auto tmp62 = (tmp36[ 0 ])[ 0 ] * tmp35;
      const auto tmp63 = tmp37 * tmp41;
      const auto tmp64 = tmp63 + tmp62;
      const auto tmp65 = (tmp54[ 0 ])[ 0 ] * tmp53;
      const auto tmp66 = (tmp54[ 0 ])[ 1 ] * tmp56;
      const auto tmp67 = tmp66 + tmp65;
      const auto tmp68 = tmp67 + tmp64;
      const auto tmp69 = tmp60 * tmp68;
      return RangeValueType{ { tmp12 }, { { tmp61, tmp69 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = -0.45 + tmp0[ 0 ];
      const auto tmp2 = std::abs( tmp1 );
      const auto tmp3 = -0.49 + tmp0[ 1 ];
      const auto tmp4 = std::abs( tmp3 );
      const auto tmp5 = (tmp4 < 0.03 ? 1.0 : 0.0) * (tmp2 < 0.11 ? 1.0 : 0.0);
      const auto tmp6 = -0.010000000000000009 * tmp5;
      const auto tmp7 = 0.4 + tmp6;
      const auto tmp8 = -0.01999999999999999 * tmp5;
      const auto tmp9 = 0.12 + tmp8;
      const auto tmp10 = -1 * tmp9;
      const auto tmp11 = 1.0 + tmp10;
      const auto tmp12 = -1.0 / tmp11;
      const auto tmp13 = -0.7000000000000002 * tmp5;
      const auto tmp14 = 2.7 + tmp13;
      const auto tmp15 = -1.0 / tmp14;
      const auto tmp16 = -1.0 + tmp15;
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp17 = evaluateCoefficient< 2 >( x );
      const auto tmp18 = -1 * tmp17[ 1 ];
      const auto tmp19 = 1 + tmp18;
      const auto tmp20 = tmp19 + tmp10;
      const auto tmp21 = tmp20 / tmp11;
      const auto tmp22 = std::pow( tmp21, tmp16 );
      const auto tmp23 = 4245.0 * tmp5;
      const auto tmp24 = 755.0 + tmp23;
      const auto tmp25 = tmp24 * tmp15;
      const auto tmp26 = tmp25 * tmp22;
      const auto tmp27 = tmp26 * tmp12;
      const auto tmp28 = tmp17[ 1 ] / tmp11;
      const auto tmp29 = std::pow( tmp28, 2 );
      const auto tmp30 = 2.0 + tmp14;
      const auto tmp31 = tmp30 / tmp14;
      const auto tmp32 = std::pow( tmp21, tmp31 );
      const auto tmp33 = -1 * tmp32;
      const auto tmp34 = 1.0 + tmp33;
      const auto tmp35 = tmp34 * tmp29;
      const auto tmp36 = tmp35 / 0.0009;
      const auto tmp37 = -9.993360000000001e-11 * tmp5;
      const auto tmp38 = 1e-10 + tmp37;
      const auto tmp39 = tmp38 * tmp36;
      const auto tmp40 = tmp39 * tmp27;
      double tmp41 = constant< 0 >();
      const auto tmp42 = tmp41 * tmp40;
      const auto tmp43 = 5e-11 * tmp5;
      const auto tmp44 = -5e-11 + tmp43;
      const auto tmp45 = tmp44 * tmp36;
      const auto tmp46 = tmp45 * tmp27;
      const auto tmp47 = tmp41 * tmp46;
      return [ tmp42, tmp47, tmp7 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp7 * (std::get< 0 >( phi ))[ 0 ] }, { { tmp42 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp47 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], tmp47 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp42 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = -0.9 + tmp0[ 0 ];
      const auto tmp2 = std::abs( tmp1 );
      const auto tmp3 = std::abs( tmp0[ 0 ] );
      const auto tmp4 = (tmp3 < 1e-08 ? 1.0 : 0.0) + (tmp2 < 1e-08 ? 1.0 : 0.0);
      const auto tmp5 = -0.45 + tmp0[ 0 ];
      const auto tmp6 = std::abs( tmp5 );
      const auto tmp7 = -0.49 + tmp0[ 1 ];
      const auto tmp8 = std::abs( tmp7 );
      const auto tmp9 = (tmp8 < 0.03 ? 1.0 : 0.0) * (tmp6 < 0.11 ? 1.0 : 0.0);
      const auto tmp10 = -0.01999999999999999 * tmp9;
      const auto tmp11 = 0.12 + tmp10;
      const auto tmp12 = -1 * tmp11;
      const auto tmp13 = 1.0 + tmp12;
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp14 = evaluateCoefficient< 2 >( x );
      const auto tmp15 = tmp14[ 1 ] / tmp13;
      const auto tmp16 = std::pow( tmp15, 2 );
      const auto tmp17 = -0.7000000000000002 * tmp9;
      const auto tmp18 = 2.7 + tmp17;
      const auto tmp19 = 2.0 + tmp18;
      const auto tmp20 = tmp19 / tmp18;
      const auto tmp21 = -1 * tmp14[ 1 ];
      const auto tmp22 = 1 + tmp21;
      const auto tmp23 = tmp22 + tmp12;
      const auto tmp24 = tmp23 / tmp13;
      const auto tmp25 = std::pow( tmp24, tmp20 );
      const auto tmp26 = -1 * tmp25;
      const auto tmp27 = 1.0 + tmp26;
      const auto tmp28 = tmp27 * tmp16;
      const auto tmp29 = tmp28 / 0.0009;
      const auto tmp30 = 5e-11 * tmp9;
      const auto tmp31 = -5e-11 + tmp30;
      const auto tmp32 = tmp31 * tmp29;
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp33 = jacobianCoefficient< 0 >( x );
      const auto tmp34 = 14322.6 + (tmp33[ 0 ])[ 1 ];
      const auto tmp35 = tmp34 * tmp32;
      const auto tmp36 = -9.993360000000001e-11 * tmp9;
      const auto tmp37 = 1e-10 + tmp36;
      const auto tmp38 = tmp37 * tmp29;
      const auto tmp39 = (tmp33[ 0 ])[ 0 ] * tmp38;
      const auto tmp40 = tmp39 + tmp35;
      const auto tmp41 = -1.0 / tmp13;
      const auto tmp42 = -1.0 / tmp18;
      const auto tmp43 = -1.0 + tmp42;
      const auto tmp44 = std::pow( tmp24, tmp43 );
      const auto tmp45 = 4245.0 * tmp9;
      const auto tmp46 = 755.0 + tmp45;
      const auto tmp47 = tmp46 * tmp42;
      const auto tmp48 = tmp47 * tmp44;
      const auto tmp49 = tmp48 * tmp41;
      const auto tmp50 = tmp32 * tmp49;
      const auto tmp51 = std::get< 1 >( u );
      const auto tmp52 = (tmp51[ 0 ])[ 1 ] * tmp50;
      const auto tmp53 = tmp38 * tmp49;
      const auto tmp54 = (tmp51[ 0 ])[ 0 ] * tmp53;
      const auto tmp55 = tmp54 + tmp52;
      const auto tmp56 = tmp55 + tmp40;
      GlobalCoordinateType tmp57 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp58 = tmp57[ 0 ] * tmp56;
      const auto tmp59 = (tmp33[ 0 ])[ 0 ] * tmp32;
      const auto tmp60 = tmp34 * tmp38;
      const auto tmp61 = tmp60 + tmp59;
      const auto tmp62 = (tmp51[ 0 ])[ 0 ] * tmp50;
      const auto tmp63 = (tmp51[ 0 ])[ 1 ] * tmp53;
      const auto tmp64 = tmp63 + tmp62;
      const auto tmp65 = tmp64 + tmp61;
      const auto tmp66 = tmp57[ 1 ] * tmp65;
      const auto tmp67 = tmp66 + tmp58;
      const auto tmp68 = tmp67 * tmp4;
      const auto tmp69 = -1 * tmp68;
      double tmp70 = constant< 0 >();
      const auto tmp71 = tmp70 * tmp69;
      const auto tmp72 = -0.65 + tmp0[ 1 ];
      const auto tmp73 = std::abs( tmp72 );
      const auto tmp74 = (tmp6 < 0.06 ? 1.0 : 0.0) * (tmp73 < 1e-08 ? 1.0 : 0.0);
      const auto tmp75 = -5.137e-05 * tmp74;
      const auto tmp76 = tmp70 * tmp75;
      auto tmp77 = facetArea_;
      auto tmp78 = cellVolume();
      const auto tmp79 = tmp78 / tmp77;
      const auto tmp80 = 0.5 / tmp13;
      const auto tmp81 = std::pow( tmp80, 2 );
      const auto tmp82 = 0.5 + tmp12;
      const auto tmp83 = tmp82 / tmp13;
      const auto tmp84 = std::pow( tmp83, tmp20 );
      const auto tmp85 = -1 * tmp84;
      const auto tmp86 = 1.0 + tmp85;
      const auto tmp87 = tmp86 * tmp81;
      const auto tmp88 = tmp87 / 0.0009;
      const auto tmp89 = std::pow( tmp83, tmp43 );
      const auto tmp90 = tmp47 * tmp89;
      const auto tmp91 = tmp90 * tmp41;
      const auto tmp92 = tmp91 * tmp88;
      const auto tmp93 = tmp57[ 1 ] * tmp31;
      const auto tmp94 = tmp57[ 0 ] * tmp37;
      const auto tmp95 = tmp94 + tmp93;
      const auto tmp96 = tmp57[ 0 ] * tmp95;
      const auto tmp97 = tmp57[ 0 ] * tmp31;
      const auto tmp98 = tmp57[ 1 ] * tmp37;
      const auto tmp99 = tmp98 + tmp97;
      const auto tmp100 = tmp57[ 1 ] * tmp99;
      const auto tmp101 = tmp100 + tmp96;
      double tmp102 = constant< 1 >();
      const auto tmp103 = tmp102 * tmp101;
      const auto tmp104 = tmp103 * tmp92;
      const auto tmp105 = tmp104 / tmp79;
      const auto tmp106 = std::get< 0 >( u );
      const auto tmp107 = tmp106[ 0 ] * tmp105;
      const auto tmp108 = tmp4 * tmp107;
      const auto tmp109 = tmp70 * tmp108;
      const auto tmp110 = tmp109 + tmp76;
      const auto tmp111 = tmp110 + tmp71;
      return RangeValueType{ { tmp111 }, { { 0, 0 } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      auto tmp0 = facetArea_;
      auto tmp1 = cellVolume();
      const auto tmp2 = tmp1 / tmp0;
      GlobalCoordinateType tmp3 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp4 = -0.45 + tmp3[ 0 ];
      const auto tmp5 = std::abs( tmp4 );
      const auto tmp6 = -0.49 + tmp3[ 1 ];
      const auto tmp7 = std::abs( tmp6 );
      const auto tmp8 = (tmp7 < 0.03 ? 1.0 : 0.0) * (tmp5 < 0.11 ? 1.0 : 0.0);
      const auto tmp9 = -0.01999999999999999 * tmp8;
      const auto tmp10 = 0.12 + tmp9;
      const auto tmp11 = -1 * tmp10;
      const auto tmp12 = 1.0 + tmp11;
      const auto tmp13 = 0.5 / tmp12;
      const auto tmp14 = std::pow( tmp13, 2 );
      const auto tmp15 = -0.7000000000000002 * tmp8;
      const auto tmp16 = 2.7 + tmp15;
      const auto tmp17 = 2.0 + tmp16;
      const auto tmp18 = tmp17 / tmp16;
      const auto tmp19 = 0.5 + tmp11;
      const auto tmp20 = tmp19 / tmp12;
      const auto tmp21 = std::pow( tmp20, tmp18 );
      const auto tmp22 = -1 * tmp21;
      const auto tmp23 = 1.0 + tmp22;
      const auto tmp24 = tmp23 * tmp14;
      const auto tmp25 = tmp24 / 0.0009;
      const auto tmp26 = -1.0 / tmp12;
      const auto tmp27 = -1.0 / tmp16;
      const auto tmp28 = -1.0 + tmp27;
      const auto tmp29 = std::pow( tmp20, tmp28 );
      const auto tmp30 = 4245.0 * tmp8;
      const auto tmp31 = 755.0 + tmp30;
      const auto tmp32 = tmp31 * tmp27;
      const auto tmp33 = tmp32 * tmp29;
      const auto tmp34 = tmp33 * tmp26;
      const auto tmp35 = tmp34 * tmp25;
      const auto tmp36 = 5e-11 * tmp8;
      const auto tmp37 = -5e-11 + tmp36;
      GlobalCoordinateType tmp38 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp39 = tmp38[ 1 ] * tmp37;
      const auto tmp40 = -9.993360000000001e-11 * tmp8;
      const auto tmp41 = 1e-10 + tmp40;
      const auto tmp42 = tmp38[ 0 ] * tmp41;
      const auto tmp43 = tmp42 + tmp39;
      const auto tmp44 = tmp38[ 0 ] * tmp43;
      const auto tmp45 = tmp38[ 0 ] * tmp37;
      const auto tmp46 = tmp38[ 1 ] * tmp41;
      const auto tmp47 = tmp46 + tmp45;
      const auto tmp48 = tmp38[ 1 ] * tmp47;
      const auto tmp49 = tmp48 + tmp44;
      double tmp50 = constant< 1 >();
      const auto tmp51 = tmp50 * tmp49;
      const auto tmp52 = tmp51 * tmp35;
      const auto tmp53 = tmp52 / tmp2;
      const auto tmp54 = -0.9 + tmp3[ 0 ];
      const auto tmp55 = std::abs( tmp54 );
      const auto tmp56 = std::abs( tmp3[ 0 ] );
      const auto tmp57 = (tmp56 < 1e-08 ? 1.0 : 0.0) + (tmp55 < 1e-08 ? 1.0 : 0.0);
      const auto tmp58 = tmp57 * tmp53;
      double tmp59 = constant< 0 >();
      const auto tmp60 = tmp59 * tmp58;
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp61 = evaluateCoefficient< 2 >( x );
      const auto tmp62 = -1 * tmp61[ 1 ];
      const auto tmp63 = 1 + tmp62;
      const auto tmp64 = tmp63 + tmp11;
      const auto tmp65 = tmp64 / tmp12;
      const auto tmp66 = std::pow( tmp65, tmp28 );
      const auto tmp67 = tmp32 * tmp66;
      const auto tmp68 = tmp67 * tmp26;
      const auto tmp69 = tmp61[ 1 ] / tmp12;
      const auto tmp70 = std::pow( tmp69, 2 );
      const auto tmp71 = std::pow( tmp65, tmp18 );
      const auto tmp72 = -1 * tmp71;
      const auto tmp73 = 1.0 + tmp72;
      const auto tmp74 = tmp73 * tmp70;
      const auto tmp75 = tmp74 / 0.0009;
      const auto tmp76 = tmp37 * tmp75;
      const auto tmp77 = tmp76 * tmp68;
      const auto tmp78 = tmp38[ 1 ] * tmp77;
      const auto tmp79 = tmp41 * tmp75;
      const auto tmp80 = tmp79 * tmp68;
      const auto tmp81 = tmp38[ 0 ] * tmp80;
      const auto tmp82 = tmp81 + tmp78;
      const auto tmp83 = tmp82 * tmp57;
      const auto tmp84 = -1 * tmp83;
      const auto tmp85 = tmp59 * tmp84;
      const auto tmp86 = tmp38[ 0 ] * tmp77;
      const auto tmp87 = tmp38[ 1 ] * tmp80;
      const auto tmp88 = tmp87 + tmp86;
      const auto tmp89 = tmp88 * tmp57;
      const auto tmp90 = -1 * tmp89;
      const auto tmp91 = tmp59 * tmp90;
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
      auto tmp0 = facetArea_;
      auto tmp1 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
      const auto tmp3 = tmp2 + tmp1;
      const auto tmp4 = 0.5 * tmp3;
      const auto tmp5 = tmp4 / tmp0;
      GlobalCoordinateType tmp6 = geometry().global( Dune::Fem::coordinate( xIn ) );
      const auto tmp7 = -0.45 + tmp6[ 0 ];
      const auto tmp8 = std::abs( tmp7 );
      const auto tmp9 = -0.49 + tmp6[ 1 ];
      const auto tmp10 = std::abs( tmp9 );
      const auto tmp11 = (tmp10 < 0.03 ? 1.0 : 0.0) * (tmp8 < 0.11 ? 1.0 : 0.0);
      const auto tmp12 = -0.01999999999999999 * tmp11;
      const auto tmp13 = 0.12 + tmp12;
      const auto tmp14 = -1 * tmp13;
      const auto tmp15 = 1.0 + tmp14;
      const auto tmp16 = 0.5 / tmp15;
      const auto tmp17 = std::pow( tmp16, 2 );
      const auto tmp18 = -0.7000000000000002 * tmp11;
      const auto tmp19 = 2.7 + tmp18;
      const auto tmp20 = 2.0 + tmp19;
      const auto tmp21 = tmp20 / tmp19;
      const auto tmp22 = 0.5 + tmp14;
      const auto tmp23 = tmp22 / tmp15;
      const auto tmp24 = std::pow( tmp23, tmp21 );
      const auto tmp25 = -1 * tmp24;
      const auto tmp26 = 1.0 + tmp25;
      const auto tmp27 = tmp26 * tmp17;
      const auto tmp28 = tmp27 / 0.0009;
      const auto tmp29 = -1.0 / tmp15;
      const auto tmp30 = -1.0 / tmp19;
      const auto tmp31 = -1.0 + tmp30;
      const auto tmp32 = std::pow( tmp23, tmp31 );
      const auto tmp33 = 4245.0 * tmp11;
      const auto tmp34 = 755.0 + tmp33;
      const auto tmp35 = tmp34 * tmp30;
      const auto tmp36 = tmp35 * tmp32;
      const auto tmp37 = tmp36 * tmp29;
      const auto tmp38 = tmp37 * tmp28;
      const auto tmp39 = std::max( tmp38, tmp38 );
      GlobalCoordinateType tmp40 = intersection_.unitOuterNormal( xIn.localPosition() );
      const auto tmp41 = -1 * tmp40[ 1 ];
      const auto tmp42 = -9.993360000000001e-11 * tmp11;
      const auto tmp43 = 1e-10 + tmp42;
      const auto tmp44 = tmp43 * tmp41;
      const auto tmp45 = -1 * tmp40[ 0 ];
      const auto tmp46 = 5e-11 * tmp11;
      const auto tmp47 = -5e-11 + tmp46;
      const auto tmp48 = tmp47 * tmp45;
      const auto tmp49 = tmp48 + tmp44;
      const auto tmp50 = tmp49 * tmp41;
      const auto tmp51 = tmp47 * tmp41;
      const auto tmp52 = tmp43 * tmp45;
      const auto tmp53 = tmp52 + tmp51;
      const auto tmp54 = tmp53 * tmp45;
      const auto tmp55 = tmp54 + tmp50;
      const auto tmp56 = tmp40[ 1 ] * tmp47;
      const auto tmp57 = tmp40[ 0 ] * tmp43;
      const auto tmp58 = tmp57 + tmp56;
      const auto tmp59 = tmp40[ 0 ] * tmp58;
      const auto tmp60 = tmp40[ 0 ] * tmp47;
      const auto tmp61 = tmp40[ 1 ] * tmp43;
      const auto tmp62 = tmp61 + tmp60;
      const auto tmp63 = tmp40[ 1 ] * tmp62;
      const auto tmp64 = tmp63 + tmp59;
      const auto tmp65 = tmp64 + tmp55;
      const auto tmp66 = 0.5 * tmp65;
      const auto tmp67 = tmp64 * tmp55;
      const auto tmp68 = tmp67 / tmp66;
      double tmp69 = constant< 1 >();
      const auto tmp70 = tmp69 * tmp68;
      const auto tmp71 = tmp70 * tmp39;
      const auto tmp72 = tmp71 / tmp5;
      const auto tmp73 = std::get< 0 >( uOut );
      const auto tmp74 = -1 * tmp73[ 0 ];
      const auto tmp75 = std::get< 0 >( uIn );
      const auto tmp76 = tmp75[ 0 ] + tmp74;
      const auto tmp77 = tmp76 * tmp72;
      double tmp78 = constant< 0 >();
      const auto tmp79 = tmp78 * tmp77;
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp80 = evaluateCoefficient< 2, Side::in >( xIn );
      const auto tmp81 = tmp80[ 1 ] / tmp15;
      const auto tmp82 = std::pow( tmp81, 2 );
      const auto tmp83 = -1 * tmp80[ 1 ];
      const auto tmp84 = 1 + tmp83;
      const auto tmp85 = tmp84 + tmp14;
      const auto tmp86 = tmp85 / tmp15;
      const auto tmp87 = std::pow( tmp86, tmp21 );
      const auto tmp88 = -1 * tmp87;
      const auto tmp89 = 1.0 + tmp88;
      const auto tmp90 = tmp89 * tmp82;
      const auto tmp91 = tmp90 / 0.0009;
      const auto tmp92 = tmp47 * tmp91;
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp93 = jacobianCoefficient< 0, Side::in >( xIn );
      const auto tmp94 = 14322.6 + (tmp93[ 0 ])[ 1 ];
      const auto tmp95 = tmp94 * tmp92;
      const auto tmp96 = tmp43 * tmp91;
      const auto tmp97 = (tmp93[ 0 ])[ 0 ] * tmp96;
      const auto tmp98 = tmp97 + tmp95;
      const auto tmp99 = std::pow( tmp86, tmp31 );
      const auto tmp100 = tmp35 * tmp99;
      const auto tmp101 = tmp100 * tmp29;
      const auto tmp102 = tmp92 * tmp101;
      const auto tmp103 = std::get< 1 >( uIn );
      const auto tmp104 = (tmp103[ 0 ])[ 1 ] * tmp102;
      const auto tmp105 = tmp96 * tmp101;
      const auto tmp106 = (tmp103[ 0 ])[ 0 ] * tmp105;
      const auto tmp107 = tmp106 + tmp104;
      const auto tmp108 = tmp107 + tmp98;
      const auto tmp109 = tmp108 * tmp55;
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp110 = evaluateCoefficient< 2, Side::out >( xOut );
      const auto tmp111 = tmp110[ 1 ] / tmp15;
      const auto tmp112 = std::pow( tmp111, 2 );
      const auto tmp113 = -1 * tmp110[ 1 ];
      const auto tmp114 = 1 + tmp113;
      const auto tmp115 = tmp114 + tmp14;
      const auto tmp116 = tmp115 / tmp15;
      const auto tmp117 = std::pow( tmp116, tmp21 );
      const auto tmp118 = -1 * tmp117;
      const auto tmp119 = 1.0 + tmp118;
      const auto tmp120 = tmp119 * tmp112;
      const auto tmp121 = tmp120 / 0.0009;
      const auto tmp122 = tmp47 * tmp121;
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp123 = jacobianCoefficient< 0, Side::out >( xOut );
      const auto tmp124 = 14322.6 + (tmp123[ 0 ])[ 1 ];
      const auto tmp125 = tmp124 * tmp122;
      const auto tmp126 = tmp43 * tmp121;
      const auto tmp127 = (tmp123[ 0 ])[ 0 ] * tmp126;
      const auto tmp128 = tmp127 + tmp125;
      const auto tmp129 = std::pow( tmp116, tmp31 );
      const auto tmp130 = tmp35 * tmp129;
      const auto tmp131 = tmp130 * tmp29;
      const auto tmp132 = tmp122 * tmp131;
      const auto tmp133 = std::get< 1 >( uOut );
      const auto tmp134 = (tmp133[ 0 ])[ 1 ] * tmp132;
      const auto tmp135 = tmp126 * tmp131;
      const auto tmp136 = (tmp133[ 0 ])[ 0 ] * tmp135;
      const auto tmp137 = tmp136 + tmp134;
      const auto tmp138 = tmp137 + tmp128;
      const auto tmp139 = tmp138 * tmp64;
      const auto tmp140 = tmp139 + tmp109;
      const auto tmp141 = tmp140 / tmp65;
      const auto tmp142 = tmp40[ 0 ] * tmp141;
      const auto tmp143 = (tmp93[ 0 ])[ 0 ] * tmp92;
      const auto tmp144 = tmp94 * tmp96;
      const auto tmp145 = tmp144 + tmp143;
      const auto tmp146 = (tmp103[ 0 ])[ 0 ] * tmp102;
      const auto tmp147 = (tmp103[ 0 ])[ 1 ] * tmp105;
      const auto tmp148 = tmp147 + tmp146;
      const auto tmp149 = tmp148 + tmp145;
      const auto tmp150 = tmp149 * tmp55;
      const auto tmp151 = (tmp123[ 0 ])[ 0 ] * tmp122;
      const auto tmp152 = tmp124 * tmp126;
      const auto tmp153 = tmp152 + tmp151;
      const auto tmp154 = (tmp133[ 0 ])[ 0 ] * tmp132;
      const auto tmp155 = (tmp133[ 0 ])[ 1 ] * tmp135;
      const auto tmp156 = tmp155 + tmp154;
      const auto tmp157 = tmp156 + tmp153;
      const auto tmp158 = tmp157 * tmp64;
      const auto tmp159 = tmp158 + tmp150;
      const auto tmp160 = tmp159 / tmp65;
      const auto tmp161 = tmp40[ 1 ] * tmp160;
      const auto tmp162 = tmp161 + tmp142;
      const auto tmp163 = -1 * tmp162;
      const auto tmp164 = tmp78 * tmp163;
      const auto tmp165 = tmp164 + tmp79;
      const auto tmp166 = -1 * tmp77;
      const auto tmp167 = tmp78 * tmp166;
      const auto tmp168 = -1 * tmp163;
      const auto tmp169 = tmp78 * tmp168;
      const auto tmp170 = tmp169 + tmp167;
      return std::make_pair( RangeValueType{ { tmp165 }, { { 0, 0 } } }, RangeValueType{ { tmp170 }, { { 0, 0 } } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      using std::max;
      using std::pow;
      auto tmp0 = facetArea_;
      auto tmp1 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
      const auto tmp3 = tmp2 + tmp1;
      const auto tmp4 = 0.5 * tmp3;
      const auto tmp5 = tmp4 / tmp0;
      GlobalCoordinateType tmp6 = geometry().global( Dune::Fem::coordinate( xIn ) );
      const auto tmp7 = -0.45 + tmp6[ 0 ];
      const auto tmp8 = std::abs( tmp7 );
      const auto tmp9 = -0.49 + tmp6[ 1 ];
      const auto tmp10 = std::abs( tmp9 );
      const auto tmp11 = (tmp10 < 0.03 ? 1.0 : 0.0) * (tmp8 < 0.11 ? 1.0 : 0.0);
      const auto tmp12 = -0.01999999999999999 * tmp11;
      const auto tmp13 = 0.12 + tmp12;
      const auto tmp14 = -1 * tmp13;
      const auto tmp15 = 1.0 + tmp14;
      const auto tmp16 = 0.5 / tmp15;
      const auto tmp17 = std::pow( tmp16, 2 );
      const auto tmp18 = -0.7000000000000002 * tmp11;
      const auto tmp19 = 2.7 + tmp18;
      const auto tmp20 = 2.0 + tmp19;
      const auto tmp21 = tmp20 / tmp19;
      const auto tmp22 = 0.5 + tmp14;
      const auto tmp23 = tmp22 / tmp15;
      const auto tmp24 = std::pow( tmp23, tmp21 );
      const auto tmp25 = -1 * tmp24;
      const auto tmp26 = 1.0 + tmp25;
      const auto tmp27 = tmp26 * tmp17;
      const auto tmp28 = tmp27 / 0.0009;
      const auto tmp29 = -1.0 / tmp15;
      const auto tmp30 = -1.0 / tmp19;
      const auto tmp31 = -1.0 + tmp30;
      const auto tmp32 = std::pow( tmp23, tmp31 );
      const auto tmp33 = 4245.0 * tmp11;
      const auto tmp34 = 755.0 + tmp33;
      const auto tmp35 = tmp34 * tmp30;
      const auto tmp36 = tmp35 * tmp32;
      const auto tmp37 = tmp36 * tmp29;
      const auto tmp38 = tmp37 * tmp28;
      const auto tmp39 = std::max( tmp38, tmp38 );
      GlobalCoordinateType tmp40 = intersection_.unitOuterNormal( xIn.localPosition() );
      const auto tmp41 = -1 * tmp40[ 1 ];
      const auto tmp42 = -9.993360000000001e-11 * tmp11;
      const auto tmp43 = 1e-10 + tmp42;
      const auto tmp44 = tmp43 * tmp41;
      const auto tmp45 = -1 * tmp40[ 0 ];
      const auto tmp46 = 5e-11 * tmp11;
      const auto tmp47 = -5e-11 + tmp46;
      const auto tmp48 = tmp47 * tmp45;
      const auto tmp49 = tmp48 + tmp44;
      const auto tmp50 = tmp49 * tmp41;
      const auto tmp51 = tmp47 * tmp41;
      const auto tmp52 = tmp43 * tmp45;
      const auto tmp53 = tmp52 + tmp51;
      const auto tmp54 = tmp53 * tmp45;
      const auto tmp55 = tmp54 + tmp50;
      const auto tmp56 = tmp40[ 1 ] * tmp47;
      const auto tmp57 = tmp40[ 0 ] * tmp43;
      const auto tmp58 = tmp57 + tmp56;
      const auto tmp59 = tmp40[ 0 ] * tmp58;
      const auto tmp60 = tmp40[ 0 ] * tmp47;
      const auto tmp61 = tmp40[ 1 ] * tmp43;
      const auto tmp62 = tmp61 + tmp60;
      const auto tmp63 = tmp40[ 1 ] * tmp62;
      const auto tmp64 = tmp63 + tmp59;
      const auto tmp65 = tmp64 + tmp55;
      const auto tmp66 = 0.5 * tmp65;
      const auto tmp67 = tmp64 * tmp55;
      const auto tmp68 = tmp67 / tmp66;
      double tmp69 = constant< 1 >();
      const auto tmp70 = tmp69 * tmp68;
      const auto tmp71 = tmp70 * tmp39;
      const auto tmp72 = tmp71 / tmp5;
      double tmp73 = constant< 0 >();
      const auto tmp74 = tmp73 * tmp72;
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp75 = evaluateCoefficient< 2, Side::in >( xIn );
      const auto tmp76 = -1 * tmp75[ 1 ];
      const auto tmp77 = 1 + tmp76;
      const auto tmp78 = tmp77 + tmp14;
      const auto tmp79 = tmp78 / tmp15;
      const auto tmp80 = std::pow( tmp79, tmp31 );
      const auto tmp81 = tmp35 * tmp80;
      const auto tmp82 = tmp81 * tmp29;
      const auto tmp83 = tmp75[ 1 ] / tmp15;
      const auto tmp84 = std::pow( tmp83, 2 );
      const auto tmp85 = std::pow( tmp79, tmp21 );
      const auto tmp86 = -1 * tmp85;
      const auto tmp87 = 1.0 + tmp86;
      const auto tmp88 = tmp87 * tmp84;
      const auto tmp89 = tmp88 / 0.0009;
      const auto tmp90 = tmp47 * tmp89;
      const auto tmp91 = tmp90 * tmp82;
      const auto tmp92 = tmp55 * tmp91;
      const auto tmp93 = tmp92 / tmp65;
      const auto tmp94 = tmp40[ 1 ] * tmp93;
      const auto tmp95 = tmp43 * tmp89;
      const auto tmp96 = tmp95 * tmp82;
      const auto tmp97 = tmp55 * tmp96;
      const auto tmp98 = tmp97 / tmp65;
      const auto tmp99 = tmp40[ 0 ] * tmp98;
      const auto tmp100 = tmp99 + tmp94;
      const auto tmp101 = -1 * tmp100;
      const auto tmp102 = tmp73 * tmp101;
      const auto tmp103 = tmp40[ 0 ] * tmp93;
      const auto tmp104 = tmp40[ 1 ] * tmp98;
      const auto tmp105 = tmp104 + tmp103;
      const auto tmp106 = -1 * tmp105;
      const auto tmp107 = tmp73 * tmp106;
      const auto tmp108 = -1 * tmp72;
      const auto tmp109 = tmp73 * tmp108;
      const auto tmp110 = -1 * tmp101;
      const auto tmp111 = tmp73 * tmp110;
      const auto tmp112 = -1 * tmp106;
      const auto tmp113 = tmp73 * tmp112;
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp114 = evaluateCoefficient< 2, Side::out >( xOut );
      const auto tmp115 = -1 * tmp114[ 1 ];
      const auto tmp116 = 1 + tmp115;
      const auto tmp117 = tmp116 + tmp14;
      const auto tmp118 = tmp117 / tmp15;
      const auto tmp119 = std::pow( tmp118, tmp31 );
      const auto tmp120 = tmp35 * tmp119;
      const auto tmp121 = tmp120 * tmp29;
      const auto tmp122 = tmp114[ 1 ] / tmp15;
      const auto tmp123 = std::pow( tmp122, 2 );
      const auto tmp124 = std::pow( tmp118, tmp21 );
      const auto tmp125 = -1 * tmp124;
      const auto tmp126 = 1.0 + tmp125;
      const auto tmp127 = tmp126 * tmp123;
      const auto tmp128 = tmp127 / 0.0009;
      const auto tmp129 = tmp47 * tmp128;
      const auto tmp130 = tmp129 * tmp121;
      const auto tmp131 = tmp64 * tmp130;
      const auto tmp132 = tmp131 / tmp65;
      const auto tmp133 = tmp40[ 1 ] * tmp132;
      const auto tmp134 = tmp43 * tmp128;
      const auto tmp135 = tmp134 * tmp121;
      const auto tmp136 = tmp64 * tmp135;
      const auto tmp137 = tmp136 / tmp65;
      const auto tmp138 = tmp40[ 0 ] * tmp137;
      const auto tmp139 = tmp138 + tmp133;
      const auto tmp140 = -1 * tmp139;
      const auto tmp141 = tmp73 * tmp140;
      const auto tmp142 = tmp40[ 0 ] * tmp132;
      const auto tmp143 = tmp40[ 1 ] * tmp137;
      const auto tmp144 = tmp143 + tmp142;
      const auto tmp145 = -1 * tmp144;
      const auto tmp146 = tmp73 * tmp145;
      const auto tmp147 = -1 * tmp108;
      const auto tmp148 = tmp73 * tmp147;
      const auto tmp149 = -1 * tmp140;
      const auto tmp150 = tmp73 * tmp149;
      const auto tmp151 = -1 * tmp145;
      const auto tmp152 = tmp73 * tmp151;
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

} // namespace Integrands_53f62a551fba616bdfa07ec7b5bdc440v1_3

PYBIND11_MODULE( integrands_53f62a551fba616bdfa07ec7b5bdc440v1_3_220094643888a7089a01c0741f938529, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > > GridPart;
  typedef Integrands_53f62a551fba616bdfa07ec7b5bdc440v1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_53f62a551fba616bdfa07ec7b5bdc440v1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_53f62a551fba616bdfa07ec7b5bdc440v1_3_220094643888a7089a01c0741f938529.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffsolution, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffsolution_old, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffiterate ) { return new Integrands( coeffsolution, coeffsolution_old, coeffiterate ); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def_property( "timeStep", [] ( Integrands &self ) -> Integrands::ContimeStep { return self.contimeStep(); }, [] ( Integrands &self, const Integrands::ContimeStep &v ) { self.contimeStep() = v; } );
      cls.def_property( "penalty", [] ( Integrands &self ) -> Integrands::Conpenalty { return self.conpenalty(); }, [] ( Integrands &self, const Integrands::Conpenalty &v ) { self.conpenalty() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_53f62a551fba616bdfa07ec7b5bdc440v1_3
