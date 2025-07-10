#ifndef GuardIntegrands_9e780fcc181c319a61b8b9e91d9c556ev1_3
#define GuardIntegrands_9e780fcc181c319a61b8b9e91d9c556ev1_3
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

namespace Integrands_9e780fcc181c319a61b8b9e91d9c556ev1_3
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffiterate_comp_1, class Coeffiterate >
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
    typedef double Conpenalty;
    typedef std::tuple< std::shared_ptr< Conpenalty > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffiterate_comp_1, Coeffiterate > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffiterate_comp_1>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffiterate>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffiterate_comp_1 &coeffiterate_comp_1, const Coeffiterate &coeffiterate, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_{{ std::tuple< Dune::Fem::ConstLocalFunction< Coeffiterate_comp_1 >, Dune::Fem::ConstLocalFunction< Coeffiterate > >( Dune::Fem::ConstLocalFunction< Coeffiterate_comp_1 >( coeffiterate_comp_1 ), Dune::Fem::ConstLocalFunction< Coeffiterate >( coeffiterate ) ), std::tuple< Dune::Fem::ConstLocalFunction< Coeffiterate_comp_1 >, Dune::Fem::ConstLocalFunction< Coeffiterate > >( Dune::Fem::ConstLocalFunction< Coeffiterate_comp_1 >( coeffiterate_comp_1 ), Dune::Fem::ConstLocalFunction< Coeffiterate >( coeffiterate ) ) }}
    {
      std::get< 0 >( constants_ ) = std::make_shared< Conpenalty >( (Conpenalty(0)) );
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
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
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
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
        const Geometry& outsideGeom = entity_[ static_cast< std::size_t >( Side::out ) ].geometry();
        cellVolume_[ static_cast< std::size_t >( Side::out ) ] = outsideGeom.volume();
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
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
      const auto tmp6 = -0.7000000000000002 * tmp5;
      const auto tmp7 = 2.7 + tmp6;
      const auto tmp8 = 3.0 * tmp7;
      const auto tmp9 = 2.0 + tmp8;
      const auto tmp10 = tmp9 / tmp7;
      const auto tmp11 = -0.01999999999999999 * tmp5;
      const auto tmp12 = 0.12 + tmp11;
      const auto tmp13 = -1 * tmp12;
      const auto tmp14 = 1.0 + tmp13;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp15 = evaluateCoefficient< 1 >( x );
      const auto tmp16 = -1 * tmp15[ 1 ];
      const auto tmp17 = 1 + tmp16;
      const auto tmp18 = tmp17 + tmp13;
      const auto tmp19 = tmp18 / tmp14;
      const auto tmp20 = std::pow( tmp19, tmp10 );
      const auto tmp21 = tmp20 / 0.001;
      const auto tmp22 = 1000.0 * tmp21;
      const auto tmp23 = tmp15[ 1 ] / tmp14;
      const auto tmp24 = std::pow( tmp23, 2 );
      const auto tmp25 = 2.0 + tmp7;
      const auto tmp26 = tmp25 / tmp7;
      const auto tmp27 = std::pow( tmp19, tmp26 );
      const auto tmp28 = -1 * tmp27;
      const auto tmp29 = 1.0 + tmp28;
      const auto tmp30 = tmp29 * tmp24;
      const auto tmp31 = tmp30 / 0.0009;
      const auto tmp32 = 1460.0 * tmp31;
      const auto tmp33 = tmp32 + tmp22;
      const auto tmp34 = -9.81 * tmp33;
      const auto tmp35 = 5e-11 * tmp5;
      const auto tmp36 = -5e-11 + tmp35;
      const auto tmp37 = -1 * tmp36;
      const auto tmp38 = tmp37 * tmp34;
      const auto tmp39 = -1.0 / tmp14;
      const auto tmp40 = -1.0 / tmp7;
      const auto tmp41 = -1.0 + tmp40;
      const auto tmp42 = std::pow( tmp19, tmp41 );
      const auto tmp43 = 4245.0 * tmp5;
      const auto tmp44 = 755.0 + tmp43;
      const auto tmp45 = tmp44 * tmp40;
      const auto tmp46 = tmp45 * tmp42;
      const auto tmp47 = tmp46 * tmp39;
      const auto tmp48 = tmp47 * tmp31;
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp49 = jacobianCoefficient< 0 >( x );
      const auto tmp50 = (tmp49[ 0 ])[ 1 ] * tmp48;
      const auto tmp51 = tmp31 + tmp21;
      const auto tmp52 = std::get< 1 >( u );
      const auto tmp53 = (tmp52[ 0 ])[ 1 ] * tmp51;
      const auto tmp54 = tmp53 + tmp50;
      const auto tmp55 = tmp54 * tmp36;
      const auto tmp56 = -9.993360000000001e-11 * tmp5;
      const auto tmp57 = 1e-10 + tmp56;
      const auto tmp58 = (tmp49[ 0 ])[ 0 ] * tmp48;
      const auto tmp59 = (tmp52[ 0 ])[ 0 ] * tmp51;
      const auto tmp60 = tmp59 + tmp58;
      const auto tmp61 = tmp60 * tmp57;
      const auto tmp62 = tmp61 + tmp55;
      const auto tmp63 = tmp62 + tmp38;
      const auto tmp64 = -1 * tmp57;
      const auto tmp65 = tmp64 * tmp34;
      const auto tmp66 = tmp60 * tmp36;
      const auto tmp67 = tmp54 * tmp57;
      const auto tmp68 = tmp67 + tmp66;
      const auto tmp69 = tmp68 + tmp65;
      return RangeValueType{ { 0 }, { { tmp63, tmp69 } } };
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
      const auto tmp6 = -0.7000000000000002 * tmp5;
      const auto tmp7 = 2.7 + tmp6;
      const auto tmp8 = 3.0 * tmp7;
      const auto tmp9 = 2.0 + tmp8;
      const auto tmp10 = tmp9 / tmp7;
      const auto tmp11 = -0.01999999999999999 * tmp5;
      const auto tmp12 = 0.12 + tmp11;
      const auto tmp13 = -1 * tmp12;
      const auto tmp14 = 1.0 + tmp13;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp15 = evaluateCoefficient< 1 >( x );
      const auto tmp16 = -1 * tmp15[ 1 ];
      const auto tmp17 = 1 + tmp16;
      const auto tmp18 = tmp17 + tmp13;
      const auto tmp19 = tmp18 / tmp14;
      const auto tmp20 = std::pow( tmp19, tmp10 );
      const auto tmp21 = tmp20 / 0.001;
      const auto tmp22 = tmp15[ 1 ] / tmp14;
      const auto tmp23 = std::pow( tmp22, 2 );
      const auto tmp24 = 2.0 + tmp7;
      const auto tmp25 = tmp24 / tmp7;
      const auto tmp26 = std::pow( tmp19, tmp25 );
      const auto tmp27 = -1 * tmp26;
      const auto tmp28 = 1.0 + tmp27;
      const auto tmp29 = tmp28 * tmp23;
      const auto tmp30 = tmp29 / 0.0009;
      const auto tmp31 = tmp30 + tmp21;
      const auto tmp32 = -9.993360000000001e-11 * tmp5;
      const auto tmp33 = 1e-10 + tmp32;
      const auto tmp34 = tmp33 * tmp31;
      const auto tmp35 = 5e-11 * tmp5;
      const auto tmp36 = -5e-11 + tmp35;
      const auto tmp37 = tmp36 * tmp31;
      return [ tmp34, tmp37 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { tmp34 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp37 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], tmp37 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp34 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
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
      const auto tmp10 = -0.7000000000000002 * tmp9;
      const auto tmp11 = 2.7 + tmp10;
      const auto tmp12 = 3.0 * tmp11;
      const auto tmp13 = 2.0 + tmp12;
      const auto tmp14 = tmp13 / tmp11;
      const auto tmp15 = -0.01999999999999999 * tmp9;
      const auto tmp16 = 0.12 + tmp15;
      const auto tmp17 = -1 * tmp16;
      const auto tmp18 = 1.0 + tmp17;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp19 = evaluateCoefficient< 1 >( x );
      const auto tmp20 = -1 * tmp19[ 1 ];
      const auto tmp21 = 1 + tmp20;
      const auto tmp22 = tmp21 + tmp17;
      const auto tmp23 = tmp22 / tmp18;
      const auto tmp24 = std::pow( tmp23, tmp14 );
      const auto tmp25 = tmp24 / 0.001;
      const auto tmp26 = 1000.0 * tmp25;
      const auto tmp27 = tmp19[ 1 ] / tmp18;
      const auto tmp28 = std::pow( tmp27, 2 );
      const auto tmp29 = 2.0 + tmp11;
      const auto tmp30 = tmp29 / tmp11;
      const auto tmp31 = std::pow( tmp23, tmp30 );
      const auto tmp32 = -1 * tmp31;
      const auto tmp33 = 1.0 + tmp32;
      const auto tmp34 = tmp33 * tmp28;
      const auto tmp35 = tmp34 / 0.0009;
      const auto tmp36 = 1460.0 * tmp35;
      const auto tmp37 = tmp36 + tmp26;
      const auto tmp38 = -9.81 * tmp37;
      const auto tmp39 = 5e-11 * tmp9;
      const auto tmp40 = -5e-11 + tmp39;
      const auto tmp41 = -1 * tmp40;
      const auto tmp42 = tmp41 * tmp38;
      const auto tmp43 = -1.0 / tmp18;
      const auto tmp44 = -1.0 / tmp11;
      const auto tmp45 = -1.0 + tmp44;
      const auto tmp46 = std::pow( tmp23, tmp45 );
      const auto tmp47 = 4245.0 * tmp9;
      const auto tmp48 = 755.0 + tmp47;
      const auto tmp49 = tmp48 * tmp44;
      const auto tmp50 = tmp49 * tmp46;
      const auto tmp51 = tmp50 * tmp43;
      const auto tmp52 = tmp51 * tmp35;
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp53 = jacobianCoefficient< 0 >( x );
      const auto tmp54 = (tmp53[ 0 ])[ 1 ] * tmp52;
      const auto tmp55 = tmp35 + tmp25;
      const auto tmp56 = std::get< 1 >( u );
      const auto tmp57 = (tmp56[ 0 ])[ 1 ] * tmp55;
      const auto tmp58 = tmp57 + tmp54;
      const auto tmp59 = tmp58 * tmp40;
      const auto tmp60 = -9.993360000000001e-11 * tmp9;
      const auto tmp61 = 1e-10 + tmp60;
      const auto tmp62 = (tmp53[ 0 ])[ 0 ] * tmp52;
      const auto tmp63 = (tmp56[ 0 ])[ 0 ] * tmp55;
      const auto tmp64 = tmp63 + tmp62;
      const auto tmp65 = tmp64 * tmp61;
      const auto tmp66 = tmp65 + tmp59;
      const auto tmp67 = tmp66 + tmp42;
      GlobalCoordinateType tmp68 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp69 = tmp68[ 0 ] * tmp67;
      const auto tmp70 = -1 * tmp61;
      const auto tmp71 = tmp70 * tmp38;
      const auto tmp72 = tmp64 * tmp40;
      const auto tmp73 = tmp58 * tmp61;
      const auto tmp74 = tmp73 + tmp72;
      const auto tmp75 = tmp74 + tmp71;
      const auto tmp76 = tmp68[ 1 ] * tmp75;
      const auto tmp77 = tmp76 + tmp69;
      const auto tmp78 = tmp77 * tmp4;
      const auto tmp79 = -1 * tmp78;
      const auto tmp80 = -0.65 + tmp0[ 1 ];
      const auto tmp81 = std::abs( tmp80 );
      const auto tmp82 = (tmp6 < 0.06 ? 1.0 : 0.0) * (tmp81 < 1e-08 ? 1.0 : 0.0);
      const auto tmp83 = -5.136999999999999e-05 * tmp82;
      auto tmp84 = facetArea_;
      auto tmp85 = cellVolume();
      const auto tmp86 = tmp85 / tmp84;
      const auto tmp87 = tmp68[ 1 ] * tmp40;
      const auto tmp88 = tmp68[ 0 ] * tmp61;
      const auto tmp89 = tmp88 + tmp87;
      const auto tmp90 = tmp68[ 0 ] * tmp89;
      const auto tmp91 = tmp68[ 0 ] * tmp40;
      const auto tmp92 = tmp68[ 1 ] * tmp61;
      const auto tmp93 = tmp92 + tmp91;
      const auto tmp94 = tmp68[ 1 ] * tmp93;
      const auto tmp95 = tmp94 + tmp90;
      double tmp96 = constant< 0 >();
      const auto tmp97 = tmp96 * tmp95;
      const auto tmp98 = 0.5 + tmp17;
      const auto tmp99 = tmp98 / tmp18;
      const auto tmp100 = std::pow( tmp99, tmp14 );
      const auto tmp101 = tmp100 / 0.001;
      const auto tmp102 = 0.5 / tmp18;
      const auto tmp103 = std::pow( tmp102, 2 );
      const auto tmp104 = std::pow( tmp99, tmp30 );
      const auto tmp105 = -1 * tmp104;
      const auto tmp106 = 1.0 + tmp105;
      const auto tmp107 = tmp106 * tmp103;
      const auto tmp108 = tmp107 / 0.0009;
      const auto tmp109 = tmp108 + tmp101;
      const auto tmp110 = tmp109 * tmp97;
      const auto tmp111 = tmp110 / tmp86;
      const auto tmp112 = -1 * tmp0[ 1 ];
      const auto tmp113 = 0.65 + tmp112;
      const auto tmp114 = 9810.0 * tmp113;
      const auto tmp115 = -1 * tmp114;
      const auto tmp116 = std::get< 0 >( u );
      const auto tmp117 = tmp116[ 0 ] + tmp115;
      const auto tmp118 = tmp117 * tmp111;
      const auto tmp119 = tmp4 * tmp118;
      const auto tmp120 = tmp119 + tmp83;
      const auto tmp121 = tmp120 + tmp79;
      return RangeValueType{ { tmp121 }, { { 0, 0 } } };
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
      const auto tmp9 = 5e-11 * tmp8;
      const auto tmp10 = -5e-11 + tmp9;
      GlobalCoordinateType tmp11 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp12 = tmp11[ 1 ] * tmp10;
      const auto tmp13 = -9.993360000000001e-11 * tmp8;
      const auto tmp14 = 1e-10 + tmp13;
      const auto tmp15 = tmp11[ 0 ] * tmp14;
      const auto tmp16 = tmp15 + tmp12;
      const auto tmp17 = tmp11[ 0 ] * tmp16;
      const auto tmp18 = tmp11[ 0 ] * tmp10;
      const auto tmp19 = tmp11[ 1 ] * tmp14;
      const auto tmp20 = tmp19 + tmp18;
      const auto tmp21 = tmp11[ 1 ] * tmp20;
      const auto tmp22 = tmp21 + tmp17;
      double tmp23 = constant< 0 >();
      const auto tmp24 = tmp23 * tmp22;
      const auto tmp25 = -0.7000000000000002 * tmp8;
      const auto tmp26 = 2.7 + tmp25;
      const auto tmp27 = 3.0 * tmp26;
      const auto tmp28 = 2.0 + tmp27;
      const auto tmp29 = tmp28 / tmp26;
      const auto tmp30 = -0.01999999999999999 * tmp8;
      const auto tmp31 = 0.12 + tmp30;
      const auto tmp32 = -1 * tmp31;
      const auto tmp33 = 1.0 + tmp32;
      const auto tmp34 = 0.5 + tmp32;
      const auto tmp35 = tmp34 / tmp33;
      const auto tmp36 = std::pow( tmp35, tmp29 );
      const auto tmp37 = tmp36 / 0.001;
      const auto tmp38 = 0.5 / tmp33;
      const auto tmp39 = std::pow( tmp38, 2 );
      const auto tmp40 = 2.0 + tmp26;
      const auto tmp41 = tmp40 / tmp26;
      const auto tmp42 = std::pow( tmp35, tmp41 );
      const auto tmp43 = -1 * tmp42;
      const auto tmp44 = 1.0 + tmp43;
      const auto tmp45 = tmp44 * tmp39;
      const auto tmp46 = tmp45 / 0.0009;
      const auto tmp47 = tmp46 + tmp37;
      const auto tmp48 = tmp47 * tmp24;
      const auto tmp49 = tmp48 / tmp2;
      const auto tmp50 = -0.9 + tmp3[ 0 ];
      const auto tmp51 = std::abs( tmp50 );
      const auto tmp52 = std::abs( tmp3[ 0 ] );
      const auto tmp53 = (tmp52 < 1e-08 ? 1.0 : 0.0) + (tmp51 < 1e-08 ? 1.0 : 0.0);
      const auto tmp54 = tmp53 * tmp49;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp55 = evaluateCoefficient< 1 >( x );
      const auto tmp56 = -1 * tmp55[ 1 ];
      const auto tmp57 = 1 + tmp56;
      const auto tmp58 = tmp57 + tmp32;
      const auto tmp59 = tmp58 / tmp33;
      const auto tmp60 = std::pow( tmp59, tmp29 );
      const auto tmp61 = tmp60 / 0.001;
      const auto tmp62 = tmp55[ 1 ] / tmp33;
      const auto tmp63 = std::pow( tmp62, 2 );
      const auto tmp64 = std::pow( tmp59, tmp41 );
      const auto tmp65 = -1 * tmp64;
      const auto tmp66 = 1.0 + tmp65;
      const auto tmp67 = tmp66 * tmp63;
      const auto tmp68 = tmp67 / 0.0009;
      const auto tmp69 = tmp68 + tmp61;
      const auto tmp70 = tmp10 * tmp69;
      const auto tmp71 = tmp11[ 1 ] * tmp70;
      const auto tmp72 = tmp14 * tmp69;
      const auto tmp73 = tmp11[ 0 ] * tmp72;
      const auto tmp74 = tmp73 + tmp71;
      const auto tmp75 = tmp74 * tmp53;
      const auto tmp76 = -1 * tmp75;
      const auto tmp77 = tmp11[ 0 ] * tmp70;
      const auto tmp78 = tmp11[ 1 ] * tmp72;
      const auto tmp79 = tmp78 + tmp77;
      const auto tmp80 = tmp79 * tmp53;
      const auto tmp81 = -1 * tmp80;
      return [ tmp54, tmp76, tmp81 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (tmp54 * (std::get< 0 >( phi ))[ 0 ] + tmp76 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp81 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { { 0, 0 } } };
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
      const auto tmp12 = -0.7000000000000002 * tmp11;
      const auto tmp13 = 2.7 + tmp12;
      const auto tmp14 = 3.0 * tmp13;
      const auto tmp15 = 2.0 + tmp14;
      const auto tmp16 = tmp15 / tmp13;
      const auto tmp17 = -0.01999999999999999 * tmp11;
      const auto tmp18 = 0.12 + tmp17;
      const auto tmp19 = -1 * tmp18;
      const auto tmp20 = 1.0 + tmp19;
      const auto tmp21 = 0.5 + tmp19;
      const auto tmp22 = tmp21 / tmp20;
      const auto tmp23 = std::pow( tmp22, tmp16 );
      const auto tmp24 = tmp23 / 0.001;
      const auto tmp25 = 0.5 / tmp20;
      const auto tmp26 = std::pow( tmp25, 2 );
      const auto tmp27 = 2.0 + tmp13;
      const auto tmp28 = tmp27 / tmp13;
      const auto tmp29 = std::pow( tmp22, tmp28 );
      const auto tmp30 = -1 * tmp29;
      const auto tmp31 = 1.0 + tmp30;
      const auto tmp32 = tmp31 * tmp26;
      const auto tmp33 = tmp32 / 0.0009;
      const auto tmp34 = tmp33 + tmp24;
      const auto tmp35 = std::max( tmp34, tmp34 );
      GlobalCoordinateType tmp36 = intersection_.unitOuterNormal( xIn.localPosition() );
      const auto tmp37 = -1 * tmp36[ 1 ];
      const auto tmp38 = -9.993360000000001e-11 * tmp11;
      const auto tmp39 = 1e-10 + tmp38;
      const auto tmp40 = tmp39 * tmp37;
      const auto tmp41 = -1 * tmp36[ 0 ];
      const auto tmp42 = 5e-11 * tmp11;
      const auto tmp43 = -5e-11 + tmp42;
      const auto tmp44 = tmp43 * tmp41;
      const auto tmp45 = tmp44 + tmp40;
      const auto tmp46 = tmp45 * tmp37;
      const auto tmp47 = tmp43 * tmp37;
      const auto tmp48 = tmp39 * tmp41;
      const auto tmp49 = tmp48 + tmp47;
      const auto tmp50 = tmp49 * tmp41;
      const auto tmp51 = tmp50 + tmp46;
      const auto tmp52 = tmp36[ 1 ] * tmp43;
      const auto tmp53 = tmp36[ 0 ] * tmp39;
      const auto tmp54 = tmp53 + tmp52;
      const auto tmp55 = tmp36[ 0 ] * tmp54;
      const auto tmp56 = tmp36[ 0 ] * tmp43;
      const auto tmp57 = tmp36[ 1 ] * tmp39;
      const auto tmp58 = tmp57 + tmp56;
      const auto tmp59 = tmp36[ 1 ] * tmp58;
      const auto tmp60 = tmp59 + tmp55;
      const auto tmp61 = tmp60 + tmp51;
      const auto tmp62 = 0.5 * tmp61;
      const auto tmp63 = tmp60 * tmp51;
      const auto tmp64 = tmp63 / tmp62;
      double tmp65 = constant< 0 >();
      const auto tmp66 = tmp65 * tmp64;
      const auto tmp67 = tmp66 * tmp35;
      const auto tmp68 = tmp67 / tmp5;
      const auto tmp69 = std::get< 0 >( uOut );
      const auto tmp70 = -1 * tmp69[ 0 ];
      const auto tmp71 = std::get< 0 >( uIn );
      const auto tmp72 = tmp71[ 0 ] + tmp70;
      const auto tmp73 = tmp72 * tmp68;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp74 = evaluateCoefficient< 1, Side::out >( xOut );
      const auto tmp75 = -1 * tmp74[ 1 ];
      const auto tmp76 = 1 + tmp75;
      const auto tmp77 = tmp76 + tmp19;
      const auto tmp78 = tmp77 / tmp20;
      const auto tmp79 = std::pow( tmp78, tmp16 );
      const auto tmp80 = tmp79 / 0.001;
      const auto tmp81 = 1000.0 * tmp80;
      const auto tmp82 = tmp74[ 1 ] / tmp20;
      const auto tmp83 = std::pow( tmp82, 2 );
      const auto tmp84 = std::pow( tmp78, tmp28 );
      const auto tmp85 = -1 * tmp84;
      const auto tmp86 = 1.0 + tmp85;
      const auto tmp87 = tmp86 * tmp83;
      const auto tmp88 = tmp87 / 0.0009;
      const auto tmp89 = 1460.0 * tmp88;
      const auto tmp90 = tmp89 + tmp81;
      const auto tmp91 = -9.81 * tmp90;
      const auto tmp92 = -1 * tmp43;
      const auto tmp93 = tmp92 * tmp91;
      const auto tmp94 = -1.0 / tmp20;
      const auto tmp95 = -1.0 / tmp13;
      const auto tmp96 = -1.0 + tmp95;
      const auto tmp97 = std::pow( tmp78, tmp96 );
      const auto tmp98 = 4245.0 * tmp11;
      const auto tmp99 = 755.0 + tmp98;
      const auto tmp100 = tmp99 * tmp95;
      const auto tmp101 = tmp100 * tmp97;
      const auto tmp102 = tmp101 * tmp94;
      const auto tmp103 = tmp102 * tmp88;
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp104 = jacobianCoefficient< 0, Side::out >( xOut );
      const auto tmp105 = (tmp104[ 0 ])[ 1 ] * tmp103;
      const auto tmp106 = tmp88 + tmp80;
      const auto tmp107 = std::get< 1 >( uOut );
      const auto tmp108 = (tmp107[ 0 ])[ 1 ] * tmp106;
      const auto tmp109 = tmp108 + tmp105;
      const auto tmp110 = tmp109 * tmp43;
      const auto tmp111 = (tmp104[ 0 ])[ 0 ] * tmp103;
      const auto tmp112 = (tmp107[ 0 ])[ 0 ] * tmp106;
      const auto tmp113 = tmp112 + tmp111;
      const auto tmp114 = tmp113 * tmp39;
      const auto tmp115 = tmp114 + tmp110;
      const auto tmp116 = tmp115 + tmp93;
      const auto tmp117 = tmp60 * tmp116;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp118 = evaluateCoefficient< 1, Side::in >( xIn );
      const auto tmp119 = -1 * tmp118[ 1 ];
      const auto tmp120 = 1 + tmp119;
      const auto tmp121 = tmp120 + tmp19;
      const auto tmp122 = tmp121 / tmp20;
      const auto tmp123 = std::pow( tmp122, tmp16 );
      const auto tmp124 = tmp123 / 0.001;
      const auto tmp125 = 1000.0 * tmp124;
      const auto tmp126 = tmp118[ 1 ] / tmp20;
      const auto tmp127 = std::pow( tmp126, 2 );
      const auto tmp128 = std::pow( tmp122, tmp28 );
      const auto tmp129 = -1 * tmp128;
      const auto tmp130 = 1.0 + tmp129;
      const auto tmp131 = tmp130 * tmp127;
      const auto tmp132 = tmp131 / 0.0009;
      const auto tmp133 = 1460.0 * tmp132;
      const auto tmp134 = tmp133 + tmp125;
      const auto tmp135 = -9.81 * tmp134;
      const auto tmp136 = tmp92 * tmp135;
      const auto tmp137 = std::pow( tmp122, tmp96 );
      const auto tmp138 = tmp100 * tmp137;
      const auto tmp139 = tmp138 * tmp94;
      const auto tmp140 = tmp139 * tmp132;
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp141 = jacobianCoefficient< 0, Side::in >( xIn );
      const auto tmp142 = (tmp141[ 0 ])[ 1 ] * tmp140;
      const auto tmp143 = tmp132 + tmp124;
      const auto tmp144 = std::get< 1 >( uIn );
      const auto tmp145 = (tmp144[ 0 ])[ 1 ] * tmp143;
      const auto tmp146 = tmp145 + tmp142;
      const auto tmp147 = tmp146 * tmp43;
      const auto tmp148 = (tmp141[ 0 ])[ 0 ] * tmp140;
      const auto tmp149 = (tmp144[ 0 ])[ 0 ] * tmp143;
      const auto tmp150 = tmp149 + tmp148;
      const auto tmp151 = tmp150 * tmp39;
      const auto tmp152 = tmp151 + tmp147;
      const auto tmp153 = tmp152 + tmp136;
      const auto tmp154 = tmp51 * tmp153;
      const auto tmp155 = tmp154 + tmp117;
      const auto tmp156 = tmp155 / tmp61;
      const auto tmp157 = tmp36[ 0 ] * tmp156;
      const auto tmp158 = -1 * tmp39;
      const auto tmp159 = tmp158 * tmp91;
      const auto tmp160 = tmp113 * tmp43;
      const auto tmp161 = tmp109 * tmp39;
      const auto tmp162 = tmp161 + tmp160;
      const auto tmp163 = tmp162 + tmp159;
      const auto tmp164 = tmp60 * tmp163;
      const auto tmp165 = tmp158 * tmp135;
      const auto tmp166 = tmp150 * tmp43;
      const auto tmp167 = tmp146 * tmp39;
      const auto tmp168 = tmp167 + tmp166;
      const auto tmp169 = tmp168 + tmp165;
      const auto tmp170 = tmp51 * tmp169;
      const auto tmp171 = tmp170 + tmp164;
      const auto tmp172 = tmp171 / tmp61;
      const auto tmp173 = tmp36[ 1 ] * tmp172;
      const auto tmp174 = tmp173 + tmp157;
      const auto tmp175 = -1 * tmp174;
      const auto tmp176 = tmp175 + tmp73;
      const auto tmp177 = -1 * tmp73;
      const auto tmp178 = -1 * tmp175;
      const auto tmp179 = tmp178 + tmp177;
      return std::make_pair( RangeValueType{ { tmp176 }, { { 0, 0 } } }, RangeValueType{ { tmp179 }, { { 0, 0 } } } );
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
      const auto tmp12 = -0.7000000000000002 * tmp11;
      const auto tmp13 = 2.7 + tmp12;
      const auto tmp14 = 3.0 * tmp13;
      const auto tmp15 = 2.0 + tmp14;
      const auto tmp16 = tmp15 / tmp13;
      const auto tmp17 = -0.01999999999999999 * tmp11;
      const auto tmp18 = 0.12 + tmp17;
      const auto tmp19 = -1 * tmp18;
      const auto tmp20 = 1.0 + tmp19;
      const auto tmp21 = 0.5 + tmp19;
      const auto tmp22 = tmp21 / tmp20;
      const auto tmp23 = std::pow( tmp22, tmp16 );
      const auto tmp24 = tmp23 / 0.001;
      const auto tmp25 = 0.5 / tmp20;
      const auto tmp26 = std::pow( tmp25, 2 );
      const auto tmp27 = 2.0 + tmp13;
      const auto tmp28 = tmp27 / tmp13;
      const auto tmp29 = std::pow( tmp22, tmp28 );
      const auto tmp30 = -1 * tmp29;
      const auto tmp31 = 1.0 + tmp30;
      const auto tmp32 = tmp31 * tmp26;
      const auto tmp33 = tmp32 / 0.0009;
      const auto tmp34 = tmp33 + tmp24;
      const auto tmp35 = std::max( tmp34, tmp34 );
      GlobalCoordinateType tmp36 = intersection_.unitOuterNormal( xIn.localPosition() );
      const auto tmp37 = -1 * tmp36[ 1 ];
      const auto tmp38 = -9.993360000000001e-11 * tmp11;
      const auto tmp39 = 1e-10 + tmp38;
      const auto tmp40 = tmp39 * tmp37;
      const auto tmp41 = -1 * tmp36[ 0 ];
      const auto tmp42 = 5e-11 * tmp11;
      const auto tmp43 = -5e-11 + tmp42;
      const auto tmp44 = tmp43 * tmp41;
      const auto tmp45 = tmp44 + tmp40;
      const auto tmp46 = tmp45 * tmp37;
      const auto tmp47 = tmp43 * tmp37;
      const auto tmp48 = tmp39 * tmp41;
      const auto tmp49 = tmp48 + tmp47;
      const auto tmp50 = tmp49 * tmp41;
      const auto tmp51 = tmp50 + tmp46;
      const auto tmp52 = tmp36[ 1 ] * tmp43;
      const auto tmp53 = tmp36[ 0 ] * tmp39;
      const auto tmp54 = tmp53 + tmp52;
      const auto tmp55 = tmp36[ 0 ] * tmp54;
      const auto tmp56 = tmp36[ 0 ] * tmp43;
      const auto tmp57 = tmp36[ 1 ] * tmp39;
      const auto tmp58 = tmp57 + tmp56;
      const auto tmp59 = tmp36[ 1 ] * tmp58;
      const auto tmp60 = tmp59 + tmp55;
      const auto tmp61 = tmp60 + tmp51;
      const auto tmp62 = 0.5 * tmp61;
      const auto tmp63 = tmp60 * tmp51;
      const auto tmp64 = tmp63 / tmp62;
      double tmp65 = constant< 0 >();
      const auto tmp66 = tmp65 * tmp64;
      const auto tmp67 = tmp66 * tmp35;
      const auto tmp68 = tmp67 / tmp5;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp69 = evaluateCoefficient< 1, Side::in >( xIn );
      const auto tmp70 = -1 * tmp69[ 1 ];
      const auto tmp71 = 1 + tmp70;
      const auto tmp72 = tmp71 + tmp19;
      const auto tmp73 = tmp72 / tmp20;
      const auto tmp74 = std::pow( tmp73, tmp16 );
      const auto tmp75 = tmp74 / 0.001;
      const auto tmp76 = tmp69[ 1 ] / tmp20;
      const auto tmp77 = std::pow( tmp76, 2 );
      const auto tmp78 = std::pow( tmp73, tmp28 );
      const auto tmp79 = -1 * tmp78;
      const auto tmp80 = 1.0 + tmp79;
      const auto tmp81 = tmp80 * tmp77;
      const auto tmp82 = tmp81 / 0.0009;
      const auto tmp83 = tmp82 + tmp75;
      const auto tmp84 = tmp43 * tmp83;
      const auto tmp85 = tmp51 * tmp84;
      const auto tmp86 = tmp85 / tmp61;
      const auto tmp87 = tmp36[ 1 ] * tmp86;
      const auto tmp88 = tmp39 * tmp83;
      const auto tmp89 = tmp51 * tmp88;
      const auto tmp90 = tmp89 / tmp61;
      const auto tmp91 = tmp36[ 0 ] * tmp90;
      const auto tmp92 = tmp91 + tmp87;
      const auto tmp93 = -1 * tmp92;
      const auto tmp94 = tmp36[ 0 ] * tmp86;
      const auto tmp95 = tmp36[ 1 ] * tmp90;
      const auto tmp96 = tmp95 + tmp94;
      const auto tmp97 = -1 * tmp96;
      const auto tmp98 = -1 * tmp68;
      const auto tmp99 = -1 * tmp93;
      const auto tmp100 = -1 * tmp97;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp101 = evaluateCoefficient< 1, Side::out >( xOut );
      const auto tmp102 = -1 * tmp101[ 1 ];
      const auto tmp103 = 1 + tmp102;
      const auto tmp104 = tmp103 + tmp19;
      const auto tmp105 = tmp104 / tmp20;
      const auto tmp106 = std::pow( tmp105, tmp16 );
      const auto tmp107 = tmp106 / 0.001;
      const auto tmp108 = tmp101[ 1 ] / tmp20;
      const auto tmp109 = std::pow( tmp108, 2 );
      const auto tmp110 = std::pow( tmp105, tmp28 );
      const auto tmp111 = -1 * tmp110;
      const auto tmp112 = 1.0 + tmp111;
      const auto tmp113 = tmp112 * tmp109;
      const auto tmp114 = tmp113 / 0.0009;
      const auto tmp115 = tmp114 + tmp107;
      const auto tmp116 = tmp43 * tmp115;
      const auto tmp117 = tmp60 * tmp116;
      const auto tmp118 = tmp117 / tmp61;
      const auto tmp119 = tmp36[ 1 ] * tmp118;
      const auto tmp120 = tmp39 * tmp115;
      const auto tmp121 = tmp60 * tmp120;
      const auto tmp122 = tmp121 / tmp61;
      const auto tmp123 = tmp36[ 0 ] * tmp122;
      const auto tmp124 = tmp123 + tmp119;
      const auto tmp125 = -1 * tmp124;
      const auto tmp126 = tmp36[ 0 ] * tmp118;
      const auto tmp127 = tmp36[ 1 ] * tmp122;
      const auto tmp128 = tmp127 + tmp126;
      const auto tmp129 = -1 * tmp128;
      const auto tmp130 = -1 * tmp98;
      const auto tmp131 = -1 * tmp125;
      const auto tmp132 = -1 * tmp129;
      return std::make_pair( [ tmp100, tmp68, tmp93, tmp97, tmp98, tmp99 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { (tmp68 * (std::get< 0 >( phiIn ))[ 0 ] + tmp93 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp97 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ] }, { { 0, 0 } } }, RangeValueType{ { (tmp98 * (std::get< 0 >( phiIn ))[ 0 ] + tmp99 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp100 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ] }, { { 0, 0 } } } );
        }, [ tmp125, tmp129, tmp130, tmp131, tmp132, tmp98 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { (tmp98 * (std::get< 0 >( phiOut ))[ 0 ] + tmp125 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp129 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ] }, { { 0, 0 } } }, RangeValueType{ { (tmp130 * (std::get< 0 >( phiOut ))[ 0 ] + tmp131 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp132 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ] }, { { 0, 0 } } } );
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

    const Conpenalty &conpenalty () const
    {
      return *std::get< 0 >( constants_ );
    }

    Conpenalty &conpenalty ()
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
    std::array< std::tuple< Dune::Fem::ConstLocalFunction< Coeffiterate_comp_1 >, Dune::Fem::ConstLocalFunction< Coeffiterate > >, 2 > coefficients_;
  };

} // namespace Integrands_9e780fcc181c319a61b8b9e91d9c556ev1_3

PYBIND11_MODULE( integrands_9e780fcc181c319a61b8b9e91d9c556ev1_3_220094643888a7089a01c0741f938529, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > > GridPart;
  typedef Integrands_9e780fcc181c319a61b8b9e91d9c556ev1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_9e780fcc181c319a61b8b9e91d9c556ev1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_9e780fcc181c319a61b8b9e91d9c556ev1_3_220094643888a7089a01c0741f938529.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 1 > > &coeffiterate_comp_1, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffiterate ) { return new Integrands( coeffiterate_comp_1, coeffiterate ); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property( "penalty", [] ( Integrands &self ) -> Integrands::Conpenalty { return self.conpenalty(); }, [] ( Integrands &self, const Integrands::Conpenalty &v ) { self.conpenalty() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_9e780fcc181c319a61b8b9e91d9c556ev1_3
