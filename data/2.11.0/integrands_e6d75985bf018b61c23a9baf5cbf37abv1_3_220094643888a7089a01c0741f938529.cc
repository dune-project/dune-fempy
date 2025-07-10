#ifndef GuardIntegrands_e6d75985bf018b61c23a9baf5cbf37abv1_3
#define GuardIntegrands_e6d75985bf018b61c23a9baf5cbf37abv1_3
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

namespace Integrands_e6d75985bf018b61c23a9baf5cbf37abv1_3
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffsolution_old >
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
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffsolution_old > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsolution_old>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffsolution_old &coeffsolution_old, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_{{ std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution_old > >( Dune::Fem::ConstLocalFunction< Coeffsolution_old >( coeffsolution_old ) ), std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution_old > >( Dune::Fem::ConstLocalFunction< Coeffsolution_old >( coeffsolution_old ) ) }}
    {
      std::get< 0 >( constants_ ) = std::make_shared< ContimeStep >( (ContimeStep(0)) );
      std::get< 1 >( constants_ ) = std::make_shared< Conpenalty >( (Conpenalty(0)) );
    }

    bool init ( const EntityType &entity )
    {
      entity_[ (static_cast< std::size_t >( Side::in )) ] = entity;
      geometry_.emplace( this->entity().geometry() );
      cellVolume_[ (static_cast< std::size_t >( Side::in )) ] = (geometry().volume());
      maxCellEdgeLength_[ (static_cast< std::size_t >( Side::in )) ] = Dune::FemPy::maxEdgeLength( (geometry()) );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
    }

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      const auto& intersectionGeom = intersection_.geometry();
      facetArea_ = (intersectionGeom.volume());
      entity_[ (static_cast< std::size_t >( Side::in )) ] = (intersection.inside());
      geometry_.emplace( this->entity().geometry() );
      cellVolume_[ (static_cast< std::size_t >( Side::in )) ] = (geometry().volume());
      maxCellEdgeLength_[ (static_cast< std::size_t >( Side::in )) ] = Dune::FemPy::maxEdgeLength( (geometry()) );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
        const Geometry& outsideGeom = entity_[ static_cast< std::size_t >( Side::out ) ].geometry();
        cellVolume_[ static_cast< std::size_t >( Side::out ) ] = outsideGeom.volume();
      maxCellEdgeLength_[ (static_cast< std::size_t >( Side::out )) ] = Dune::FemPy::maxEdgeLength( (outsideGeom) );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
      }
      return true;
    }
    typedef std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 > > RangeValueType;
    static constexpr bool _nonlinear = true;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      double tmp0 = constant< 0 >();
      GlobalCoordinateType tmp1 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp2 = -0.45 + tmp1[ 0 ];
      const auto tmp3 = std::abs( tmp2 );
      const auto tmp4 = -0.49 + tmp1[ 1 ];
      const auto tmp5 = std::abs( tmp4 );
      const auto tmp6 = (tmp5 < 0.03 ? 1.0 : 0.0) * (tmp3 < 0.11 ? 1.0 : 0.0);
      const auto tmp7 = -0.010000000000000009 * tmp6;
      const auto tmp8 = 0.4 + tmp7;
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp9 = evaluateCoefficient< 0 >( x );
      const auto tmp10 = -1 * tmp9[ 1 ];
      const auto tmp11 = std::get< 0 >( u );
      const auto tmp12 = tmp11[ 1 ] + tmp10;
      const auto tmp13 = tmp12 * tmp8;
      const auto tmp14 = tmp13 / tmp0;
      const auto tmp15 = -0.01999999999999999 * tmp6;
      const auto tmp16 = 0.12 + tmp15;
      const auto tmp17 = -1 * tmp16;
      const auto tmp18 = 1.0 + tmp17;
      const auto tmp19 = tmp11[ 1 ] / tmp18;
      const auto tmp20 = std::pow( tmp19, 2 );
      const auto tmp21 = -0.7000000000000002 * tmp6;
      const auto tmp22 = 2.7 + tmp21;
      const auto tmp23 = 2.0 + tmp22;
      const auto tmp24 = tmp23 / tmp22;
      const auto tmp25 = -1 * tmp11[ 1 ];
      const auto tmp26 = 1 + tmp25;
      const auto tmp27 = tmp26 + tmp17;
      const auto tmp28 = tmp27 / tmp18;
      const auto tmp29 = std::pow( tmp28, tmp24 );
      const auto tmp30 = -1 * tmp29;
      const auto tmp31 = 1.0 + tmp30;
      const auto tmp32 = tmp31 * tmp20;
      const auto tmp33 = tmp32 / 0.0009;
      const auto tmp34 = 5e-11 * tmp6;
      const auto tmp35 = -5e-11 + tmp34;
      const auto tmp36 = tmp35 * tmp33;
      const auto tmp37 = std::get< 2 >( u );
      const auto tmp38 = ((tmp37[ 0 ])[ 1 ])[ 0 ] * tmp36;
      const auto tmp39 = -1 + tmp24;
      const auto tmp40 = std::pow( tmp28, tmp39 );
      const auto tmp41 = std::get< 1 >( u );
      const auto tmp42 = -1 * (tmp41[ 1 ])[ 0 ];
      const auto tmp43 = tmp42 / tmp18;
      const auto tmp44 = tmp43 * tmp24;
      const auto tmp45 = tmp44 * tmp40;
      const auto tmp46 = -1 * tmp45;
      const auto tmp47 = tmp46 * tmp20;
      const auto tmp48 = (tmp41[ 1 ])[ 0 ] / tmp18;
      const auto tmp49 = 2 * tmp48;
      const auto tmp50 = tmp49 * tmp19;
      const auto tmp51 = tmp31 * tmp50;
      const auto tmp52 = tmp51 + tmp47;
      const auto tmp53 = tmp52 / 0.0009;
      const auto tmp54 = tmp35 * tmp53;
      const auto tmp55 = 14322.6 + (tmp41[ 0 ])[ 1 ];
      const auto tmp56 = tmp55 * tmp54;
      const auto tmp57 = tmp56 + tmp38;
      const auto tmp58 = -9.993360000000001e-11 * tmp6;
      const auto tmp59 = 1e-10 + tmp58;
      const auto tmp60 = tmp59 * tmp33;
      const auto tmp61 = ((tmp37[ 0 ])[ 0 ])[ 0 ] * tmp60;
      const auto tmp62 = tmp59 * tmp53;
      const auto tmp63 = (tmp41[ 0 ])[ 0 ] * tmp62;
      const auto tmp64 = tmp63 + tmp61;
      const auto tmp65 = tmp64 + tmp57;
      const auto tmp66 = -1.0 / tmp18;
      const auto tmp67 = -1.0 / tmp22;
      const auto tmp68 = -1.0 + tmp67;
      const auto tmp69 = std::pow( tmp28, tmp68 );
      const auto tmp70 = 4245.0 * tmp6;
      const auto tmp71 = 755.0 + tmp70;
      const auto tmp72 = tmp71 * tmp67;
      const auto tmp73 = tmp72 * tmp69;
      const auto tmp74 = tmp73 * tmp66;
      const auto tmp75 = tmp36 * tmp74;
      const auto tmp76 = ((tmp37[ 1 ])[ 1 ])[ 0 ] * tmp75;
      const auto tmp77 = tmp54 * tmp74;
      const auto tmp78 = -1 + tmp68;
      const auto tmp79 = std::pow( tmp28, tmp78 );
      const auto tmp80 = tmp68 * tmp43;
      const auto tmp81 = tmp80 * tmp79;
      const auto tmp82 = tmp72 * tmp81;
      const auto tmp83 = tmp82 * tmp66;
      const auto tmp84 = tmp36 * tmp83;
      const auto tmp85 = tmp84 + tmp77;
      const auto tmp86 = (tmp41[ 1 ])[ 1 ] * tmp85;
      const auto tmp87 = tmp86 + tmp76;
      const auto tmp88 = tmp60 * tmp74;
      const auto tmp89 = ((tmp37[ 1 ])[ 0 ])[ 0 ] * tmp88;
      const auto tmp90 = tmp62 * tmp74;
      const auto tmp91 = tmp60 * tmp83;
      const auto tmp92 = tmp91 + tmp90;
      const auto tmp93 = (tmp41[ 1 ])[ 0 ] * tmp92;
      const auto tmp94 = tmp93 + tmp89;
      const auto tmp95 = tmp94 + tmp87;
      const auto tmp96 = tmp95 + tmp65;
      const auto tmp97 = ((tmp37[ 0 ])[ 0 ])[ 1 ] * tmp36;
      const auto tmp98 = -1 * (tmp41[ 1 ])[ 1 ];
      const auto tmp99 = tmp98 / tmp18;
      const auto tmp100 = tmp99 * tmp24;
      const auto tmp101 = tmp100 * tmp40;
      const auto tmp102 = -1 * tmp101;
      const auto tmp103 = tmp102 * tmp20;
      const auto tmp104 = (tmp41[ 1 ])[ 1 ] / tmp18;
      const auto tmp105 = 2 * tmp104;
      const auto tmp106 = tmp105 * tmp19;
      const auto tmp107 = tmp31 * tmp106;
      const auto tmp108 = tmp107 + tmp103;
      const auto tmp109 = tmp108 / 0.0009;
      const auto tmp110 = tmp35 * tmp109;
      const auto tmp111 = (tmp41[ 0 ])[ 0 ] * tmp110;
      const auto tmp112 = tmp111 + tmp97;
      const auto tmp113 = ((tmp37[ 0 ])[ 1 ])[ 1 ] * tmp60;
      const auto tmp114 = tmp59 * tmp109;
      const auto tmp115 = tmp55 * tmp114;
      const auto tmp116 = tmp115 + tmp113;
      const auto tmp117 = tmp116 + tmp112;
      const auto tmp118 = ((tmp37[ 1 ])[ 0 ])[ 1 ] * tmp75;
      const auto tmp119 = tmp110 * tmp74;
      const auto tmp120 = tmp68 * tmp99;
      const auto tmp121 = tmp120 * tmp79;
      const auto tmp122 = tmp72 * tmp121;
      const auto tmp123 = tmp122 * tmp66;
      const auto tmp124 = tmp36 * tmp123;
      const auto tmp125 = tmp124 + tmp119;
      const auto tmp126 = (tmp41[ 1 ])[ 0 ] * tmp125;
      const auto tmp127 = tmp126 + tmp118;
      const auto tmp128 = ((tmp37[ 1 ])[ 1 ])[ 1 ] * tmp88;
      const auto tmp129 = tmp114 * tmp74;
      const auto tmp130 = tmp60 * tmp123;
      const auto tmp131 = tmp130 + tmp129;
      const auto tmp132 = (tmp41[ 1 ])[ 1 ] * tmp131;
      const auto tmp133 = tmp132 + tmp128;
      const auto tmp134 = tmp133 + tmp127;
      const auto tmp135 = tmp134 + tmp117;
      const auto tmp136 = tmp135 + tmp96;
      const auto tmp137 = -1 * tmp136;
      const auto tmp138 = tmp137 + tmp14;
      const auto tmp139 = std::pow( tmp138, 2 );
      auto tmp140 = maxCellEdgeLength();
      const auto tmp141 = std::pow( tmp140, 2 );
      const auto tmp142 = tmp141 * tmp139;
      return RangeValueType{ { tmp142 } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      auto tmp0 = maxCellEdgeLength();
      const auto tmp1 = std::pow( tmp0, 2 );
      double tmp2 = constant< 0 >();
      GlobalCoordinateType tmp3 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp4 = -0.45 + tmp3[ 0 ];
      const auto tmp5 = std::abs( tmp4 );
      const auto tmp6 = -0.49 + tmp3[ 1 ];
      const auto tmp7 = std::abs( tmp6 );
      const auto tmp8 = (tmp7 < 0.03 ? 1.0 : 0.0) * (tmp5 < 0.11 ? 1.0 : 0.0);
      const auto tmp9 = -0.010000000000000009 * tmp8;
      const auto tmp10 = 0.4 + tmp9;
      const auto tmp11 = tmp10 / tmp2;
      const auto tmp12 = -0.01999999999999999 * tmp8;
      const auto tmp13 = 0.12 + tmp12;
      const auto tmp14 = -1 * tmp13;
      const auto tmp15 = 1.0 + tmp14;
      const auto tmp16 = std::get< 0 >( u );
      const auto tmp17 = tmp16[ 1 ] / tmp15;
      const auto tmp18 = std::pow( tmp17, 2 );
      const auto tmp19 = -0.7000000000000002 * tmp8;
      const auto tmp20 = 2.7 + tmp19;
      const auto tmp21 = 2.0 + tmp20;
      const auto tmp22 = tmp21 / tmp20;
      const auto tmp23 = -1 + tmp22;
      const auto tmp24 = -1 * tmp16[ 1 ];
      const auto tmp25 = 1 + tmp24;
      const auto tmp26 = tmp25 + tmp14;
      const auto tmp27 = tmp26 / tmp15;
      const auto tmp28 = std::pow( tmp27, tmp23 );
      const auto tmp29 = -1 / tmp15;
      const auto tmp30 = tmp29 * tmp22;
      const auto tmp31 = tmp30 * tmp28;
      const auto tmp32 = -1 * tmp31;
      const auto tmp33 = tmp32 * tmp18;
      const auto tmp34 = 1 / tmp15;
      const auto tmp35 = 2 * tmp34;
      const auto tmp36 = tmp35 * tmp17;
      const auto tmp37 = std::pow( tmp27, tmp22 );
      const auto tmp38 = -1 * tmp37;
      const auto tmp39 = 1.0 + tmp38;
      const auto tmp40 = tmp39 * tmp36;
      const auto tmp41 = tmp40 + tmp33;
      const auto tmp42 = tmp41 / 0.0009;
      const auto tmp43 = 5e-11 * tmp8;
      const auto tmp44 = -5e-11 + tmp43;
      const auto tmp45 = tmp44 * tmp42;
      const auto tmp46 = std::get< 2 >( u );
      const auto tmp47 = ((tmp46[ 0 ])[ 1 ])[ 0 ] * tmp45;
      const auto tmp48 = -1 + tmp23;
      const auto tmp49 = std::pow( tmp27, tmp48 );
      const auto tmp50 = tmp23 * tmp29;
      const auto tmp51 = tmp50 * tmp49;
      const auto tmp52 = std::get< 1 >( u );
      const auto tmp53 = -1 * (tmp52[ 1 ])[ 0 ];
      const auto tmp54 = tmp53 / tmp15;
      const auto tmp55 = tmp54 * tmp22;
      const auto tmp56 = tmp55 * tmp51;
      const auto tmp57 = -1 * tmp56;
      const auto tmp58 = tmp57 * tmp18;
      const auto tmp59 = tmp55 * tmp28;
      const auto tmp60 = -1 * tmp59;
      const auto tmp61 = tmp60 * tmp36;
      const auto tmp62 = tmp61 + tmp58;
      const auto tmp63 = (tmp52[ 1 ])[ 0 ] / tmp15;
      const auto tmp64 = 2 * tmp63;
      const auto tmp65 = tmp64 * tmp17;
      const auto tmp66 = tmp32 * tmp65;
      const auto tmp67 = tmp64 * tmp34;
      const auto tmp68 = tmp39 * tmp67;
      const auto tmp69 = tmp68 + tmp66;
      const auto tmp70 = tmp69 + tmp62;
      const auto tmp71 = tmp70 / 0.0009;
      const auto tmp72 = tmp44 * tmp71;
      const auto tmp73 = 14322.6 + (tmp52[ 0 ])[ 1 ];
      const auto tmp74 = tmp73 * tmp72;
      const auto tmp75 = tmp74 + tmp47;
      const auto tmp76 = -9.993360000000001e-11 * tmp8;
      const auto tmp77 = 1e-10 + tmp76;
      const auto tmp78 = tmp77 * tmp42;
      const auto tmp79 = ((tmp46[ 0 ])[ 0 ])[ 0 ] * tmp78;
      const auto tmp80 = tmp77 * tmp71;
      const auto tmp81 = (tmp52[ 0 ])[ 0 ] * tmp80;
      const auto tmp82 = tmp81 + tmp79;
      const auto tmp83 = tmp82 + tmp75;
      const auto tmp84 = -1.0 / tmp15;
      const auto tmp85 = -1.0 / tmp20;
      const auto tmp86 = -1.0 + tmp85;
      const auto tmp87 = std::pow( tmp27, tmp86 );
      const auto tmp88 = 4245.0 * tmp8;
      const auto tmp89 = 755.0 + tmp88;
      const auto tmp90 = tmp89 * tmp85;
      const auto tmp91 = tmp90 * tmp87;
      const auto tmp92 = tmp91 * tmp84;
      const auto tmp93 = tmp45 * tmp92;
      const auto tmp94 = -1 + tmp86;
      const auto tmp95 = std::pow( tmp27, tmp94 );
      const auto tmp96 = tmp86 * tmp29;
      const auto tmp97 = tmp96 * tmp95;
      const auto tmp98 = tmp90 * tmp97;
      const auto tmp99 = tmp98 * tmp84;
      const auto tmp100 = tmp39 * tmp18;
      const auto tmp101 = tmp100 / 0.0009;
      const auto tmp102 = tmp44 * tmp101;
      const auto tmp103 = tmp102 * tmp99;
      const auto tmp104 = tmp103 + tmp93;
      const auto tmp105 = ((tmp46[ 1 ])[ 1 ])[ 0 ] * tmp104;
      const auto tmp106 = tmp72 * tmp92;
      const auto tmp107 = tmp60 * tmp18;
      const auto tmp108 = tmp39 * tmp65;
      const auto tmp109 = tmp108 + tmp107;
      const auto tmp110 = tmp109 / 0.0009;
      const auto tmp111 = tmp44 * tmp110;
      const auto tmp112 = tmp111 * tmp99;
      const auto tmp113 = tmp112 + tmp106;
      const auto tmp114 = tmp86 * tmp54;
      const auto tmp115 = tmp114 * tmp95;
      const auto tmp116 = tmp90 * tmp115;
      const auto tmp117 = tmp116 * tmp84;
      const auto tmp118 = tmp45 * tmp117;
      const auto tmp119 = -1 + tmp94;
      const auto tmp120 = std::pow( tmp27, tmp119 );
      const auto tmp121 = tmp94 * tmp29;
      const auto tmp122 = tmp121 * tmp120;
      const auto tmp123 = tmp114 * tmp122;
      const auto tmp124 = tmp123 * tmp90;
      const auto tmp125 = tmp124 * tmp84;
      const auto tmp126 = tmp102 * tmp125;
      const auto tmp127 = tmp126 + tmp118;
      const auto tmp128 = tmp127 + tmp113;
      const auto tmp129 = (tmp52[ 1 ])[ 1 ] * tmp128;
      const auto tmp130 = tmp129 + tmp105;
      const auto tmp131 = tmp78 * tmp92;
      const auto tmp132 = tmp77 * tmp101;
      const auto tmp133 = tmp132 * tmp99;
      const auto tmp134 = tmp133 + tmp131;
      const auto tmp135 = ((tmp46[ 1 ])[ 0 ])[ 0 ] * tmp134;
      const auto tmp136 = tmp80 * tmp92;
      const auto tmp137 = tmp77 * tmp110;
      const auto tmp138 = tmp137 * tmp99;
      const auto tmp139 = tmp138 + tmp136;
      const auto tmp140 = tmp78 * tmp117;
      const auto tmp141 = tmp132 * tmp125;
      const auto tmp142 = tmp141 + tmp140;
      const auto tmp143 = tmp142 + tmp139;
      const auto tmp144 = (tmp52[ 1 ])[ 0 ] * tmp143;
      const auto tmp145 = tmp144 + tmp135;
      const auto tmp146 = tmp145 + tmp130;
      const auto tmp147 = tmp146 + tmp83;
      const auto tmp148 = ((tmp46[ 0 ])[ 0 ])[ 1 ] * tmp45;
      const auto tmp149 = -1 * (tmp52[ 1 ])[ 1 ];
      const auto tmp150 = tmp149 / tmp15;
      const auto tmp151 = tmp150 * tmp22;
      const auto tmp152 = tmp151 * tmp51;
      const auto tmp153 = -1 * tmp152;
      const auto tmp154 = tmp153 * tmp18;
      const auto tmp155 = tmp151 * tmp28;
      const auto tmp156 = -1 * tmp155;
      const auto tmp157 = tmp156 * tmp36;
      const auto tmp158 = tmp157 + tmp154;
      const auto tmp159 = (tmp52[ 1 ])[ 1 ] / tmp15;
      const auto tmp160 = 2 * tmp159;
      const auto tmp161 = tmp160 * tmp17;
      const auto tmp162 = tmp32 * tmp161;
      const auto tmp163 = tmp160 * tmp34;
      const auto tmp164 = tmp39 * tmp163;
      const auto tmp165 = tmp164 + tmp162;
      const auto tmp166 = tmp165 + tmp158;
      const auto tmp167 = tmp166 / 0.0009;
      const auto tmp168 = tmp44 * tmp167;
      const auto tmp169 = (tmp52[ 0 ])[ 0 ] * tmp168;
      const auto tmp170 = tmp169 + tmp148;
      const auto tmp171 = ((tmp46[ 0 ])[ 1 ])[ 1 ] * tmp78;
      const auto tmp172 = tmp77 * tmp167;
      const auto tmp173 = tmp73 * tmp172;
      const auto tmp174 = tmp173 + tmp171;
      const auto tmp175 = tmp174 + tmp170;
      const auto tmp176 = ((tmp46[ 1 ])[ 0 ])[ 1 ] * tmp104;
      const auto tmp177 = tmp168 * tmp92;
      const auto tmp178 = tmp156 * tmp18;
      const auto tmp179 = tmp39 * tmp161;
      const auto tmp180 = tmp179 + tmp178;
      const auto tmp181 = tmp180 / 0.0009;
      const auto tmp182 = tmp44 * tmp181;
      const auto tmp183 = tmp182 * tmp99;
      const auto tmp184 = tmp183 + tmp177;
      const auto tmp185 = tmp86 * tmp150;
      const auto tmp186 = tmp185 * tmp95;
      const auto tmp187 = tmp90 * tmp186;
      const auto tmp188 = tmp187 * tmp84;
      const auto tmp189 = tmp45 * tmp188;
      const auto tmp190 = tmp185 * tmp122;
      const auto tmp191 = tmp190 * tmp90;
      const auto tmp192 = tmp191 * tmp84;
      const auto tmp193 = tmp102 * tmp192;
      const auto tmp194 = tmp193 + tmp189;
      const auto tmp195 = tmp194 + tmp184;
      const auto tmp196 = (tmp52[ 1 ])[ 0 ] * tmp195;
      const auto tmp197 = tmp196 + tmp176;
      const auto tmp198 = ((tmp46[ 1 ])[ 1 ])[ 1 ] * tmp134;
      const auto tmp199 = tmp172 * tmp92;
      const auto tmp200 = tmp77 * tmp181;
      const auto tmp201 = tmp200 * tmp99;
      const auto tmp202 = tmp201 + tmp199;
      const auto tmp203 = tmp78 * tmp188;
      const auto tmp204 = tmp132 * tmp192;
      const auto tmp205 = tmp204 + tmp203;
      const auto tmp206 = tmp205 + tmp202;
      const auto tmp207 = (tmp52[ 1 ])[ 1 ] * tmp206;
      const auto tmp208 = tmp207 + tmp198;
      const auto tmp209 = tmp208 + tmp197;
      const auto tmp210 = tmp209 + tmp175;
      const auto tmp211 = tmp210 + tmp147;
      const auto tmp212 = -1 * tmp211;
      const auto tmp213 = tmp212 + tmp11;
      const auto tmp214 = 2 * tmp213;
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp215 = evaluateCoefficient< 0 >( x );
      const auto tmp216 = -1 * tmp215[ 1 ];
      const auto tmp217 = tmp16[ 1 ] + tmp216;
      const auto tmp218 = tmp217 * tmp10;
      const auto tmp219 = tmp218 / tmp2;
      const auto tmp220 = ((tmp46[ 0 ])[ 1 ])[ 0 ] * tmp102;
      const auto tmp221 = tmp73 * tmp111;
      const auto tmp222 = tmp221 + tmp220;
      const auto tmp223 = ((tmp46[ 0 ])[ 0 ])[ 0 ] * tmp132;
      const auto tmp224 = (tmp52[ 0 ])[ 0 ] * tmp137;
      const auto tmp225 = tmp224 + tmp223;
      const auto tmp226 = tmp225 + tmp222;
      const auto tmp227 = tmp102 * tmp92;
      const auto tmp228 = ((tmp46[ 1 ])[ 1 ])[ 0 ] * tmp227;
      const auto tmp229 = tmp111 * tmp92;
      const auto tmp230 = tmp102 * tmp117;
      const auto tmp231 = tmp230 + tmp229;
      const auto tmp232 = (tmp52[ 1 ])[ 1 ] * tmp231;
      const auto tmp233 = tmp232 + tmp228;
      const auto tmp234 = tmp132 * tmp92;
      const auto tmp235 = ((tmp46[ 1 ])[ 0 ])[ 0 ] * tmp234;
      const auto tmp236 = tmp137 * tmp92;
      const auto tmp237 = tmp132 * tmp117;
      const auto tmp238 = tmp237 + tmp236;
      const auto tmp239 = (tmp52[ 1 ])[ 0 ] * tmp238;
      const auto tmp240 = tmp239 + tmp235;
      const auto tmp241 = tmp240 + tmp233;
      const auto tmp242 = tmp241 + tmp226;
      const auto tmp243 = ((tmp46[ 0 ])[ 0 ])[ 1 ] * tmp102;
      const auto tmp244 = (tmp52[ 0 ])[ 0 ] * tmp182;
      const auto tmp245 = tmp244 + tmp243;
      const auto tmp246 = ((tmp46[ 0 ])[ 1 ])[ 1 ] * tmp132;
      const auto tmp247 = tmp73 * tmp200;
      const auto tmp248 = tmp247 + tmp246;
      const auto tmp249 = tmp248 + tmp245;
      const auto tmp250 = ((tmp46[ 1 ])[ 0 ])[ 1 ] * tmp227;
      const auto tmp251 = tmp182 * tmp92;
      const auto tmp252 = tmp102 * tmp188;
      const auto tmp253 = tmp252 + tmp251;
      const auto tmp254 = (tmp52[ 1 ])[ 0 ] * tmp253;
      const auto tmp255 = tmp254 + tmp250;
      const auto tmp256 = ((tmp46[ 1 ])[ 1 ])[ 1 ] * tmp234;
      const auto tmp257 = tmp200 * tmp92;
      const auto tmp258 = tmp132 * tmp188;
      const auto tmp259 = tmp258 + tmp257;
      const auto tmp260 = (tmp52[ 1 ])[ 1 ] * tmp259;
      const auto tmp261 = tmp260 + tmp256;
      const auto tmp262 = tmp261 + tmp255;
      const auto tmp263 = tmp262 + tmp249;
      const auto tmp264 = tmp263 + tmp242;
      const auto tmp265 = -1 * tmp264;
      const auto tmp266 = tmp265 + tmp219;
      const auto tmp267 = tmp266 * tmp214;
      const auto tmp268 = tmp267 * tmp1;
      const auto tmp269 = tmp137 + tmp182;
      const auto tmp270 = -1 * tmp269;
      const auto tmp271 = 2 * tmp270;
      const auto tmp272 = tmp266 * tmp271;
      const auto tmp273 = tmp272 * tmp1;
      const auto tmp274 = tmp111 + tmp200;
      const auto tmp275 = -1 * tmp274;
      const auto tmp276 = 2 * tmp275;
      const auto tmp277 = tmp266 * tmp276;
      const auto tmp278 = tmp277 * tmp1;
      const auto tmp279 = tmp73 * tmp45;
      const auto tmp280 = (tmp52[ 0 ])[ 0 ] * tmp78;
      const auto tmp281 = tmp280 + tmp279;
      const auto tmp282 = (tmp52[ 1 ])[ 1 ] * tmp104;
      const auto tmp283 = (tmp52[ 1 ])[ 0 ] * tmp134;
      const auto tmp284 = tmp238 + tmp283;
      const auto tmp285 = tmp284 + tmp282;
      const auto tmp286 = tmp285 + tmp281;
      const auto tmp287 = tmp286 + tmp253;
      const auto tmp288 = -1 * tmp287;
      const auto tmp289 = 2 * tmp288;
      const auto tmp290 = tmp266 * tmp289;
      const auto tmp291 = tmp290 * tmp1;
      const auto tmp292 = (tmp52[ 0 ])[ 0 ] * tmp45;
      const auto tmp293 = tmp73 * tmp78;
      const auto tmp294 = tmp293 + tmp292;
      const auto tmp295 = (tmp52[ 1 ])[ 0 ] * tmp104;
      const auto tmp296 = (tmp52[ 1 ])[ 1 ] * tmp134;
      const auto tmp297 = tmp259 + tmp296;
      const auto tmp298 = tmp297 + tmp295;
      const auto tmp299 = tmp298 + tmp294;
      const auto tmp300 = tmp299 + tmp231;
      const auto tmp301 = -1 * tmp300;
      const auto tmp302 = 2 * tmp301;
      const auto tmp303 = tmp266 * tmp302;
      const auto tmp304 = tmp303 * tmp1;
      const auto tmp305 = -1 * tmp132;
      const auto tmp306 = 2 * tmp305;
      const auto tmp307 = tmp266 * tmp306;
      const auto tmp308 = tmp307 * tmp1;
      const auto tmp309 = -1 * tmp102;
      const auto tmp310 = 2 * tmp309;
      const auto tmp311 = tmp266 * tmp310;
      const auto tmp312 = tmp311 * tmp1;
      const auto tmp313 = -1 * tmp234;
      const auto tmp314 = 2 * tmp313;
      const auto tmp315 = tmp266 * tmp314;
      const auto tmp316 = tmp315 * tmp1;
      const auto tmp317 = -1 * tmp227;
      const auto tmp318 = 2 * tmp317;
      const auto tmp319 = tmp266 * tmp318;
      const auto tmp320 = tmp319 * tmp1;
      return [ tmp268, tmp273, tmp278, tmp291, tmp304, tmp308, tmp312, tmp316, tmp320 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (((((((((((tmp268 * (std::get< 0 >( phi ))[ 1 ] + tmp273 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp278 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ]) + tmp291 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ]) + tmp304 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ]) + tmp308 * (((std::get< 2 >( phi ))[ 0 ])[ 0 ])[ 0 ]) + tmp312 * (((std::get< 2 >( phi ))[ 0 ])[ 0 ])[ 1 ]) + tmp312 * (((std::get< 2 >( phi ))[ 0 ])[ 1 ])[ 0 ]) + tmp308 * (((std::get< 2 >( phi ))[ 0 ])[ 1 ])[ 1 ]) + tmp316 * (((std::get< 2 >( phi ))[ 1 ])[ 0 ])[ 0 ]) + tmp320 * (((std::get< 2 >( phi ))[ 1 ])[ 0 ])[ 1 ]) + tmp320 * (((std::get< 2 >( phi ))[ 1 ])[ 1 ])[ 0 ]) + tmp316 * (((std::get< 2 >( phi ))[ 1 ])[ 1 ])[ 1 ] } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = -0.65 + tmp0[ 1 ];
      const auto tmp2 = std::abs( tmp1 );
      const auto tmp3 = -0.45 + tmp0[ 0 ];
      const auto tmp4 = std::abs( tmp3 );
      const auto tmp5 = (tmp4 < 0.06 ? 1.0 : 0.0) * (tmp2 < 1e-08 ? 1.0 : 0.0);
      const auto tmp6 = -0.49 + tmp0[ 1 ];
      const auto tmp7 = std::abs( tmp6 );
      const auto tmp8 = (tmp7 < 0.03 ? 1.0 : 0.0) * (tmp4 < 0.11 ? 1.0 : 0.0);
      const auto tmp9 = -0.01999999999999999 * tmp8;
      const auto tmp10 = 0.12 + tmp9;
      const auto tmp11 = -1 * tmp10;
      const auto tmp12 = 1.0 + tmp11;
      const auto tmp13 = std::get< 0 >( u );
      const auto tmp14 = tmp13[ 1 ] / tmp12;
      const auto tmp15 = std::pow( tmp14, 2 );
      const auto tmp16 = -0.7000000000000002 * tmp8;
      const auto tmp17 = 2.7 + tmp16;
      const auto tmp18 = 2.0 + tmp17;
      const auto tmp19 = tmp18 / tmp17;
      const auto tmp20 = -1 * tmp13[ 1 ];
      const auto tmp21 = 1 + tmp20;
      const auto tmp22 = tmp21 + tmp11;
      const auto tmp23 = tmp22 / tmp12;
      const auto tmp24 = std::pow( tmp23, tmp19 );
      const auto tmp25 = -1 * tmp24;
      const auto tmp26 = 1.0 + tmp25;
      const auto tmp27 = tmp26 * tmp15;
      const auto tmp28 = tmp27 / 0.0009;
      const auto tmp29 = 5e-11 * tmp8;
      const auto tmp30 = -5e-11 + tmp29;
      const auto tmp31 = tmp30 * tmp28;
      const auto tmp32 = std::get< 1 >( u );
      const auto tmp33 = 14322.6 + (tmp32[ 0 ])[ 1 ];
      const auto tmp34 = tmp33 * tmp31;
      const auto tmp35 = -9.993360000000001e-11 * tmp8;
      const auto tmp36 = 1e-10 + tmp35;
      const auto tmp37 = tmp36 * tmp28;
      const auto tmp38 = (tmp32[ 0 ])[ 0 ] * tmp37;
      const auto tmp39 = tmp38 + tmp34;
      const auto tmp40 = -1.0 / tmp12;
      const auto tmp41 = -1.0 / tmp17;
      const auto tmp42 = -1.0 + tmp41;
      const auto tmp43 = std::pow( tmp23, tmp42 );
      const auto tmp44 = 4245.0 * tmp8;
      const auto tmp45 = 755.0 + tmp44;
      const auto tmp46 = tmp45 * tmp41;
      const auto tmp47 = tmp46 * tmp43;
      const auto tmp48 = tmp47 * tmp40;
      const auto tmp49 = tmp31 * tmp48;
      const auto tmp50 = (tmp32[ 1 ])[ 1 ] * tmp49;
      const auto tmp51 = tmp37 * tmp48;
      const auto tmp52 = (tmp32[ 1 ])[ 0 ] * tmp51;
      const auto tmp53 = tmp52 + tmp50;
      const auto tmp54 = tmp53 + tmp39;
      GlobalCoordinateType tmp55 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp56 = tmp55[ 0 ] * tmp54;
      const auto tmp57 = (tmp32[ 0 ])[ 0 ] * tmp31;
      const auto tmp58 = tmp33 * tmp37;
      const auto tmp59 = tmp58 + tmp57;
      const auto tmp60 = (tmp32[ 1 ])[ 0 ] * tmp49;
      const auto tmp61 = (tmp32[ 1 ])[ 1 ] * tmp51;
      const auto tmp62 = tmp61 + tmp60;
      const auto tmp63 = tmp62 + tmp59;
      const auto tmp64 = tmp55[ 1 ] * tmp63;
      const auto tmp65 = tmp64 + tmp56;
      const auto tmp66 = -5.137e-05 + tmp65;
      const auto tmp67 = std::pow( tmp66, 2 );
      auto tmp68 = facetArea_;
      auto tmp69 = cellVolume();
      const auto tmp70 = tmp69 / tmp68;
      const auto tmp71 = tmp70 * tmp67;
      const auto tmp72 = tmp71 * tmp5;
      const auto tmp73 = std::pow( tmp20, 2 );
      const auto tmp74 = 0.5 / tmp12;
      const auto tmp75 = std::pow( tmp74, 2 );
      const auto tmp76 = 0.5 + tmp11;
      const auto tmp77 = tmp76 / tmp12;
      const auto tmp78 = std::pow( tmp77, tmp19 );
      const auto tmp79 = -1 * tmp78;
      const auto tmp80 = 1.0 + tmp79;
      const auto tmp81 = tmp80 * tmp75;
      const auto tmp82 = tmp81 / 0.0009;
      const auto tmp83 = std::pow( tmp77, tmp42 );
      const auto tmp84 = tmp46 * tmp83;
      const auto tmp85 = tmp84 * tmp40;
      const auto tmp86 = tmp85 * tmp82;
      const auto tmp87 = tmp55[ 1 ] * tmp30;
      const auto tmp88 = tmp55[ 0 ] * tmp36;
      const auto tmp89 = tmp88 + tmp87;
      const auto tmp90 = tmp55[ 0 ] * tmp89;
      const auto tmp91 = tmp55[ 0 ] * tmp30;
      const auto tmp92 = tmp55[ 1 ] * tmp36;
      const auto tmp93 = tmp92 + tmp91;
      const auto tmp94 = tmp55[ 1 ] * tmp93;
      const auto tmp95 = tmp94 + tmp90;
      double tmp96 = constant< 1 >();
      const auto tmp97 = tmp96 * tmp95;
      const auto tmp98 = tmp97 * tmp86;
      const auto tmp99 = std::pow( tmp98, 2 );
      const auto tmp100 = tmp99 / tmp70;
      const auto tmp101 = tmp100 * tmp73;
      const auto tmp102 = -0.9 + tmp0[ 0 ];
      const auto tmp103 = std::abs( tmp102 );
      const auto tmp104 = std::abs( tmp0[ 0 ] );
      const auto tmp105 = (tmp104 < 1e-08 ? 1.0 : 0.0) + (tmp103 < 1e-08 ? 1.0 : 0.0);
      const auto tmp106 = tmp105 * tmp101;
      const auto tmp107 = tmp106 + tmp72;
      return RangeValueType{ { tmp107 } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      using std::pow;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = -0.65 + tmp0[ 1 ];
      const auto tmp2 = std::abs( tmp1 );
      const auto tmp3 = -0.45 + tmp0[ 0 ];
      const auto tmp4 = std::abs( tmp3 );
      const auto tmp5 = (tmp4 < 0.06 ? 1.0 : 0.0) * (tmp2 < 1e-08 ? 1.0 : 0.0);
      auto tmp6 = facetArea_;
      auto tmp7 = cellVolume();
      const auto tmp8 = tmp7 / tmp6;
      const auto tmp9 = -0.49 + tmp0[ 1 ];
      const auto tmp10 = std::abs( tmp9 );
      const auto tmp11 = (tmp10 < 0.03 ? 1.0 : 0.0) * (tmp4 < 0.11 ? 1.0 : 0.0);
      const auto tmp12 = -0.01999999999999999 * tmp11;
      const auto tmp13 = 0.12 + tmp12;
      const auto tmp14 = -1 * tmp13;
      const auto tmp15 = 1.0 + tmp14;
      const auto tmp16 = std::get< 0 >( u );
      const auto tmp17 = tmp16[ 1 ] / tmp15;
      const auto tmp18 = std::pow( tmp17, 2 );
      const auto tmp19 = -0.7000000000000002 * tmp11;
      const auto tmp20 = 2.7 + tmp19;
      const auto tmp21 = 2.0 + tmp20;
      const auto tmp22 = tmp21 / tmp20;
      const auto tmp23 = -1 + tmp22;
      const auto tmp24 = -1 * tmp16[ 1 ];
      const auto tmp25 = 1 + tmp24;
      const auto tmp26 = tmp25 + tmp14;
      const auto tmp27 = tmp26 / tmp15;
      const auto tmp28 = std::pow( tmp27, tmp23 );
      const auto tmp29 = -1 / tmp15;
      const auto tmp30 = tmp29 * tmp22;
      const auto tmp31 = tmp30 * tmp28;
      const auto tmp32 = -1 * tmp31;
      const auto tmp33 = tmp32 * tmp18;
      const auto tmp34 = 1 / tmp15;
      const auto tmp35 = 2 * tmp34;
      const auto tmp36 = tmp35 * tmp17;
      const auto tmp37 = std::pow( tmp27, tmp22 );
      const auto tmp38 = -1 * tmp37;
      const auto tmp39 = 1.0 + tmp38;
      const auto tmp40 = tmp39 * tmp36;
      const auto tmp41 = tmp40 + tmp33;
      const auto tmp42 = tmp41 / 0.0009;
      const auto tmp43 = 5e-11 * tmp11;
      const auto tmp44 = -5e-11 + tmp43;
      const auto tmp45 = tmp44 * tmp42;
      const auto tmp46 = std::get< 1 >( u );
      const auto tmp47 = 14322.6 + (tmp46[ 0 ])[ 1 ];
      const auto tmp48 = tmp47 * tmp45;
      const auto tmp49 = -9.993360000000001e-11 * tmp11;
      const auto tmp50 = 1e-10 + tmp49;
      const auto tmp51 = tmp50 * tmp42;
      const auto tmp52 = (tmp46[ 0 ])[ 0 ] * tmp51;
      const auto tmp53 = tmp52 + tmp48;
      const auto tmp54 = -1.0 / tmp15;
      const auto tmp55 = -1.0 / tmp20;
      const auto tmp56 = -1.0 + tmp55;
      const auto tmp57 = std::pow( tmp27, tmp56 );
      const auto tmp58 = 4245.0 * tmp11;
      const auto tmp59 = 755.0 + tmp58;
      const auto tmp60 = tmp59 * tmp55;
      const auto tmp61 = tmp60 * tmp57;
      const auto tmp62 = tmp61 * tmp54;
      const auto tmp63 = tmp45 * tmp62;
      const auto tmp64 = -1 + tmp56;
      const auto tmp65 = std::pow( tmp27, tmp64 );
      const auto tmp66 = tmp56 * tmp29;
      const auto tmp67 = tmp66 * tmp65;
      const auto tmp68 = tmp60 * tmp67;
      const auto tmp69 = tmp68 * tmp54;
      const auto tmp70 = tmp39 * tmp18;
      const auto tmp71 = tmp70 / 0.0009;
      const auto tmp72 = tmp44 * tmp71;
      const auto tmp73 = tmp72 * tmp69;
      const auto tmp74 = tmp73 + tmp63;
      const auto tmp75 = (tmp46[ 1 ])[ 1 ] * tmp74;
      const auto tmp76 = tmp51 * tmp62;
      const auto tmp77 = tmp50 * tmp71;
      const auto tmp78 = tmp77 * tmp69;
      const auto tmp79 = tmp78 + tmp76;
      const auto tmp80 = (tmp46[ 1 ])[ 0 ] * tmp79;
      const auto tmp81 = tmp80 + tmp75;
      const auto tmp82 = tmp81 + tmp53;
      GlobalCoordinateType tmp83 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp84 = tmp83[ 0 ] * tmp82;
      const auto tmp85 = (tmp46[ 0 ])[ 0 ] * tmp45;
      const auto tmp86 = tmp47 * tmp51;
      const auto tmp87 = tmp86 + tmp85;
      const auto tmp88 = (tmp46[ 1 ])[ 0 ] * tmp74;
      const auto tmp89 = (tmp46[ 1 ])[ 1 ] * tmp79;
      const auto tmp90 = tmp89 + tmp88;
      const auto tmp91 = tmp90 + tmp87;
      const auto tmp92 = tmp83[ 1 ] * tmp91;
      const auto tmp93 = tmp92 + tmp84;
      const auto tmp94 = 2 * tmp93;
      const auto tmp95 = tmp47 * tmp72;
      const auto tmp96 = (tmp46[ 0 ])[ 0 ] * tmp77;
      const auto tmp97 = tmp96 + tmp95;
      const auto tmp98 = tmp72 * tmp62;
      const auto tmp99 = (tmp46[ 1 ])[ 1 ] * tmp98;
      const auto tmp100 = tmp77 * tmp62;
      const auto tmp101 = (tmp46[ 1 ])[ 0 ] * tmp100;
      const auto tmp102 = tmp101 + tmp99;
      const auto tmp103 = tmp102 + tmp97;
      const auto tmp104 = tmp83[ 0 ] * tmp103;
      const auto tmp105 = (tmp46[ 0 ])[ 0 ] * tmp72;
      const auto tmp106 = tmp47 * tmp77;
      const auto tmp107 = tmp106 + tmp105;
      const auto tmp108 = (tmp46[ 1 ])[ 0 ] * tmp98;
      const auto tmp109 = (tmp46[ 1 ])[ 1 ] * tmp100;
      const auto tmp110 = tmp109 + tmp108;
      const auto tmp111 = tmp110 + tmp107;
      const auto tmp112 = tmp83[ 1 ] * tmp111;
      const auto tmp113 = tmp112 + tmp104;
      const auto tmp114 = -5.137e-05 + tmp113;
      const auto tmp115 = tmp114 * tmp94;
      const auto tmp116 = tmp115 * tmp8;
      const auto tmp117 = tmp116 * tmp5;
      const auto tmp118 = 0.5 / tmp15;
      const auto tmp119 = std::pow( tmp118, 2 );
      const auto tmp120 = 0.5 + tmp14;
      const auto tmp121 = tmp120 / tmp15;
      const auto tmp122 = std::pow( tmp121, tmp22 );
      const auto tmp123 = -1 * tmp122;
      const auto tmp124 = 1.0 + tmp123;
      const auto tmp125 = tmp124 * tmp119;
      const auto tmp126 = tmp125 / 0.0009;
      const auto tmp127 = std::pow( tmp121, tmp56 );
      const auto tmp128 = tmp60 * tmp127;
      const auto tmp129 = tmp128 * tmp54;
      const auto tmp130 = tmp129 * tmp126;
      const auto tmp131 = tmp83[ 1 ] * tmp44;
      const auto tmp132 = tmp83[ 0 ] * tmp50;
      const auto tmp133 = tmp132 + tmp131;
      const auto tmp134 = tmp83[ 0 ] * tmp133;
      const auto tmp135 = tmp83[ 0 ] * tmp44;
      const auto tmp136 = tmp83[ 1 ] * tmp50;
      const auto tmp137 = tmp136 + tmp135;
      const auto tmp138 = tmp83[ 1 ] * tmp137;
      const auto tmp139 = tmp138 + tmp134;
      double tmp140 = constant< 1 >();
      const auto tmp141 = tmp140 * tmp139;
      const auto tmp142 = tmp141 * tmp130;
      const auto tmp143 = std::pow( tmp142, 2 );
      const auto tmp144 = tmp143 / tmp8;
      const auto tmp145 = -2 * tmp24;
      const auto tmp146 = tmp145 * tmp144;
      const auto tmp147 = -0.9 + tmp0[ 0 ];
      const auto tmp148 = std::abs( tmp147 );
      const auto tmp149 = std::abs( tmp0[ 0 ] );
      const auto tmp150 = (tmp149 < 1e-08 ? 1.0 : 0.0) + (tmp148 < 1e-08 ? 1.0 : 0.0);
      const auto tmp151 = tmp150 * tmp146;
      const auto tmp152 = tmp151 + tmp117;
      const auto tmp153 = tmp83[ 1 ] * tmp72;
      const auto tmp154 = tmp83[ 0 ] * tmp77;
      const auto tmp155 = tmp154 + tmp153;
      const auto tmp156 = 2 * tmp155;
      const auto tmp157 = tmp114 * tmp156;
      const auto tmp158 = tmp157 * tmp8;
      const auto tmp159 = tmp158 * tmp5;
      const auto tmp160 = tmp83[ 0 ] * tmp72;
      const auto tmp161 = tmp83[ 1 ] * tmp77;
      const auto tmp162 = tmp161 + tmp160;
      const auto tmp163 = 2 * tmp162;
      const auto tmp164 = tmp114 * tmp163;
      const auto tmp165 = tmp164 * tmp8;
      const auto tmp166 = tmp165 * tmp5;
      const auto tmp167 = tmp83[ 1 ] * tmp98;
      const auto tmp168 = tmp83[ 0 ] * tmp100;
      const auto tmp169 = tmp168 + tmp167;
      const auto tmp170 = 2 * tmp169;
      const auto tmp171 = tmp114 * tmp170;
      const auto tmp172 = tmp171 * tmp8;
      const auto tmp173 = tmp172 * tmp5;
      const auto tmp174 = tmp83[ 0 ] * tmp98;
      const auto tmp175 = tmp83[ 1 ] * tmp100;
      const auto tmp176 = tmp175 + tmp174;
      const auto tmp177 = 2 * tmp176;
      const auto tmp178 = tmp114 * tmp177;
      const auto tmp179 = tmp178 * tmp8;
      const auto tmp180 = tmp179 * tmp5;
      return [ tmp152, tmp159, tmp166, tmp173, tmp180 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (((tmp152 * (std::get< 0 >( phi ))[ 1 ] + tmp159 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp166 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ]) + tmp173 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ]) + tmp180 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ] } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      using std::max;
      using std::pow;
      GlobalCoordinateType tmp0 = intersection_.unitOuterNormal( xIn.localPosition() );
      GlobalCoordinateType tmp1 = geometry().global( Dune::Fem::coordinate( xIn ) );
      const auto tmp2 = -0.45 + tmp1[ 0 ];
      const auto tmp3 = std::abs( tmp2 );
      const auto tmp4 = -0.49 + tmp1[ 1 ];
      const auto tmp5 = std::abs( tmp4 );
      const auto tmp6 = (tmp5 < 0.03 ? 1.0 : 0.0) * (tmp3 < 0.11 ? 1.0 : 0.0);
      const auto tmp7 = -0.01999999999999999 * tmp6;
      const auto tmp8 = 0.12 + tmp7;
      const auto tmp9 = -1 * tmp8;
      const auto tmp10 = 1.0 + tmp9;
      const auto tmp11 = std::get< 0 >( uOut );
      const auto tmp12 = tmp11[ 1 ] / tmp10;
      const auto tmp13 = std::pow( tmp12, 2 );
      const auto tmp14 = -0.7000000000000002 * tmp6;
      const auto tmp15 = 2.7 + tmp14;
      const auto tmp16 = 2.0 + tmp15;
      const auto tmp17 = tmp16 / tmp15;
      const auto tmp18 = -1 * tmp11[ 1 ];
      const auto tmp19 = 1 + tmp18;
      const auto tmp20 = tmp19 + tmp9;
      const auto tmp21 = tmp20 / tmp10;
      const auto tmp22 = std::pow( tmp21, tmp17 );
      const auto tmp23 = -1 * tmp22;
      const auto tmp24 = 1.0 + tmp23;
      const auto tmp25 = tmp24 * tmp13;
      const auto tmp26 = tmp25 / 0.0009;
      const auto tmp27 = 5e-11 * tmp6;
      const auto tmp28 = -5e-11 + tmp27;
      const auto tmp29 = tmp28 * tmp26;
      const auto tmp30 = std::get< 1 >( uOut );
      const auto tmp31 = (tmp30[ 0 ])[ 0 ] * tmp29;
      const auto tmp32 = -9.993360000000001e-11 * tmp6;
      const auto tmp33 = 1e-10 + tmp32;
      const auto tmp34 = tmp33 * tmp26;
      const auto tmp35 = 14322.6 + (tmp30[ 0 ])[ 1 ];
      const auto tmp36 = tmp35 * tmp34;
      const auto tmp37 = tmp36 + tmp31;
      const auto tmp38 = -1.0 / tmp10;
      const auto tmp39 = -1.0 / tmp15;
      const auto tmp40 = -1.0 + tmp39;
      const auto tmp41 = std::pow( tmp21, tmp40 );
      const auto tmp42 = 4245.0 * tmp6;
      const auto tmp43 = 755.0 + tmp42;
      const auto tmp44 = tmp43 * tmp39;
      const auto tmp45 = tmp44 * tmp41;
      const auto tmp46 = tmp45 * tmp38;
      const auto tmp47 = tmp29 * tmp46;
      const auto tmp48 = (tmp30[ 1 ])[ 0 ] * tmp47;
      const auto tmp49 = tmp34 * tmp46;
      const auto tmp50 = (tmp30[ 1 ])[ 1 ] * tmp49;
      const auto tmp51 = tmp50 + tmp48;
      const auto tmp52 = tmp51 + tmp37;
      const auto tmp53 = -1 * tmp52;
      const auto tmp54 = std::get< 0 >( uIn );
      const auto tmp55 = tmp54[ 1 ] / tmp10;
      const auto tmp56 = std::pow( tmp55, 2 );
      const auto tmp57 = -1 * tmp54[ 1 ];
      const auto tmp58 = 1 + tmp57;
      const auto tmp59 = tmp58 + tmp9;
      const auto tmp60 = tmp59 / tmp10;
      const auto tmp61 = std::pow( tmp60, tmp17 );
      const auto tmp62 = -1 * tmp61;
      const auto tmp63 = 1.0 + tmp62;
      const auto tmp64 = tmp63 * tmp56;
      const auto tmp65 = tmp64 / 0.0009;
      const auto tmp66 = tmp28 * tmp65;
      const auto tmp67 = std::get< 1 >( uIn );
      const auto tmp68 = (tmp67[ 0 ])[ 0 ] * tmp66;
      const auto tmp69 = tmp33 * tmp65;
      const auto tmp70 = 14322.6 + (tmp67[ 0 ])[ 1 ];
      const auto tmp71 = tmp70 * tmp69;
      const auto tmp72 = tmp71 + tmp68;
      const auto tmp73 = std::pow( tmp60, tmp40 );
      const auto tmp74 = tmp44 * tmp73;
      const auto tmp75 = tmp74 * tmp38;
      const auto tmp76 = tmp66 * tmp75;
      const auto tmp77 = (tmp67[ 1 ])[ 0 ] * tmp76;
      const auto tmp78 = tmp69 * tmp75;
      const auto tmp79 = (tmp67[ 1 ])[ 1 ] * tmp78;
      const auto tmp80 = tmp79 + tmp77;
      const auto tmp81 = tmp80 + tmp72;
      const auto tmp82 = tmp81 + tmp53;
      const auto tmp83 = tmp82 * tmp0[ 1 ];
      const auto tmp84 = tmp35 * tmp29;
      const auto tmp85 = (tmp30[ 0 ])[ 0 ] * tmp34;
      const auto tmp86 = tmp85 + tmp84;
      const auto tmp87 = (tmp30[ 1 ])[ 1 ] * tmp47;
      const auto tmp88 = (tmp30[ 1 ])[ 0 ] * tmp49;
      const auto tmp89 = tmp88 + tmp87;
      const auto tmp90 = tmp89 + tmp86;
      const auto tmp91 = -1 * tmp90;
      const auto tmp92 = tmp70 * tmp66;
      const auto tmp93 = (tmp67[ 0 ])[ 0 ] * tmp69;
      const auto tmp94 = tmp93 + tmp92;
      const auto tmp95 = (tmp67[ 1 ])[ 1 ] * tmp76;
      const auto tmp96 = (tmp67[ 1 ])[ 0 ] * tmp78;
      const auto tmp97 = tmp96 + tmp95;
      const auto tmp98 = tmp97 + tmp94;
      const auto tmp99 = tmp98 + tmp91;
      const auto tmp100 = tmp99 * tmp0[ 0 ];
      const auto tmp101 = tmp100 + tmp83;
      const auto tmp102 = std::pow( tmp101, 2 );
      auto tmp103 = facetArea_;
      auto tmp104 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
      auto tmp105 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
      const auto tmp106 = tmp105 + tmp104;
      const auto tmp107 = 0.5 * tmp106;
      const auto tmp108 = tmp107 / tmp103;
      const auto tmp109 = tmp108 * tmp102;
      const auto tmp110 = 0.5 * tmp109;
      const auto tmp111 = tmp54[ 1 ] + tmp18;
      const auto tmp112 = std::pow( tmp111, 2 );
      const auto tmp113 = 0.5 / tmp10;
      const auto tmp114 = std::pow( tmp113, 2 );
      const auto tmp115 = 0.5 + tmp9;
      const auto tmp116 = tmp115 / tmp10;
      const auto tmp117 = std::pow( tmp116, tmp17 );
      const auto tmp118 = -1 * tmp117;
      const auto tmp119 = 1.0 + tmp118;
      const auto tmp120 = tmp119 * tmp114;
      const auto tmp121 = tmp120 / 0.0009;
      const auto tmp122 = std::pow( tmp116, tmp40 );
      const auto tmp123 = tmp44 * tmp122;
      const auto tmp124 = tmp123 * tmp38;
      const auto tmp125 = tmp124 * tmp121;
      const auto tmp126 = std::max( tmp125, tmp125 );
      const auto tmp127 = -1 * tmp0[ 1 ];
      const auto tmp128 = tmp33 * tmp127;
      const auto tmp129 = -1 * tmp0[ 0 ];
      const auto tmp130 = tmp28 * tmp129;
      const auto tmp131 = tmp130 + tmp128;
      const auto tmp132 = tmp131 * tmp127;
      const auto tmp133 = tmp28 * tmp127;
      const auto tmp134 = tmp33 * tmp129;
      const auto tmp135 = tmp134 + tmp133;
      const auto tmp136 = tmp135 * tmp129;
      const auto tmp137 = tmp136 + tmp132;
      const auto tmp138 = tmp0[ 1 ] * tmp28;
      const auto tmp139 = tmp0[ 0 ] * tmp33;
      const auto tmp140 = tmp139 + tmp138;
      const auto tmp141 = tmp0[ 0 ] * tmp140;
      const auto tmp142 = tmp0[ 0 ] * tmp28;
      const auto tmp143 = tmp0[ 1 ] * tmp33;
      const auto tmp144 = tmp143 + tmp142;
      const auto tmp145 = tmp0[ 1 ] * tmp144;
      const auto tmp146 = tmp145 + tmp141;
      const auto tmp147 = tmp146 + tmp137;
      const auto tmp148 = 0.5 * tmp147;
      const auto tmp149 = tmp146 * tmp137;
      const auto tmp150 = tmp149 / tmp148;
      double tmp151 = constant< 1 >();
      const auto tmp152 = tmp151 * tmp150;
      const auto tmp153 = tmp152 * tmp126;
      const auto tmp154 = std::pow( tmp153, 2 );
      const auto tmp155 = tmp154 / tmp108;
      const auto tmp156 = tmp155 * tmp112;
      const auto tmp157 = 0.5 * tmp156;
      const auto tmp158 = tmp157 + tmp110;
      return std::make_pair( RangeValueType{ { tmp158 } }, RangeValueType{ { tmp158 } } );
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
      const auto tmp72 = std::pow( tmp71, 2 );
      const auto tmp73 = tmp72 / tmp5;
      const auto tmp74 = std::get< 0 >( uOut );
      const auto tmp75 = -1 * tmp74[ 1 ];
      const auto tmp76 = std::get< 0 >( uIn );
      const auto tmp77 = tmp76[ 1 ] + tmp75;
      const auto tmp78 = 2 * tmp77;
      const auto tmp79 = tmp78 * tmp73;
      const auto tmp80 = 0.5 * tmp79;
      const auto tmp81 = tmp76[ 1 ] / tmp15;
      const auto tmp82 = std::pow( tmp81, 2 );
      const auto tmp83 = -1 + tmp21;
      const auto tmp84 = -1 * tmp76[ 1 ];
      const auto tmp85 = 1 + tmp84;
      const auto tmp86 = tmp85 + tmp14;
      const auto tmp87 = tmp86 / tmp15;
      const auto tmp88 = std::pow( tmp87, tmp83 );
      const auto tmp89 = -1 / tmp15;
      const auto tmp90 = tmp89 * tmp21;
      const auto tmp91 = tmp90 * tmp88;
      const auto tmp92 = -1 * tmp91;
      const auto tmp93 = tmp92 * tmp82;
      const auto tmp94 = 1 / tmp15;
      const auto tmp95 = 2 * tmp94;
      const auto tmp96 = tmp95 * tmp81;
      const auto tmp97 = std::pow( tmp87, tmp21 );
      const auto tmp98 = -1 * tmp97;
      const auto tmp99 = 1.0 + tmp98;
      const auto tmp100 = tmp99 * tmp96;
      const auto tmp101 = tmp100 + tmp93;
      const auto tmp102 = tmp101 / 0.0009;
      const auto tmp103 = tmp47 * tmp102;
      const auto tmp104 = std::get< 1 >( uIn );
      const auto tmp105 = 14322.6 + (tmp104[ 0 ])[ 1 ];
      const auto tmp106 = tmp105 * tmp103;
      const auto tmp107 = tmp43 * tmp102;
      const auto tmp108 = (tmp104[ 0 ])[ 0 ] * tmp107;
      const auto tmp109 = tmp108 + tmp106;
      const auto tmp110 = std::pow( tmp87, tmp31 );
      const auto tmp111 = tmp35 * tmp110;
      const auto tmp112 = tmp111 * tmp29;
      const auto tmp113 = tmp103 * tmp112;
      const auto tmp114 = -1 + tmp31;
      const auto tmp115 = std::pow( tmp87, tmp114 );
      const auto tmp116 = tmp31 * tmp89;
      const auto tmp117 = tmp116 * tmp115;
      const auto tmp118 = tmp35 * tmp117;
      const auto tmp119 = tmp118 * tmp29;
      const auto tmp120 = tmp99 * tmp82;
      const auto tmp121 = tmp120 / 0.0009;
      const auto tmp122 = tmp47 * tmp121;
      const auto tmp123 = tmp122 * tmp119;
      const auto tmp124 = tmp123 + tmp113;
      const auto tmp125 = (tmp104[ 1 ])[ 1 ] * tmp124;
      const auto tmp126 = tmp107 * tmp112;
      const auto tmp127 = tmp43 * tmp121;
      const auto tmp128 = tmp127 * tmp119;
      const auto tmp129 = tmp128 + tmp126;
      const auto tmp130 = (tmp104[ 1 ])[ 0 ] * tmp129;
      const auto tmp131 = tmp130 + tmp125;
      const auto tmp132 = tmp131 + tmp109;
      const auto tmp133 = tmp40[ 0 ] * tmp132;
      const auto tmp134 = (tmp104[ 0 ])[ 0 ] * tmp103;
      const auto tmp135 = tmp105 * tmp107;
      const auto tmp136 = tmp135 + tmp134;
      const auto tmp137 = (tmp104[ 1 ])[ 0 ] * tmp124;
      const auto tmp138 = (tmp104[ 1 ])[ 1 ] * tmp129;
      const auto tmp139 = tmp138 + tmp137;
      const auto tmp140 = tmp139 + tmp136;
      const auto tmp141 = tmp40[ 1 ] * tmp140;
      const auto tmp142 = tmp141 + tmp133;
      const auto tmp143 = 2 * tmp142;
      const auto tmp144 = tmp74[ 1 ] / tmp15;
      const auto tmp145 = std::pow( tmp144, 2 );
      const auto tmp146 = 1 + tmp75;
      const auto tmp147 = tmp146 + tmp14;
      const auto tmp148 = tmp147 / tmp15;
      const auto tmp149 = std::pow( tmp148, tmp21 );
      const auto tmp150 = -1 * tmp149;
      const auto tmp151 = 1.0 + tmp150;
      const auto tmp152 = tmp151 * tmp145;
      const auto tmp153 = tmp152 / 0.0009;
      const auto tmp154 = tmp47 * tmp153;
      const auto tmp155 = std::get< 1 >( uOut );
      const auto tmp156 = 14322.6 + (tmp155[ 0 ])[ 1 ];
      const auto tmp157 = tmp156 * tmp154;
      const auto tmp158 = tmp43 * tmp153;
      const auto tmp159 = (tmp155[ 0 ])[ 0 ] * tmp158;
      const auto tmp160 = tmp159 + tmp157;
      const auto tmp161 = std::pow( tmp148, tmp31 );
      const auto tmp162 = tmp35 * tmp161;
      const auto tmp163 = tmp162 * tmp29;
      const auto tmp164 = tmp154 * tmp163;
      const auto tmp165 = (tmp155[ 1 ])[ 1 ] * tmp164;
      const auto tmp166 = tmp158 * tmp163;
      const auto tmp167 = (tmp155[ 1 ])[ 0 ] * tmp166;
      const auto tmp168 = tmp167 + tmp165;
      const auto tmp169 = tmp168 + tmp160;
      const auto tmp170 = -1 * tmp169;
      const auto tmp171 = tmp105 * tmp122;
      const auto tmp172 = (tmp104[ 0 ])[ 0 ] * tmp127;
      const auto tmp173 = tmp172 + tmp171;
      const auto tmp174 = tmp122 * tmp112;
      const auto tmp175 = (tmp104[ 1 ])[ 1 ] * tmp174;
      const auto tmp176 = tmp127 * tmp112;
      const auto tmp177 = (tmp104[ 1 ])[ 0 ] * tmp176;
      const auto tmp178 = tmp177 + tmp175;
      const auto tmp179 = tmp178 + tmp173;
      const auto tmp180 = tmp179 + tmp170;
      const auto tmp181 = tmp40[ 0 ] * tmp180;
      const auto tmp182 = (tmp155[ 0 ])[ 0 ] * tmp154;
      const auto tmp183 = tmp156 * tmp158;
      const auto tmp184 = tmp183 + tmp182;
      const auto tmp185 = (tmp155[ 1 ])[ 0 ] * tmp164;
      const auto tmp186 = (tmp155[ 1 ])[ 1 ] * tmp166;
      const auto tmp187 = tmp186 + tmp185;
      const auto tmp188 = tmp187 + tmp184;
      const auto tmp189 = -1 * tmp188;
      const auto tmp190 = (tmp104[ 0 ])[ 0 ] * tmp122;
      const auto tmp191 = tmp105 * tmp127;
      const auto tmp192 = tmp191 + tmp190;
      const auto tmp193 = (tmp104[ 1 ])[ 0 ] * tmp174;
      const auto tmp194 = (tmp104[ 1 ])[ 1 ] * tmp176;
      const auto tmp195 = tmp194 + tmp193;
      const auto tmp196 = tmp195 + tmp192;
      const auto tmp197 = tmp196 + tmp189;
      const auto tmp198 = tmp40[ 1 ] * tmp197;
      const auto tmp199 = tmp198 + tmp181;
      const auto tmp200 = tmp199 * tmp143;
      const auto tmp201 = tmp200 * tmp5;
      const auto tmp202 = 0.5 * tmp201;
      const auto tmp203 = tmp202 + tmp80;
      const auto tmp204 = tmp40[ 1 ] * tmp122;
      const auto tmp205 = tmp40[ 0 ] * tmp127;
      const auto tmp206 = tmp205 + tmp204;
      const auto tmp207 = 2 * tmp206;
      const auto tmp208 = tmp199 * tmp207;
      const auto tmp209 = tmp208 * tmp5;
      const auto tmp210 = 0.5 * tmp209;
      const auto tmp211 = tmp40[ 0 ] * tmp122;
      const auto tmp212 = tmp40[ 1 ] * tmp127;
      const auto tmp213 = tmp212 + tmp211;
      const auto tmp214 = 2 * tmp213;
      const auto tmp215 = tmp199 * tmp214;
      const auto tmp216 = tmp215 * tmp5;
      const auto tmp217 = 0.5 * tmp216;
      const auto tmp218 = tmp40[ 1 ] * tmp174;
      const auto tmp219 = tmp40[ 0 ] * tmp176;
      const auto tmp220 = tmp219 + tmp218;
      const auto tmp221 = 2 * tmp220;
      const auto tmp222 = tmp199 * tmp221;
      const auto tmp223 = tmp222 * tmp5;
      const auto tmp224 = 0.5 * tmp223;
      const auto tmp225 = tmp40[ 0 ] * tmp174;
      const auto tmp226 = tmp40[ 1 ] * tmp176;
      const auto tmp227 = tmp226 + tmp225;
      const auto tmp228 = 2 * tmp227;
      const auto tmp229 = tmp199 * tmp228;
      const auto tmp230 = tmp229 * tmp5;
      const auto tmp231 = 0.5 * tmp230;
      const auto tmp232 = -2 * tmp77;
      const auto tmp233 = tmp232 * tmp73;
      const auto tmp234 = 0.5 * tmp233;
      const auto tmp235 = std::pow( tmp148, tmp83 );
      const auto tmp236 = tmp90 * tmp235;
      const auto tmp237 = -1 * tmp236;
      const auto tmp238 = tmp237 * tmp145;
      const auto tmp239 = tmp95 * tmp144;
      const auto tmp240 = tmp151 * tmp239;
      const auto tmp241 = tmp240 + tmp238;
      const auto tmp242 = tmp241 / 0.0009;
      const auto tmp243 = tmp47 * tmp242;
      const auto tmp244 = tmp156 * tmp243;
      const auto tmp245 = tmp43 * tmp242;
      const auto tmp246 = (tmp155[ 0 ])[ 0 ] * tmp245;
      const auto tmp247 = tmp246 + tmp244;
      const auto tmp248 = tmp243 * tmp163;
      const auto tmp249 = std::pow( tmp148, tmp114 );
      const auto tmp250 = tmp116 * tmp249;
      const auto tmp251 = tmp35 * tmp250;
      const auto tmp252 = tmp251 * tmp29;
      const auto tmp253 = tmp154 * tmp252;
      const auto tmp254 = tmp253 + tmp248;
      const auto tmp255 = (tmp155[ 1 ])[ 1 ] * tmp254;
      const auto tmp256 = tmp245 * tmp163;
      const auto tmp257 = tmp158 * tmp252;
      const auto tmp258 = tmp257 + tmp256;
      const auto tmp259 = (tmp155[ 1 ])[ 0 ] * tmp258;
      const auto tmp260 = tmp259 + tmp255;
      const auto tmp261 = tmp260 + tmp247;
      const auto tmp262 = -1 * tmp261;
      const auto tmp263 = tmp40[ 0 ] * tmp262;
      const auto tmp264 = (tmp155[ 0 ])[ 0 ] * tmp243;
      const auto tmp265 = tmp156 * tmp245;
      const auto tmp266 = tmp265 + tmp264;
      const auto tmp267 = (tmp155[ 1 ])[ 0 ] * tmp254;
      const auto tmp268 = (tmp155[ 1 ])[ 1 ] * tmp258;
      const auto tmp269 = tmp268 + tmp267;
      const auto tmp270 = tmp269 + tmp266;
      const auto tmp271 = -1 * tmp270;
      const auto tmp272 = tmp40[ 1 ] * tmp271;
      const auto tmp273 = tmp272 + tmp263;
      const auto tmp274 = 2 * tmp273;
      const auto tmp275 = tmp199 * tmp274;
      const auto tmp276 = tmp275 * tmp5;
      const auto tmp277 = 0.5 * tmp276;
      const auto tmp278 = tmp277 + tmp234;
      const auto tmp279 = -1 * tmp154;
      const auto tmp280 = tmp40[ 1 ] * tmp279;
      const auto tmp281 = -1 * tmp158;
      const auto tmp282 = tmp40[ 0 ] * tmp281;
      const auto tmp283 = tmp282 + tmp280;
      const auto tmp284 = 2 * tmp283;
      const auto tmp285 = tmp199 * tmp284;
      const auto tmp286 = tmp285 * tmp5;
      const auto tmp287 = 0.5 * tmp286;
      const auto tmp288 = tmp40[ 0 ] * tmp279;
      const auto tmp289 = tmp40[ 1 ] * tmp281;
      const auto tmp290 = tmp289 + tmp288;
      const auto tmp291 = 2 * tmp290;
      const auto tmp292 = tmp199 * tmp291;
      const auto tmp293 = tmp292 * tmp5;
      const auto tmp294 = 0.5 * tmp293;
      const auto tmp295 = -1 * tmp164;
      const auto tmp296 = tmp40[ 1 ] * tmp295;
      const auto tmp297 = -1 * tmp166;
      const auto tmp298 = tmp40[ 0 ] * tmp297;
      const auto tmp299 = tmp298 + tmp296;
      const auto tmp300 = 2 * tmp299;
      const auto tmp301 = tmp199 * tmp300;
      const auto tmp302 = tmp301 * tmp5;
      const auto tmp303 = 0.5 * tmp302;
      const auto tmp304 = tmp40[ 0 ] * tmp295;
      const auto tmp305 = tmp40[ 1 ] * tmp297;
      const auto tmp306 = tmp305 + tmp304;
      const auto tmp307 = 2 * tmp306;
      const auto tmp308 = tmp199 * tmp307;
      const auto tmp309 = tmp308 * tmp5;
      const auto tmp310 = 0.5 * tmp309;
      return std::make_pair( [ tmp203, tmp210, tmp217, tmp224, tmp231 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { (((tmp203 * (std::get< 0 >( phiIn ))[ 1 ] + tmp210 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp217 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ]) + tmp224 * ((std::get< 1 >( phiIn ))[ 1 ])[ 0 ]) + tmp231 * ((std::get< 1 >( phiIn ))[ 1 ])[ 1 ] } }, RangeValueType{ { (((tmp203 * (std::get< 0 >( phiIn ))[ 1 ] + tmp210 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp217 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ]) + tmp224 * ((std::get< 1 >( phiIn ))[ 1 ])[ 0 ]) + tmp231 * ((std::get< 1 >( phiIn ))[ 1 ])[ 1 ] } } );
        }, [ tmp278, tmp287, tmp294, tmp303, tmp310 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { (((tmp278 * (std::get< 0 >( phiOut ))[ 1 ] + tmp287 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp294 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ]) + tmp303 * ((std::get< 1 >( phiOut ))[ 1 ])[ 0 ]) + tmp310 * ((std::get< 1 >( phiOut ))[ 1 ])[ 1 ] } }, RangeValueType{ { (((tmp278 * (std::get< 0 >( phiOut ))[ 1 ] + tmp287 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp294 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ]) + tmp303 * ((std::get< 1 >( phiOut ))[ 1 ])[ 0 ]) + tmp310 * ((std::get< 1 >( phiOut ))[ 1 ])[ 1 ] } } );
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

    const ctype maxCellEdgeLength () const
    {
      return maxCellEdgeLength_[ (static_cast< std::size_t >( Side::in )) ];
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
    std::array< ctype, 2 > maxCellEdgeLength_;
    ctype facetArea_;
    ConstantTupleType constants_;
    std::array< std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution_old > >, 2 > coefficients_;
  };

} // namespace Integrands_e6d75985bf018b61c23a9baf5cbf37abv1_3

PYBIND11_MODULE( integrands_e6d75985bf018b61c23a9baf5cbf37abv1_3_220094643888a7089a01c0741f938529, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > > GridPart;
  typedef Integrands_e6d75985bf018b61c23a9baf5cbf37abv1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_e6d75985bf018b61c23a9baf5cbf37abv1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_e6d75985bf018b61c23a9baf5cbf37abv1_3_220094643888a7089a01c0741f938529.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffsolution_old ) { return new Integrands( coeffsolution_old ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property( "timeStep", [] ( Integrands &self ) -> Integrands::ContimeStep { return self.contimeStep(); }, [] ( Integrands &self, const Integrands::ContimeStep &v ) { self.contimeStep() = v; } );
      cls.def_property( "penalty", [] ( Integrands &self ) -> Integrands::Conpenalty { return self.conpenalty(); }, [] ( Integrands &self, const Integrands::Conpenalty &v ) { self.conpenalty() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_e6d75985bf018b61c23a9baf5cbf37abv1_3
