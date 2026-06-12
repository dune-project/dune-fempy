#ifndef GuardIntegrands_2f94bce2e1cee34beba1852696343846v1_4
#define GuardIntegrands_2f94bce2e1cee34beba1852696343846v1_4
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
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_2f94bce2e1cee34beba1852696343846v1_4
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffprevious >
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
    typedef double Condt;
    typedef double Cont;
    typedef std::tuple< std::shared_ptr< Condt >, std::shared_ptr< Cont > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffprevious > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffprevious>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffprevious &coeffprevious, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< Coeffprevious >( coeffprevious ) )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Condt >( (Condt(0)) );
      std::get< 1 >( constants_ ) = std::make_shared< Cont >( (Cont(0)) );
    }

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< Condt >( (constant<0>()) );
      std::get< 1 >( constants_ ) = std::make_shared< Cont >( (constant<1>()) );
    }

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      geometry_.emplace( this->entity().geometry() );
      std::get< 0 >( coefficients_ ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_ ).unbind();
    }

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
      using std::exp;
      using std::pow;
      using std::sqrt;
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp1 = evaluateCoefficient< 0 >( x );
        assert(tmp1 == tmp1);
      const auto tmp2 = -1 * tmp1[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = std::get< 0 >( u );
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3[ 0 ] + tmp2;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 / tmp0;
        assert(tmp5 == tmp5);
      double tmp6 = constant< 1 >();
        assert(tmp6 == tmp6);
      const auto tmp7 = -2 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = std::exp( tmp7 );
        assert(tmp8 == tmp8);
      const auto tmp9 = -2 * tmp8;
        assert(tmp9 == tmp9);
      GlobalCoordinateType tmp10 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp10 == tmp10);
      const auto tmp11 = std::pow( tmp10[ 0 ], 3 );
        assert(tmp11 == tmp11);
      const auto tmp12 = std::pow( tmp10[ 1 ], 3 );
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp13 + tmp11;
        assert(tmp14 == tmp14);
      const auto tmp15 = 0.3333333333333333 * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = std::pow( tmp10[ 1 ], 2 );
        assert(tmp17 == tmp17);
      const auto tmp18 = std::pow( tmp10[ 0 ], 2 );
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp18 + tmp17;
        assert(tmp19 == tmp19);
      const auto tmp20 = 0.5 * tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 + tmp16;
        assert(tmp21 == tmp21);
      const auto tmp22 = 1 + tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = -1 + tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 * tmp9;
        assert(tmp24 == tmp24);
      const auto tmp25 = 3 * tmp18;
        assert(tmp25 == tmp25);
      const auto tmp26 = 0.3333333333333333 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = -1 * tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = 2 * tmp10[ 0 ];
        assert(tmp28 == tmp28);
      const auto tmp29 = 0.5 * tmp28;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 + tmp27;
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp30 * tmp8;
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 * tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = 3 * tmp17;
        assert(tmp33 == tmp33);
      const auto tmp34 = -1 * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = 0.3333333333333333 * tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = -1 * tmp35;
        assert(tmp36 == tmp36);
      const auto tmp37 = 2 * tmp10[ 1 ];
        assert(tmp37 == tmp37);
      const auto tmp38 = 0.5 * tmp37;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 + tmp36;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 * tmp8;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 * tmp40;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp41 + tmp32;
        assert(tmp42 == tmp42);
      const auto tmp43 = std::sqrt( tmp42 );
        assert(tmp43 == tmp43);
      const auto tmp44 = 4 * tmp43;
        assert(tmp44 == tmp44);
      const auto tmp45 = 1 + tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = std::sqrt( tmp45 );
        assert(tmp46 == tmp46);
      const auto tmp47 = 1 + tmp46;
        assert(tmp47 == tmp47);
      const auto tmp48 = 2 * tmp46;
        assert(tmp48 == tmp48);
      const auto tmp49 = 2 * tmp43;
        assert(tmp49 == tmp49);
      const auto tmp50 = 3 * tmp28;
        assert(tmp50 == tmp50);
      const auto tmp51 = 0.3333333333333333 * tmp50;
        assert(tmp51 == tmp51);
      const auto tmp52 = -1 * tmp51;
        assert(tmp52 == tmp52);
      const auto tmp53 = 1.0 + tmp52;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 * tmp8;
        assert(tmp54 == tmp54);
      const auto tmp55 = tmp54 * tmp31;
        assert(tmp55 == tmp55);
      const auto tmp56 = tmp55 + tmp55;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp56 / tmp49;
        assert(tmp57 == tmp57);
      const auto tmp58 = 4 * tmp57;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp58 / tmp48;
        assert(tmp59 == tmp59);
      const auto tmp60 = 2 / tmp47;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp60 * tmp59;
        assert(tmp61 == tmp61);
      const auto tmp62 = -1 * tmp61;
        assert(tmp62 == tmp62);
      const auto tmp63 = tmp62 / tmp47;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp31 * tmp63;
        assert(tmp64 == tmp64);
      const auto tmp65 = tmp54 * tmp60;
        assert(tmp65 == tmp65);
      const auto tmp66 = tmp65 + tmp64;
        assert(tmp66 == tmp66);
      const auto tmp67 = 3 * tmp37;
        assert(tmp67 == tmp67);
      const auto tmp68 = -1 * tmp67;
        assert(tmp68 == tmp68);
      const auto tmp69 = 0.3333333333333333 * tmp68;
        assert(tmp69 == tmp69);
      const auto tmp70 = -1 * tmp69;
        assert(tmp70 == tmp70);
      const auto tmp71 = 1.0 + tmp70;
        assert(tmp71 == tmp71);
      const auto tmp72 = tmp71 * tmp8;
        assert(tmp72 == tmp72);
      const auto tmp73 = tmp72 * tmp40;
        assert(tmp73 == tmp73);
      const auto tmp74 = tmp73 + tmp73;
        assert(tmp74 == tmp74);
      const auto tmp75 = tmp74 / tmp49;
        assert(tmp75 == tmp75);
      const auto tmp76 = 4 * tmp75;
        assert(tmp76 == tmp76);
      const auto tmp77 = tmp76 / tmp48;
        assert(tmp77 == tmp77);
      const auto tmp78 = tmp60 * tmp77;
        assert(tmp78 == tmp78);
      const auto tmp79 = -1 * tmp78;
        assert(tmp79 == tmp79);
      const auto tmp80 = tmp79 / tmp47;
        assert(tmp80 == tmp80);
      const auto tmp81 = tmp40 * tmp80;
        assert(tmp81 == tmp81);
      const auto tmp82 = tmp72 * tmp60;
        assert(tmp82 == tmp82);
      const auto tmp83 = tmp82 + tmp81;
        assert(tmp83 == tmp83);
      const auto tmp84 = tmp83 + tmp66;
        assert(tmp84 == tmp84);
      const auto tmp85 = -1 * tmp84;
        assert(tmp85 == tmp85);
      const auto tmp86 = tmp85 + tmp24;
        assert(tmp86 == tmp86);
      const auto tmp87 = tmp0 + tmp6;
        assert(tmp87 == tmp87);
      const auto tmp88 = -2 * tmp87;
        assert(tmp88 == tmp88);
      const auto tmp89 = std::exp( tmp88 );
        assert(tmp89 == tmp89);
      const auto tmp90 = -2 * tmp89;
        assert(tmp90 == tmp90);
      const auto tmp91 = tmp23 * tmp90;
        assert(tmp91 == tmp91);
      const auto tmp92 = tmp30 * tmp89;
        assert(tmp92 == tmp92);
      const auto tmp93 = tmp92 * tmp92;
        assert(tmp93 == tmp93);
      const auto tmp94 = tmp39 * tmp89;
        assert(tmp94 == tmp94);
      const auto tmp95 = tmp94 * tmp94;
        assert(tmp95 == tmp95);
      const auto tmp96 = tmp95 + tmp93;
        assert(tmp96 == tmp96);
      const auto tmp97 = std::sqrt( tmp96 );
        assert(tmp97 == tmp97);
      const auto tmp98 = 4 * tmp97;
        assert(tmp98 == tmp98);
      const auto tmp99 = 1 + tmp98;
        assert(tmp99 == tmp99);
      const auto tmp100 = std::sqrt( tmp99 );
        assert(tmp100 == tmp100);
      const auto tmp101 = 1 + tmp100;
        assert(tmp101 == tmp101);
      const auto tmp102 = 2 * tmp100;
        assert(tmp102 == tmp102);
      const auto tmp103 = 2 * tmp97;
        assert(tmp103 == tmp103);
      const auto tmp104 = tmp53 * tmp89;
        assert(tmp104 == tmp104);
      const auto tmp105 = tmp104 * tmp92;
        assert(tmp105 == tmp105);
      const auto tmp106 = tmp105 + tmp105;
        assert(tmp106 == tmp106);
      const auto tmp107 = tmp106 / tmp103;
        assert(tmp107 == tmp107);
      const auto tmp108 = 4 * tmp107;
        assert(tmp108 == tmp108);
      const auto tmp109 = tmp108 / tmp102;
        assert(tmp109 == tmp109);
      const auto tmp110 = 2 / tmp101;
        assert(tmp110 == tmp110);
      const auto tmp111 = tmp110 * tmp109;
        assert(tmp111 == tmp111);
      const auto tmp112 = -1 * tmp111;
        assert(tmp112 == tmp112);
      const auto tmp113 = tmp112 / tmp101;
        assert(tmp113 == tmp113);
      const auto tmp114 = tmp92 * tmp113;
        assert(tmp114 == tmp114);
      const auto tmp115 = tmp104 * tmp110;
        assert(tmp115 == tmp115);
      const auto tmp116 = tmp115 + tmp114;
        assert(tmp116 == tmp116);
      const auto tmp117 = tmp71 * tmp89;
        assert(tmp117 == tmp117);
      const auto tmp118 = tmp117 * tmp94;
        assert(tmp118 == tmp118);
      const auto tmp119 = tmp118 + tmp118;
        assert(tmp119 == tmp119);
      const auto tmp120 = tmp119 / tmp103;
        assert(tmp120 == tmp120);
      const auto tmp121 = 4 * tmp120;
        assert(tmp121 == tmp121);
      const auto tmp122 = tmp121 / tmp102;
        assert(tmp122 == tmp122);
      const auto tmp123 = tmp110 * tmp122;
        assert(tmp123 == tmp123);
      const auto tmp124 = -1 * tmp123;
        assert(tmp124 == tmp124);
      const auto tmp125 = tmp124 / tmp101;
        assert(tmp125 == tmp125);
      const auto tmp126 = tmp94 * tmp125;
        assert(tmp126 == tmp126);
      const auto tmp127 = tmp117 * tmp110;
        assert(tmp127 == tmp127);
      const auto tmp128 = tmp127 + tmp126;
        assert(tmp128 == tmp128);
      const auto tmp129 = tmp128 + tmp116;
        assert(tmp129 == tmp129);
      const auto tmp130 = -1 * tmp129;
        assert(tmp130 == tmp130);
      const auto tmp131 = tmp130 + tmp91;
        assert(tmp131 == tmp131);
      const auto tmp132 = tmp131 + tmp86;
        assert(tmp132 == tmp132);
      const auto tmp133 = 0.5 * tmp132;
        assert(tmp133 == tmp133);
      const auto tmp134 = -1 * tmp133;
        assert(tmp134 == tmp134);
      const auto tmp135 = tmp134 + tmp5;
        assert(tmp135 == tmp135);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp136 = jacobianCoefficient< 0 >( x );
        assert(tmp136 == tmp136);
      const auto tmp137 = (tmp136[ 0 ])[ 1 ] * (tmp136[ 0 ])[ 1 ];
        assert(tmp137 == tmp137);
      const auto tmp138 = (tmp136[ 0 ])[ 0 ] * (tmp136[ 0 ])[ 0 ];
        assert(tmp138 == tmp138);
      const auto tmp139 = tmp138 + tmp137;
        assert(tmp139 == tmp139);
      const auto tmp140 = std::sqrt( tmp139 );
        assert(tmp140 == tmp140);
      const auto tmp141 = 4 * tmp140;
        assert(tmp141 == tmp141);
      const auto tmp142 = 1 + tmp141;
        assert(tmp142 == tmp142);
      const auto tmp143 = std::sqrt( tmp142 );
        assert(tmp143 == tmp143);
      const auto tmp144 = 1 + tmp143;
        assert(tmp144 == tmp144);
      const auto tmp145 = 2 / tmp144;
        assert(tmp145 == tmp145);
      const auto tmp146 = (tmp136[ 0 ])[ 0 ] * tmp145;
        assert(tmp146 == tmp146);
      const auto tmp147 = 0.5 * tmp146;
        assert(tmp147 == tmp147);
      const auto tmp148 = std::get< 1 >( u );
        assert(tmp148 == tmp148);
      const auto tmp149 = (tmp148[ 0 ])[ 1 ] * (tmp148[ 0 ])[ 1 ];
        assert(tmp149 == tmp149);
      const auto tmp150 = (tmp148[ 0 ])[ 0 ] * (tmp148[ 0 ])[ 0 ];
        assert(tmp150 == tmp150);
      const auto tmp151 = tmp150 + tmp149;
        assert(tmp151 == tmp151);
      const auto tmp152 = std::sqrt( tmp151 );
        assert(tmp152 == tmp152);
      const auto tmp153 = 4 * tmp152;
        assert(tmp153 == tmp153);
      const auto tmp154 = 1 + tmp153;
        assert(tmp154 == tmp154);
      const auto tmp155 = std::sqrt( tmp154 );
        assert(tmp155 == tmp155);
      const auto tmp156 = 1 + tmp155;
        assert(tmp156 == tmp156);
      const auto tmp157 = 2 / tmp156;
        assert(tmp157 == tmp157);
      const auto tmp158 = (tmp148[ 0 ])[ 0 ] * tmp157;
        assert(tmp158 == tmp158);
      const auto tmp159 = 0.5 * tmp158;
        assert(tmp159 == tmp159);
      const auto tmp160 = tmp159 + tmp147;
        assert(tmp160 == tmp160);
      const auto tmp161 = (tmp136[ 0 ])[ 1 ] * tmp145;
        assert(tmp161 == tmp161);
      const auto tmp162 = 0.5 * tmp161;
        assert(tmp162 == tmp162);
      const auto tmp163 = (tmp148[ 0 ])[ 1 ] * tmp157;
        assert(tmp163 == tmp163);
      const auto tmp164 = 0.5 * tmp163;
        assert(tmp164 == tmp164);
      const auto tmp165 = tmp164 + tmp162;
        assert(tmp165 == tmp165);
      return RangeValueType{ { tmp135 }, { { tmp160, tmp165 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::sqrt;
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = 1 / tmp0;
        assert(tmp1 == tmp1);
      const auto tmp2 = std::get< 1 >( u );
        assert(tmp2 == tmp2);
      const auto tmp3 = (tmp2[ 0 ])[ 1 ] * (tmp2[ 0 ])[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = (tmp2[ 0 ])[ 0 ] * (tmp2[ 0 ])[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 + tmp3;
        assert(tmp5 == tmp5);
      const auto tmp6 = std::sqrt( tmp5 );
        assert(tmp6 == tmp6);
      const auto tmp7 = 4 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = 1 + tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = std::sqrt( tmp8 );
        assert(tmp9 == tmp9);
      const auto tmp10 = 1 + tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = 2 / tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = 2 * tmp9;
        assert(tmp12 == tmp12);
      const auto tmp13 = 2 * tmp6;
        assert(tmp13 == tmp13);
      const auto tmp14 = (tmp2[ 0 ])[ 0 ] + (tmp2[ 0 ])[ 0 ];
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 / tmp13;
        assert(tmp15 == tmp15);
      const auto tmp16 = 4 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 / tmp12;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp11 * tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = -1 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 / tmp10;
        assert(tmp20 == tmp20);
      const auto tmp21 = (tmp2[ 0 ])[ 0 ] * tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 + tmp11;
        assert(tmp22 == tmp22);
      const auto tmp23 = 0.5 * tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = (tmp2[ 0 ])[ 1 ] + (tmp2[ 0 ])[ 1 ];
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 / tmp13;
        assert(tmp25 == tmp25);
      const auto tmp26 = 4 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 / tmp12;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp11 * tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = -1 * tmp28;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 / tmp10;
        assert(tmp30 == tmp30);
      const auto tmp31 = (tmp2[ 0 ])[ 0 ] * tmp30;
        assert(tmp31 == tmp31);
      const auto tmp32 = 0.5 * tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = (tmp2[ 0 ])[ 1 ] * tmp20;
        assert(tmp33 == tmp33);
      const auto tmp34 = 0.5 * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = (tmp2[ 0 ])[ 1 ] * tmp30;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 + tmp11;
        assert(tmp36 == tmp36);
      const auto tmp37 = 0.5 * tmp36;
        assert(tmp37 == tmp37);
      return [ tmp1, tmp23, tmp32, tmp34, tmp37 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp1 * (std::get< 0 >( phi ))[ 0 ] }, { { tmp23 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp32 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], tmp34 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp37 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::exp;
      using std::pow;
      using std::sqrt;
      double tmp0 = constant< 1 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = -2 * tmp0;
        assert(tmp1 == tmp1);
      const auto tmp2 = std::exp( tmp1 );
        assert(tmp2 == tmp2);
      GlobalCoordinateType tmp3 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp3 == tmp3);
      const auto tmp4 = std::pow( tmp3[ 0 ], 2 );
        assert(tmp4 == tmp4);
      const auto tmp5 = 3 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = 0.3333333333333333 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = -1 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = 2 * tmp3[ 0 ];
        assert(tmp8 == tmp8);
      const auto tmp9 = 0.5 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 + tmp7;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * tmp2;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = std::pow( tmp3[ 1 ], 2 );
        assert(tmp13 == tmp13);
      const auto tmp14 = 3 * tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = -1 * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = 0.3333333333333333 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = -1 * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = 2 * tmp3[ 1 ];
        assert(tmp18 == tmp18);
      const auto tmp19 = 0.5 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 + tmp17;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 * tmp2;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 * tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 + tmp12;
        assert(tmp23 == tmp23);
      const auto tmp24 = std::sqrt( tmp23 );
        assert(tmp24 == tmp24);
      const auto tmp25 = 4 * tmp24;
        assert(tmp25 == tmp25);
      const auto tmp26 = 1 + tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = std::sqrt( tmp26 );
        assert(tmp27 == tmp27);
      const auto tmp28 = 1 + tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = 2 / tmp28;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp11 * tmp29;
        assert(tmp30 == tmp30);
      double tmp31 = constant< 0 >();
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 + tmp0;
        assert(tmp32 == tmp32);
      const auto tmp33 = -2 * tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = std::exp( tmp33 );
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp10 * tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 * tmp35;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp20 * tmp34;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp37 * tmp37;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 + tmp36;
        assert(tmp39 == tmp39);
      const auto tmp40 = std::sqrt( tmp39 );
        assert(tmp40 == tmp40);
      const auto tmp41 = 4 * tmp40;
        assert(tmp41 == tmp41);
      const auto tmp42 = 1 + tmp41;
        assert(tmp42 == tmp42);
      const auto tmp43 = std::sqrt( tmp42 );
        assert(tmp43 == tmp43);
      const auto tmp44 = 1 + tmp43;
        assert(tmp44 == tmp44);
      const auto tmp45 = 2 / tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp35 * tmp45;
        assert(tmp46 == tmp46);
      const auto tmp47 = tmp46 + tmp30;
        assert(tmp47 == tmp47);
      GlobalCoordinateType tmp48 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp48[ 0 ] * tmp47;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp21 * tmp29;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp37 * tmp45;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 + tmp50;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp48[ 1 ] * tmp52;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 + tmp49;
        assert(tmp54 == tmp54);
      const auto tmp55 = 0.5 * tmp54;
        assert(tmp55 == tmp55);
      const auto tmp56 = -1 * tmp55;
        assert(tmp56 == tmp56);
      return RangeValueType{ { tmp56 }, { { 0, 0 } } };
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

    const Condt &condt () const
    {
      return *std::get< 0 >( constants_ );
    }

    Condt &condt ()
    {
      return *std::get< 0 >( constants_ );
    }

    const Cont &cont () const
    {
      return *std::get< 1 >( constants_ );
    }

    Cont &cont ()
    {
      return *std::get< 1 >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_;
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::RangeType evaluateCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::RangeType result;
      std::get< i >( coefficients_ ).evaluate( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::JacobianRangeType result;
      std::get< i >( coefficients_ ).jacobian( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::HessianRangeType result;
      std::get< i >( coefficients_ ).hessian( x, result );
      return result;
    }
    EntityType entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    ConstantTupleType constants_;
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffprevious > > coefficients_;
  };

} // namespace Integrands_2f94bce2e1cee34beba1852696343846v1_4

PYBIND11_MODULE( integrands_2f94bce2e1cee34beba1852696343846v1_4_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_2f94bce2e1cee34beba1852696343846v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_2f94bce2e1cee34beba1852696343846v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_2f94bce2e1cee34beba1852696343846v1_4_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffprevious ) { return new Integrands( coeffprevious ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property( "dt", [] ( Integrands &self ) -> Integrands::Condt { return self.condt(); }, [] ( Integrands &self, const Integrands::Condt &v ) { self.condt() = v; } );
      cls.def_property( "t", [] ( Integrands &self ) -> Integrands::Cont { return self.cont(); }, [] ( Integrands &self, const Integrands::Cont &v ) { self.cont() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_2f94bce2e1cee34beba1852696343846v1_4
