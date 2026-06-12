#ifndef GUARD_730d3cf3b7a88095397c27acd374c89d
#define GUARD_730d3cf3b7a88095397c27acd374c89d

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/common/exceptions.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_730d3cf3b7a88095397c27acd374c89d
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class CoeffgradSol, class Coeffsol >
struct UFLLocalFunction
  : public Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<3>>
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridViewType GridView;
  typedef typename GridView::ctype ctype;
  typedef Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<3>> BaseType;
  typedef Dune::Fem::GridFunctionSpace<GridPartType,Dune::Dim<3>> FunctionSpaceType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef typename EntityType::Geometry Geometry;
  typedef typename Geometry::GlobalCoordinate GlobalCoordinateType;
  typedef Dune::Fem::IntersectionSide Side;
  typedef std::tuple<> ConstantTupleType;
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 3 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< CoeffgradSol, Coeffsol > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<CoeffgradSol>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsol>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const CoeffgradSol &coeffgradSol, const Coeffsol &coeffsol, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< CoeffgradSol >( coeffgradSol ), Dune::Fem::ConstLocalFunction< Coeffsol >( coeffsol ) )
  {}

  void decouple ()
  {}

  void bind ( const IntersectionType &intersection, Side side )
  {
    BaseType::bind(intersection, side);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
    std::get< 1 >( coefficients_ ).bind( this->entity() );
  }

  void bind ( const EntityType &entity )
  {
    BaseType::bind(entity);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
    std::get< 1 >( coefficients_ ).bind( this->entity() );
  }

  void unbind ()
  {
    BaseType::unbind();
    std::get< 0 >( coefficients_ ).unbind();
    std::get< 1 >( coefficients_ ).unbind();
  }

  template< class Point >
  void evaluate ( const Point &x, typename FunctionSpaceType::RangeType &result ) const
  {
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 3.141592653589793 * tmp0[ 2 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::sin( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = 3.141592653589793 * tmp0[ 1 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = std::sin( tmp3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = 3.141592653589793 * tmp0[ 0 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = std::sin( tmp5 );
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 * tmp4;
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 * tmp2;
      assert(tmp8 == tmp8);
    const auto tmp9 = -1 * tmp8;
      assert(tmp9 == tmp9);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp10 = evaluateCoefficient< 1 >( x );
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10[ 0 ] + tmp9;
      assert(tmp11 == tmp11);
    const auto tmp12 = std::pow( tmp11, 2 );
      assert(tmp12 == tmp12);
    const auto tmp13 = std::cos( tmp1 );
      assert(tmp13 == tmp13);
    const auto tmp14 = 3.141592653589793 * tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 * tmp7;
      assert(tmp15 == tmp15);
    const auto tmp16 = -1 * tmp15;
      assert(tmp16 == tmp16);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp17 = evaluateCoefficient< 0 >( x );
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp17[ 2 ] + tmp16;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 * tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = std::cos( tmp5 );
      assert(tmp20 == tmp20);
    const auto tmp21 = 3.141592653589793 * tmp20;
      assert(tmp21 == tmp21);
    const auto tmp22 = tmp21 * tmp4;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 * tmp2;
      assert(tmp23 == tmp23);
    const auto tmp24 = -1 * tmp23;
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp17[ 0 ] + tmp24;
      assert(tmp25 == tmp25);
    const auto tmp26 = tmp25 * tmp25;
      assert(tmp26 == tmp26);
    const auto tmp27 = std::cos( tmp3 );
      assert(tmp27 == tmp27);
    const auto tmp28 = 3.141592653589793 * tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp28 * tmp6;
      assert(tmp29 == tmp29);
    const auto tmp30 = tmp29 * tmp2;
      assert(tmp30 == tmp30);
    const auto tmp31 = -1 * tmp30;
      assert(tmp31 == tmp31);
    const auto tmp32 = tmp17[ 1 ] + tmp31;
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp32 * tmp32;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 + tmp26;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp34 + tmp19;
      assert(tmp35 == tmp35);
    const auto tmp36 = -1 * tmp2;
      assert(tmp36 == tmp36);
    const auto tmp37 = 3.141592653589793 * tmp36;
      assert(tmp37 == tmp37);
    const auto tmp38 = 3.141592653589793 * tmp37;
      assert(tmp38 == tmp38);
    const auto tmp39 = tmp38 * tmp7;
      assert(tmp39 == tmp39);
    const auto tmp40 = -1 * tmp39;
      assert(tmp40 == tmp40);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp41 = jacobianCoefficient< 0 >( x );
      assert(tmp41 == tmp41);
    const auto tmp42 = (tmp41[ 2 ])[ 2 ] + tmp40;
      assert(tmp42 == tmp42);
    const auto tmp43 = -1 * tmp6;
      assert(tmp43 == tmp43);
    const auto tmp44 = 3.141592653589793 * tmp43;
      assert(tmp44 == tmp44);
    const auto tmp45 = 3.141592653589793 * tmp44;
      assert(tmp45 == tmp45);
    const auto tmp46 = tmp45 * tmp4;
      assert(tmp46 == tmp46);
    const auto tmp47 = tmp46 * tmp2;
      assert(tmp47 == tmp47);
    const auto tmp48 = -1 * tmp47;
      assert(tmp48 == tmp48);
    const auto tmp49 = (tmp41[ 0 ])[ 0 ] + tmp48;
      assert(tmp49 == tmp49);
    const auto tmp50 = -1 * tmp4;
      assert(tmp50 == tmp50);
    const auto tmp51 = 3.141592653589793 * tmp50;
      assert(tmp51 == tmp51);
    const auto tmp52 = 3.141592653589793 * tmp51;
      assert(tmp52 == tmp52);
    const auto tmp53 = tmp52 * tmp6;
      assert(tmp53 == tmp53);
    const auto tmp54 = tmp53 * tmp2;
      assert(tmp54 == tmp54);
    const auto tmp55 = -1 * tmp54;
      assert(tmp55 == tmp55);
    const auto tmp56 = (tmp41[ 1 ])[ 1 ] + tmp55;
      assert(tmp56 == tmp56);
    const auto tmp57 = tmp56 + tmp49;
      assert(tmp57 == tmp57);
    const auto tmp58 = tmp57 + tmp42;
      assert(tmp58 == tmp58);
    const auto tmp59 = std::pow( tmp58, 2 );
      assert(tmp59 == tmp59);
    result[ 0 ] = tmp12;
    result[ 1 ] = tmp35;
    result[ 2 ] = tmp59;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 3.141592653589793 * tmp0[ 2 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::sin( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = 3.141592653589793 * tmp0[ 1 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = std::sin( tmp3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = 3.141592653589793 * tmp0[ 0 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = std::cos( tmp5 );
      assert(tmp6 == tmp6);
    const auto tmp7 = 3.141592653589793 * tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 * tmp4;
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp8 * tmp2;
      assert(tmp9 == tmp9);
    const auto tmp10 = -1 * tmp9;
      assert(tmp10 == tmp10);
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp11 = jacobianCoefficient< 1 >( x );
      assert(tmp11 == tmp11);
    const auto tmp12 = (tmp11[ 0 ])[ 0 ] + tmp10;
      assert(tmp12 == tmp12);
    const auto tmp13 = 2 * tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = std::sin( tmp5 );
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 * tmp4;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 * tmp2;
      assert(tmp16 == tmp16);
    const auto tmp17 = -1 * tmp16;
      assert(tmp17 == tmp17);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp18 = evaluateCoefficient< 1 >( x );
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18[ 0 ] + tmp17;
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp19 * tmp13;
      assert(tmp20 == tmp20);
    const auto tmp21 = std::cos( tmp3 );
      assert(tmp21 == tmp21);
    const auto tmp22 = 3.141592653589793 * tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 * tmp14;
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp23 * tmp2;
      assert(tmp24 == tmp24);
    const auto tmp25 = -1 * tmp24;
      assert(tmp25 == tmp25);
    const auto tmp26 = (tmp11[ 0 ])[ 1 ] + tmp25;
      assert(tmp26 == tmp26);
    const auto tmp27 = 2 * tmp26;
      assert(tmp27 == tmp27);
    const auto tmp28 = tmp19 * tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = std::cos( tmp1 );
      assert(tmp29 == tmp29);
    const auto tmp30 = 3.141592653589793 * tmp29;
      assert(tmp30 == tmp30);
    const auto tmp31 = tmp30 * tmp15;
      assert(tmp31 == tmp31);
    const auto tmp32 = -1 * tmp31;
      assert(tmp32 == tmp32);
    const auto tmp33 = (tmp11[ 0 ])[ 2 ] + tmp32;
      assert(tmp33 == tmp33);
    const auto tmp34 = 2 * tmp33;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp19 * tmp34;
      assert(tmp35 == tmp35);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp36 = evaluateCoefficient< 0 >( x );
      assert(tmp36 == tmp36);
    const auto tmp37 = tmp36[ 2 ] + tmp32;
      assert(tmp37 == tmp37);
    const auto tmp38 = tmp30 * tmp8;
      assert(tmp38 == tmp38);
    const auto tmp39 = -1 * tmp38;
      assert(tmp39 == tmp39);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp40 = jacobianCoefficient< 0 >( x );
      assert(tmp40 == tmp40);
    const auto tmp41 = (tmp40[ 2 ])[ 0 ] + tmp39;
      assert(tmp41 == tmp41);
    const auto tmp42 = tmp41 * tmp37;
      assert(tmp42 == tmp42);
    const auto tmp43 = tmp37 * tmp41;
      assert(tmp43 == tmp43);
    const auto tmp44 = tmp43 + tmp42;
      assert(tmp44 == tmp44);
    const auto tmp45 = tmp36[ 0 ] + tmp10;
      assert(tmp45 == tmp45);
    const auto tmp46 = -1 * tmp14;
      assert(tmp46 == tmp46);
    const auto tmp47 = 3.141592653589793 * tmp46;
      assert(tmp47 == tmp47);
    const auto tmp48 = 3.141592653589793 * tmp47;
      assert(tmp48 == tmp48);
    const auto tmp49 = tmp48 * tmp4;
      assert(tmp49 == tmp49);
    const auto tmp50 = tmp49 * tmp2;
      assert(tmp50 == tmp50);
    const auto tmp51 = -1 * tmp50;
      assert(tmp51 == tmp51);
    const auto tmp52 = (tmp40[ 0 ])[ 0 ] + tmp51;
      assert(tmp52 == tmp52);
    const auto tmp53 = tmp52 * tmp45;
      assert(tmp53 == tmp53);
    const auto tmp54 = tmp45 * tmp52;
      assert(tmp54 == tmp54);
    const auto tmp55 = tmp54 + tmp53;
      assert(tmp55 == tmp55);
    const auto tmp56 = tmp36[ 1 ] + tmp25;
      assert(tmp56 == tmp56);
    const auto tmp57 = tmp7 * tmp22;
      assert(tmp57 == tmp57);
    const auto tmp58 = tmp57 * tmp2;
      assert(tmp58 == tmp58);
    const auto tmp59 = -1 * tmp58;
      assert(tmp59 == tmp59);
    const auto tmp60 = (tmp40[ 1 ])[ 0 ] + tmp59;
      assert(tmp60 == tmp60);
    const auto tmp61 = tmp60 * tmp56;
      assert(tmp61 == tmp61);
    const auto tmp62 = tmp56 * tmp60;
      assert(tmp62 == tmp62);
    const auto tmp63 = tmp62 + tmp61;
      assert(tmp63 == tmp63);
    const auto tmp64 = tmp63 + tmp55;
      assert(tmp64 == tmp64);
    const auto tmp65 = tmp64 + tmp44;
      assert(tmp65 == tmp65);
    const auto tmp66 = tmp30 * tmp23;
      assert(tmp66 == tmp66);
    const auto tmp67 = -1 * tmp66;
      assert(tmp67 == tmp67);
    const auto tmp68 = (tmp40[ 2 ])[ 1 ] + tmp67;
      assert(tmp68 == tmp68);
    const auto tmp69 = tmp68 * tmp37;
      assert(tmp69 == tmp69);
    const auto tmp70 = tmp37 * tmp68;
      assert(tmp70 == tmp70);
    const auto tmp71 = tmp70 + tmp69;
      assert(tmp71 == tmp71);
    const auto tmp72 = (tmp40[ 0 ])[ 1 ] + tmp59;
      assert(tmp72 == tmp72);
    const auto tmp73 = tmp72 * tmp45;
      assert(tmp73 == tmp73);
    const auto tmp74 = tmp45 * tmp72;
      assert(tmp74 == tmp74);
    const auto tmp75 = tmp74 + tmp73;
      assert(tmp75 == tmp75);
    const auto tmp76 = -1 * tmp4;
      assert(tmp76 == tmp76);
    const auto tmp77 = 3.141592653589793 * tmp76;
      assert(tmp77 == tmp77);
    const auto tmp78 = 3.141592653589793 * tmp77;
      assert(tmp78 == tmp78);
    const auto tmp79 = tmp78 * tmp14;
      assert(tmp79 == tmp79);
    const auto tmp80 = tmp79 * tmp2;
      assert(tmp80 == tmp80);
    const auto tmp81 = -1 * tmp80;
      assert(tmp81 == tmp81);
    const auto tmp82 = (tmp40[ 1 ])[ 1 ] + tmp81;
      assert(tmp82 == tmp82);
    const auto tmp83 = tmp82 * tmp56;
      assert(tmp83 == tmp83);
    const auto tmp84 = tmp56 * tmp82;
      assert(tmp84 == tmp84);
    const auto tmp85 = tmp84 + tmp83;
      assert(tmp85 == tmp85);
    const auto tmp86 = tmp85 + tmp75;
      assert(tmp86 == tmp86);
    const auto tmp87 = tmp86 + tmp71;
      assert(tmp87 == tmp87);
    const auto tmp88 = -1 * tmp2;
      assert(tmp88 == tmp88);
    const auto tmp89 = 3.141592653589793 * tmp88;
      assert(tmp89 == tmp89);
    const auto tmp90 = 3.141592653589793 * tmp89;
      assert(tmp90 == tmp90);
    const auto tmp91 = tmp90 * tmp15;
      assert(tmp91 == tmp91);
    const auto tmp92 = -1 * tmp91;
      assert(tmp92 == tmp92);
    const auto tmp93 = (tmp40[ 2 ])[ 2 ] + tmp92;
      assert(tmp93 == tmp93);
    const auto tmp94 = tmp93 * tmp37;
      assert(tmp94 == tmp94);
    const auto tmp95 = tmp37 * tmp93;
      assert(tmp95 == tmp95);
    const auto tmp96 = tmp95 + tmp94;
      assert(tmp96 == tmp96);
    const auto tmp97 = (tmp40[ 0 ])[ 2 ] + tmp39;
      assert(tmp97 == tmp97);
    const auto tmp98 = tmp97 * tmp45;
      assert(tmp98 == tmp98);
    const auto tmp99 = tmp45 * tmp97;
      assert(tmp99 == tmp99);
    const auto tmp100 = tmp99 + tmp98;
      assert(tmp100 == tmp100);
    const auto tmp101 = (tmp40[ 1 ])[ 2 ] + tmp67;
      assert(tmp101 == tmp101);
    const auto tmp102 = tmp101 * tmp56;
      assert(tmp102 == tmp102);
    const auto tmp103 = tmp56 * tmp101;
      assert(tmp103 == tmp103);
    const auto tmp104 = tmp103 + tmp102;
      assert(tmp104 == tmp104);
    const auto tmp105 = tmp104 + tmp100;
      assert(tmp105 == tmp105);
    const auto tmp106 = tmp105 + tmp96;
      assert(tmp106 == tmp106);
    const auto tmp107 = tmp90 * tmp8;
      assert(tmp107 == tmp107);
    const auto tmp108 = -1 * tmp107;
      assert(tmp108 == tmp108);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp109 = hessianCoefficient< 0 >( x );
      assert(tmp109 == tmp109);
    const auto tmp110 = ((tmp109[ 2 ])[ 2 ])[ 0 ] + tmp108;
      assert(tmp110 == tmp110);
    const auto tmp111 = -1 * tmp7;
      assert(tmp111 == tmp111);
    const auto tmp112 = 3.141592653589793 * tmp111;
      assert(tmp112 == tmp112);
    const auto tmp113 = 3.141592653589793 * tmp112;
      assert(tmp113 == tmp113);
    const auto tmp114 = tmp113 * tmp4;
      assert(tmp114 == tmp114);
    const auto tmp115 = tmp114 * tmp2;
      assert(tmp115 == tmp115);
    const auto tmp116 = -1 * tmp115;
      assert(tmp116 == tmp116);
    const auto tmp117 = ((tmp109[ 0 ])[ 0 ])[ 0 ] + tmp116;
      assert(tmp117 == tmp117);
    const auto tmp118 = tmp78 * tmp7;
      assert(tmp118 == tmp118);
    const auto tmp119 = tmp118 * tmp2;
      assert(tmp119 == tmp119);
    const auto tmp120 = -1 * tmp119;
      assert(tmp120 == tmp120);
    const auto tmp121 = ((tmp109[ 1 ])[ 1 ])[ 0 ] + tmp120;
      assert(tmp121 == tmp121);
    const auto tmp122 = tmp121 + tmp117;
      assert(tmp122 == tmp122);
    const auto tmp123 = tmp122 + tmp110;
      assert(tmp123 == tmp123);
    const auto tmp124 = 2 * tmp123;
      assert(tmp124 == tmp124);
    const auto tmp125 = tmp82 + tmp52;
      assert(tmp125 == tmp125);
    const auto tmp126 = tmp125 + tmp93;
      assert(tmp126 == tmp126);
    const auto tmp127 = tmp126 * tmp124;
      assert(tmp127 == tmp127);
    const auto tmp128 = tmp90 * tmp23;
      assert(tmp128 == tmp128);
    const auto tmp129 = -1 * tmp128;
      assert(tmp129 == tmp129);
    const auto tmp130 = ((tmp109[ 2 ])[ 2 ])[ 1 ] + tmp129;
      assert(tmp130 == tmp130);
    const auto tmp131 = -1 * tmp22;
      assert(tmp131 == tmp131);
    const auto tmp132 = 3.141592653589793 * tmp131;
      assert(tmp132 == tmp132);
    const auto tmp133 = 3.141592653589793 * tmp132;
      assert(tmp133 == tmp133);
    const auto tmp134 = tmp133 * tmp14;
      assert(tmp134 == tmp134);
    const auto tmp135 = tmp134 * tmp2;
      assert(tmp135 == tmp135);
    const auto tmp136 = -1 * tmp135;
      assert(tmp136 == tmp136);
    const auto tmp137 = ((tmp109[ 1 ])[ 1 ])[ 1 ] + tmp136;
      assert(tmp137 == tmp137);
    const auto tmp138 = tmp48 * tmp22;
      assert(tmp138 == tmp138);
    const auto tmp139 = tmp138 * tmp2;
      assert(tmp139 == tmp139);
    const auto tmp140 = -1 * tmp139;
      assert(tmp140 == tmp140);
    const auto tmp141 = ((tmp109[ 0 ])[ 0 ])[ 1 ] + tmp140;
      assert(tmp141 == tmp141);
    const auto tmp142 = tmp141 + tmp137;
      assert(tmp142 == tmp142);
    const auto tmp143 = tmp142 + tmp130;
      assert(tmp143 == tmp143);
    const auto tmp144 = 2 * tmp143;
      assert(tmp144 == tmp144);
    const auto tmp145 = tmp126 * tmp144;
      assert(tmp145 == tmp145);
    const auto tmp146 = -1 * tmp30;
      assert(tmp146 == tmp146);
    const auto tmp147 = 3.141592653589793 * tmp146;
      assert(tmp147 == tmp147);
    const auto tmp148 = 3.141592653589793 * tmp147;
      assert(tmp148 == tmp148);
    const auto tmp149 = tmp148 * tmp15;
      assert(tmp149 == tmp149);
    const auto tmp150 = -1 * tmp149;
      assert(tmp150 == tmp150);
    const auto tmp151 = ((tmp109[ 2 ])[ 2 ])[ 2 ] + tmp150;
      assert(tmp151 == tmp151);
    const auto tmp152 = tmp30 * tmp49;
      assert(tmp152 == tmp152);
    const auto tmp153 = -1 * tmp152;
      assert(tmp153 == tmp153);
    const auto tmp154 = ((tmp109[ 0 ])[ 0 ])[ 2 ] + tmp153;
      assert(tmp154 == tmp154);
    const auto tmp155 = tmp30 * tmp79;
      assert(tmp155 == tmp155);
    const auto tmp156 = -1 * tmp155;
      assert(tmp156 == tmp156);
    const auto tmp157 = ((tmp109[ 1 ])[ 1 ])[ 2 ] + tmp156;
      assert(tmp157 == tmp157);
    const auto tmp158 = tmp157 + tmp154;
      assert(tmp158 == tmp158);
    const auto tmp159 = tmp158 + tmp151;
      assert(tmp159 == tmp159);
    const auto tmp160 = 2 * tmp159;
      assert(tmp160 == tmp160);
    const auto tmp161 = tmp126 * tmp160;
      assert(tmp161 == tmp161);
    (result[ 0 ])[ 0 ] = tmp20;
    (result[ 0 ])[ 1 ] = tmp28;
    (result[ 0 ])[ 2 ] = tmp35;
    (result[ 1 ])[ 0 ] = tmp65;
    (result[ 1 ])[ 1 ] = tmp87;
    (result[ 1 ])[ 2 ] = tmp106;
    (result[ 2 ])[ 0 ] = tmp127;
    (result[ 2 ])[ 1 ] = tmp145;
    (result[ 2 ])[ 2 ] = tmp161;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    DUNE_THROW(Dune::NotImplemented,"hessian method could not be generated for local function (TooHighDerivative('CodeGenerator does not allow for third order derivatives, yet.'))");
    result=typename FunctionSpaceType::HessianRangeType(0);
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
  ConstantTupleType constants_;
  std::tuple< Dune::Fem::ConstLocalFunction< CoeffgradSol >, Dune::Fem::ConstLocalFunction< Coeffsol > > coefficients_;
};

} // namespace UFLLocalFunctions_730d3cf3b7a88095397c27acd374c89d

PYBIND11_MODULE( localfunction_730d3cf3b7a88095397c27acd374c89d_4eef12a61e9d9c71c0346a2a711eb17c, module )
{
  typedef UFLLocalFunctions_730d3cf3b7a88095397c27acd374c89d::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView>, Dune::FieldVector< double, 3 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_730d3cf3b7a88095397c27acd374c89d::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView>, Dune::FieldVector< double, 3 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_730d3cf3b7a88095397c27acd374c89d_4eef12a61e9d9c71c0346a2a711eb17c.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView>, Dune::FieldVector< double, 3 > > &coeffgradSol, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffsol ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView>(gridView),name,order, coeffgradSol, coeffsol); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
