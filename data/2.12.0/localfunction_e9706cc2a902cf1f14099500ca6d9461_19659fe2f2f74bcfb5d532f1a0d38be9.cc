#ifndef GUARD_e9706cc2a902cf1f14099500ca6d9461
#define GUARD_e9706cc2a902cf1f14099500ca6d9461

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

namespace UFLLocalFunctions_e9706cc2a902cf1f14099500ca6d9461
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffsolution_vector_dg, class Coeffanalytical_solution_x, class Coeffanalytical_solution_y, class Coeffanalytical_solution_pressure >
struct UFLLocalFunction
  : public Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<2>>
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridViewType GridView;
  typedef typename GridView::ctype ctype;
  typedef Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<2>> BaseType;
  typedef Dune::Fem::GridFunctionSpace<GridPartType,Dune::Dim<2>> FunctionSpaceType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef typename EntityType::Geometry Geometry;
  typedef typename Geometry::GlobalCoordinate GlobalCoordinateType;
  typedef Dune::Fem::IntersectionSide Side;
  typedef double Conavg_pressure;
  typedef std::tuple< std::shared_ptr< Conavg_pressure > > ConstantTupleType;
  template< std::size_t i >
  using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 3 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< Coeffsolution_vector_dg, Coeffanalytical_solution_x, Coeffanalytical_solution_y, Coeffanalytical_solution_pressure > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsolution_vector_dg>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_x>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_y>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_pressure>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffsolution_vector_dg &coeffsolution_vector_dg, const Coeffanalytical_solution_x &coeffanalytical_solution_x, const Coeffanalytical_solution_y &coeffanalytical_solution_y, const Coeffanalytical_solution_pressure &coeffanalytical_solution_pressure, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffsolution_vector_dg >( coeffsolution_vector_dg ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >( coeffanalytical_solution_x ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y >( coeffanalytical_solution_y ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_pressure >( coeffanalytical_solution_pressure ) )
  {
    std::get< 0 >( constants_ ) = std::make_shared< Conavg_pressure >( (Conavg_pressure(0)) );
  }

  void decouple ()
  {
    std::get< 0 >( constants_ ) = std::make_shared< Conavg_pressure >( (constant<0>()) );
  }

  void bind ( const IntersectionType &intersection, Side side )
  {
    BaseType::bind(intersection, side);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
    std::get< 1 >( coefficients_ ).bind( this->entity() );
    std::get< 2 >( coefficients_ ).bind( this->entity() );
    std::get< 3 >( coefficients_ ).bind( this->entity() );
  }

  void bind ( const EntityType &entity )
  {
    BaseType::bind(entity);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
    std::get< 1 >( coefficients_ ).bind( this->entity() );
    std::get< 2 >( coefficients_ ).bind( this->entity() );
    std::get< 3 >( coefficients_ ).bind( this->entity() );
  }

  void unbind ()
  {
    BaseType::unbind();
    std::get< 0 >( coefficients_ ).unbind();
    std::get< 1 >( coefficients_ ).unbind();
    std::get< 2 >( coefficients_ ).unbind();
    std::get< 3 >( coefficients_ ).unbind();
  }

  template< class Point >
  void evaluate ( const Point &x, typename FunctionSpaceType::RangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 3 >::RangeType tmp0 = evaluateCoefficient< 3 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 * tmp0[ 0 ];
      assert(tmp1 == tmp1);
    double tmp2 = constant< 0 >();
      assert(tmp2 == tmp2);
    const auto tmp3 = -1 * tmp2;
      assert(tmp3 == tmp3);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp4 = evaluateCoefficient< 0 >( x );
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp4[ 2 ] + tmp3;
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp5 + tmp1;
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 * tmp6;
      assert(tmp7 == tmp7);
    typename CoefficientFunctionSpaceType< 2 >::RangeType tmp8 = evaluateCoefficient< 2 >( x );
      assert(tmp8 == tmp8);
    const auto tmp9 = -1 * tmp8[ 0 ];
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp4[ 1 ] + tmp9;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * tmp10;
      assert(tmp11 == tmp11);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp12 = evaluateCoefficient< 1 >( x );
      assert(tmp12 == tmp12);
    const auto tmp13 = -1 * tmp12[ 0 ];
      assert(tmp13 == tmp13);
    const auto tmp14 = tmp4[ 0 ] + tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 * tmp14;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 + tmp11;
      assert(tmp16 == tmp16);
    result[ 0 ] = tmp7;
    result[ 1 ] = tmp16;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 3 >::JacobianRangeType tmp0 = jacobianCoefficient< 3 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 * (tmp0[ 0 ])[ 0 ];
      assert(tmp1 == tmp1);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp2 = jacobianCoefficient< 0 >( x );
      assert(tmp2 == tmp2);
    const auto tmp3 = (tmp2[ 2 ])[ 0 ] + tmp1;
      assert(tmp3 == tmp3);
    typename CoefficientFunctionSpaceType< 3 >::RangeType tmp4 = evaluateCoefficient< 3 >( x );
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * tmp4[ 0 ];
      assert(tmp5 == tmp5);
    double tmp6 = constant< 0 >();
      assert(tmp6 == tmp6);
    const auto tmp7 = -1 * tmp6;
      assert(tmp7 == tmp7);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp8 = evaluateCoefficient< 0 >( x );
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp8[ 2 ] + tmp7;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 + tmp5;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * tmp3;
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 + tmp11;
      assert(tmp12 == tmp12);
    const auto tmp13 = -1 * (tmp0[ 0 ])[ 1 ];
      assert(tmp13 == tmp13);
    const auto tmp14 = (tmp2[ 2 ])[ 1 ] + tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp10 * tmp14;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 + tmp15;
      assert(tmp16 == tmp16);
    typename CoefficientFunctionSpaceType< 2 >::JacobianRangeType tmp17 = jacobianCoefficient< 2 >( x );
      assert(tmp17 == tmp17);
    const auto tmp18 = -1 * (tmp17[ 0 ])[ 0 ];
      assert(tmp18 == tmp18);
    const auto tmp19 = (tmp2[ 1 ])[ 0 ] + tmp18;
      assert(tmp19 == tmp19);
    typename CoefficientFunctionSpaceType< 2 >::RangeType tmp20 = evaluateCoefficient< 2 >( x );
      assert(tmp20 == tmp20);
    const auto tmp21 = -1 * tmp20[ 0 ];
      assert(tmp21 == tmp21);
    const auto tmp22 = tmp8[ 1 ] + tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 * tmp19;
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp23 + tmp23;
      assert(tmp24 == tmp24);
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp25 = jacobianCoefficient< 1 >( x );
      assert(tmp25 == tmp25);
    const auto tmp26 = -1 * (tmp25[ 0 ])[ 0 ];
      assert(tmp26 == tmp26);
    const auto tmp27 = (tmp2[ 0 ])[ 0 ] + tmp26;
      assert(tmp27 == tmp27);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp28 = evaluateCoefficient< 1 >( x );
      assert(tmp28 == tmp28);
    const auto tmp29 = -1 * tmp28[ 0 ];
      assert(tmp29 == tmp29);
    const auto tmp30 = tmp8[ 0 ] + tmp29;
      assert(tmp30 == tmp30);
    const auto tmp31 = tmp30 * tmp27;
      assert(tmp31 == tmp31);
    const auto tmp32 = tmp31 + tmp31;
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp32 + tmp24;
      assert(tmp33 == tmp33);
    const auto tmp34 = -1 * (tmp17[ 0 ])[ 1 ];
      assert(tmp34 == tmp34);
    const auto tmp35 = (tmp2[ 1 ])[ 1 ] + tmp34;
      assert(tmp35 == tmp35);
    const auto tmp36 = tmp22 * tmp35;
      assert(tmp36 == tmp36);
    const auto tmp37 = tmp36 + tmp36;
      assert(tmp37 == tmp37);
    const auto tmp38 = -1 * (tmp25[ 0 ])[ 1 ];
      assert(tmp38 == tmp38);
    const auto tmp39 = (tmp2[ 0 ])[ 1 ] + tmp38;
      assert(tmp39 == tmp39);
    const auto tmp40 = tmp30 * tmp39;
      assert(tmp40 == tmp40);
    const auto tmp41 = tmp40 + tmp40;
      assert(tmp41 == tmp41);
    const auto tmp42 = tmp41 + tmp37;
      assert(tmp42 == tmp42);
    (result[ 0 ])[ 0 ] = tmp12;
    (result[ 0 ])[ 1 ] = tmp16;
    (result[ 1 ])[ 0 ] = tmp33;
    (result[ 1 ])[ 1 ] = tmp42;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 3 >::HessianRangeType tmp0 = hessianCoefficient< 3 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 * ((tmp0[ 0 ])[ 0 ])[ 0 ];
      assert(tmp1 == tmp1);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp2 = hessianCoefficient< 0 >( x );
      assert(tmp2 == tmp2);
    const auto tmp3 = ((tmp2[ 2 ])[ 0 ])[ 0 ] + tmp1;
      assert(tmp3 == tmp3);
    typename CoefficientFunctionSpaceType< 3 >::RangeType tmp4 = evaluateCoefficient< 3 >( x );
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * tmp4[ 0 ];
      assert(tmp5 == tmp5);
    double tmp6 = constant< 0 >();
      assert(tmp6 == tmp6);
    const auto tmp7 = -1 * tmp6;
      assert(tmp7 == tmp7);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp8 = evaluateCoefficient< 0 >( x );
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp8[ 2 ] + tmp7;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 + tmp5;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * tmp3;
      assert(tmp11 == tmp11);
    typename CoefficientFunctionSpaceType< 3 >::JacobianRangeType tmp12 = jacobianCoefficient< 3 >( x );
      assert(tmp12 == tmp12);
    const auto tmp13 = -1 * (tmp12[ 0 ])[ 0 ];
      assert(tmp13 == tmp13);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp14 = jacobianCoefficient< 0 >( x );
      assert(tmp14 == tmp14);
    const auto tmp15 = (tmp14[ 2 ])[ 0 ] + tmp13;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 * tmp15;
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp16 + tmp11;
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp17 + tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = -1 * (tmp12[ 0 ])[ 1 ];
      assert(tmp19 == tmp19);
    const auto tmp20 = (tmp14[ 2 ])[ 1 ] + tmp19;
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp15 * tmp20;
      assert(tmp21 == tmp21);
    const auto tmp22 = -1 * ((tmp0[ 0 ])[ 0 ])[ 1 ];
      assert(tmp22 == tmp22);
    const auto tmp23 = ((tmp2[ 2 ])[ 0 ])[ 1 ] + tmp22;
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp10 * tmp23;
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp24 + tmp21;
      assert(tmp25 == tmp25);
    const auto tmp26 = tmp25 + tmp25;
      assert(tmp26 == tmp26);
    const auto tmp27 = -1 * ((tmp0[ 0 ])[ 1 ])[ 0 ];
      assert(tmp27 == tmp27);
    const auto tmp28 = ((tmp2[ 2 ])[ 1 ])[ 0 ] + tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp10 * tmp28;
      assert(tmp29 == tmp29);
    const auto tmp30 = tmp21 + tmp29;
      assert(tmp30 == tmp30);
    const auto tmp31 = tmp30 + tmp30;
      assert(tmp31 == tmp31);
    const auto tmp32 = -1 * ((tmp0[ 0 ])[ 1 ])[ 1 ];
      assert(tmp32 == tmp32);
    const auto tmp33 = ((tmp2[ 2 ])[ 1 ])[ 1 ] + tmp32;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp10 * tmp33;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp20 * tmp20;
      assert(tmp35 == tmp35);
    const auto tmp36 = tmp35 + tmp34;
      assert(tmp36 == tmp36);
    const auto tmp37 = tmp36 + tmp36;
      assert(tmp37 == tmp37);
    typename CoefficientFunctionSpaceType< 2 >::HessianRangeType tmp38 = hessianCoefficient< 2 >( x );
      assert(tmp38 == tmp38);
    const auto tmp39 = -1 * ((tmp38[ 0 ])[ 0 ])[ 0 ];
      assert(tmp39 == tmp39);
    const auto tmp40 = ((tmp2[ 1 ])[ 0 ])[ 0 ] + tmp39;
      assert(tmp40 == tmp40);
    typename CoefficientFunctionSpaceType< 2 >::RangeType tmp41 = evaluateCoefficient< 2 >( x );
      assert(tmp41 == tmp41);
    const auto tmp42 = -1 * tmp41[ 0 ];
      assert(tmp42 == tmp42);
    const auto tmp43 = tmp8[ 1 ] + tmp42;
      assert(tmp43 == tmp43);
    const auto tmp44 = tmp43 * tmp40;
      assert(tmp44 == tmp44);
    typename CoefficientFunctionSpaceType< 2 >::JacobianRangeType tmp45 = jacobianCoefficient< 2 >( x );
      assert(tmp45 == tmp45);
    const auto tmp46 = -1 * (tmp45[ 0 ])[ 0 ];
      assert(tmp46 == tmp46);
    const auto tmp47 = (tmp14[ 1 ])[ 0 ] + tmp46;
      assert(tmp47 == tmp47);
    const auto tmp48 = tmp47 * tmp47;
      assert(tmp48 == tmp48);
    const auto tmp49 = tmp48 + tmp44;
      assert(tmp49 == tmp49);
    const auto tmp50 = tmp49 + tmp49;
      assert(tmp50 == tmp50);
    typename CoefficientFunctionSpaceType< 1 >::HessianRangeType tmp51 = hessianCoefficient< 1 >( x );
      assert(tmp51 == tmp51);
    const auto tmp52 = -1 * ((tmp51[ 0 ])[ 0 ])[ 0 ];
      assert(tmp52 == tmp52);
    const auto tmp53 = ((tmp2[ 0 ])[ 0 ])[ 0 ] + tmp52;
      assert(tmp53 == tmp53);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp54 = evaluateCoefficient< 1 >( x );
      assert(tmp54 == tmp54);
    const auto tmp55 = -1 * tmp54[ 0 ];
      assert(tmp55 == tmp55);
    const auto tmp56 = tmp8[ 0 ] + tmp55;
      assert(tmp56 == tmp56);
    const auto tmp57 = tmp56 * tmp53;
      assert(tmp57 == tmp57);
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp58 = jacobianCoefficient< 1 >( x );
      assert(tmp58 == tmp58);
    const auto tmp59 = -1 * (tmp58[ 0 ])[ 0 ];
      assert(tmp59 == tmp59);
    const auto tmp60 = (tmp14[ 0 ])[ 0 ] + tmp59;
      assert(tmp60 == tmp60);
    const auto tmp61 = tmp60 * tmp60;
      assert(tmp61 == tmp61);
    const auto tmp62 = tmp61 + tmp57;
      assert(tmp62 == tmp62);
    const auto tmp63 = tmp62 + tmp62;
      assert(tmp63 == tmp63);
    const auto tmp64 = tmp63 + tmp50;
      assert(tmp64 == tmp64);
    const auto tmp65 = -1 * (tmp45[ 0 ])[ 1 ];
      assert(tmp65 == tmp65);
    const auto tmp66 = (tmp14[ 1 ])[ 1 ] + tmp65;
      assert(tmp66 == tmp66);
    const auto tmp67 = tmp47 * tmp66;
      assert(tmp67 == tmp67);
    const auto tmp68 = -1 * ((tmp38[ 0 ])[ 0 ])[ 1 ];
      assert(tmp68 == tmp68);
    const auto tmp69 = ((tmp2[ 1 ])[ 0 ])[ 1 ] + tmp68;
      assert(tmp69 == tmp69);
    const auto tmp70 = tmp43 * tmp69;
      assert(tmp70 == tmp70);
    const auto tmp71 = tmp70 + tmp67;
      assert(tmp71 == tmp71);
    const auto tmp72 = tmp71 + tmp71;
      assert(tmp72 == tmp72);
    const auto tmp73 = -1 * (tmp58[ 0 ])[ 1 ];
      assert(tmp73 == tmp73);
    const auto tmp74 = (tmp14[ 0 ])[ 1 ] + tmp73;
      assert(tmp74 == tmp74);
    const auto tmp75 = tmp60 * tmp74;
      assert(tmp75 == tmp75);
    const auto tmp76 = -1 * ((tmp51[ 0 ])[ 0 ])[ 1 ];
      assert(tmp76 == tmp76);
    const auto tmp77 = ((tmp2[ 0 ])[ 0 ])[ 1 ] + tmp76;
      assert(tmp77 == tmp77);
    const auto tmp78 = tmp56 * tmp77;
      assert(tmp78 == tmp78);
    const auto tmp79 = tmp78 + tmp75;
      assert(tmp79 == tmp79);
    const auto tmp80 = tmp79 + tmp79;
      assert(tmp80 == tmp80);
    const auto tmp81 = tmp80 + tmp72;
      assert(tmp81 == tmp81);
    const auto tmp82 = -1 * ((tmp38[ 0 ])[ 1 ])[ 0 ];
      assert(tmp82 == tmp82);
    const auto tmp83 = ((tmp2[ 1 ])[ 1 ])[ 0 ] + tmp82;
      assert(tmp83 == tmp83);
    const auto tmp84 = tmp43 * tmp83;
      assert(tmp84 == tmp84);
    const auto tmp85 = tmp67 + tmp84;
      assert(tmp85 == tmp85);
    const auto tmp86 = tmp85 + tmp85;
      assert(tmp86 == tmp86);
    const auto tmp87 = -1 * ((tmp51[ 0 ])[ 1 ])[ 0 ];
      assert(tmp87 == tmp87);
    const auto tmp88 = ((tmp2[ 0 ])[ 1 ])[ 0 ] + tmp87;
      assert(tmp88 == tmp88);
    const auto tmp89 = tmp56 * tmp88;
      assert(tmp89 == tmp89);
    const auto tmp90 = tmp75 + tmp89;
      assert(tmp90 == tmp90);
    const auto tmp91 = tmp90 + tmp90;
      assert(tmp91 == tmp91);
    const auto tmp92 = tmp91 + tmp86;
      assert(tmp92 == tmp92);
    const auto tmp93 = -1 * ((tmp38[ 0 ])[ 1 ])[ 1 ];
      assert(tmp93 == tmp93);
    const auto tmp94 = ((tmp2[ 1 ])[ 1 ])[ 1 ] + tmp93;
      assert(tmp94 == tmp94);
    const auto tmp95 = tmp43 * tmp94;
      assert(tmp95 == tmp95);
    const auto tmp96 = tmp66 * tmp66;
      assert(tmp96 == tmp96);
    const auto tmp97 = tmp96 + tmp95;
      assert(tmp97 == tmp97);
    const auto tmp98 = tmp97 + tmp97;
      assert(tmp98 == tmp98);
    const auto tmp99 = -1 * ((tmp51[ 0 ])[ 1 ])[ 1 ];
      assert(tmp99 == tmp99);
    const auto tmp100 = ((tmp2[ 0 ])[ 1 ])[ 1 ] + tmp99;
      assert(tmp100 == tmp100);
    const auto tmp101 = tmp56 * tmp100;
      assert(tmp101 == tmp101);
    const auto tmp102 = tmp74 * tmp74;
      assert(tmp102 == tmp102);
    const auto tmp103 = tmp102 + tmp101;
      assert(tmp103 == tmp103);
    const auto tmp104 = tmp103 + tmp103;
      assert(tmp104 == tmp104);
    const auto tmp105 = tmp104 + tmp98;
      assert(tmp105 == tmp105);
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp18;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp26;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp31;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp37;
    ((result[ 1 ])[ 0 ])[ 0 ] = tmp64;
    ((result[ 1 ])[ 0 ])[ 1 ] = tmp81;
    ((result[ 1 ])[ 1 ])[ 0 ] = tmp92;
    ((result[ 1 ])[ 1 ])[ 1 ] = tmp105;
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

  const Conavg_pressure &conavg_pressure () const
  {
    return *std::get< 0 >( constants_ );
  }

  Conavg_pressure &conavg_pressure ()
  {
    return *std::get< 0 >( constants_ );
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution_vector_dg >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_pressure > > coefficients_;
};

} // namespace UFLLocalFunctions_e9706cc2a902cf1f14099500ca6d9461

PYBIND11_MODULE( localfunction_e9706cc2a902cf1f14099500ca6d9461_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_e9706cc2a902cf1f14099500ca6d9461::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 3 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_e9706cc2a902cf1f14099500ca6d9461::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 3 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_e9706cc2a902cf1f14099500ca6d9461_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 3 > > &coeffsolution_vector_dg, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_x, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_y, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_pressure ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffsolution_vector_dg, coeffanalytical_solution_x, coeffanalytical_solution_y, coeffanalytical_solution_pressure); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >(), pybind11::keep_alive< 1, 5 >(), pybind11::keep_alive< 1, 6 >() );
      cls.def_property( "avg_pressure", [] ( LocalFunctionType &self ) -> UFLLocalFunctions_e9706cc2a902cf1f14099500ca6d9461::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 3 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >::Conavg_pressure { return self.conavg_pressure(); }, [] ( LocalFunctionType &self, const UFLLocalFunctions_e9706cc2a902cf1f14099500ca6d9461::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 3 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >::Conavg_pressure &v ) { self.conavg_pressure() = v; } );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
