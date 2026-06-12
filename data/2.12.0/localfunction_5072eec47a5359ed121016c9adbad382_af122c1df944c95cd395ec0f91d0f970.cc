#ifndef GUARD_5072eec47a5359ed121016c9adbad382
#define GUARD_5072eec47a5359ed121016c9adbad382

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/common/exceptions.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_5072eec47a5359ed121016c9adbad382
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffsolution >
struct UFLLocalFunction
  : public Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<1>>
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridViewType GridView;
  typedef typename GridView::ctype ctype;
  typedef Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<1>> BaseType;
  typedef Dune::Fem::GridFunctionSpace<GridPartType,Dune::Dim<1>> FunctionSpaceType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef typename EntityType::Geometry Geometry;
  typedef typename Geometry::GlobalCoordinate GlobalCoordinateType;
  typedef Dune::Fem::IntersectionSide Side;
  typedef std::tuple<> ConstantTupleType;
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< Coeffsolution > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsolution>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffsolution &coeffsolution, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffsolution >( coeffsolution ) )
  {}

  void decouple ()
  {}

  void bind ( const IntersectionType &intersection, Side side )
  {
    BaseType::bind(intersection, side);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
  }

  void bind ( const EntityType &entity )
  {
    BaseType::bind(entity);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
  }

  void unbind ()
  {
    BaseType::unbind();
    std::get< 0 >( coefficients_ ).unbind();
  }

  template< class Point >
  void evaluate ( const Point &x, typename FunctionSpaceType::RangeType &result ) const
  {
    using std::atan2;
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = std::atan2( tmp0[ 1 ], tmp0[ 0 ] );
      assert(tmp1 == tmp1);
    const auto tmp2 = (tmp0[ 1 ] < 0.0 ? 6.283185307179586 : 0.0) + tmp1;
      assert(tmp2 == tmp2);
    const auto tmp3 = 0.5625 * tmp2;
      assert(tmp3 == tmp3);
    const auto tmp4 = std::sin( tmp3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp0[ 1 ] * tmp0[ 1 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp0[ 0 ] * tmp0[ 0 ];
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 + tmp5;
      assert(tmp7 == tmp7);
    const auto tmp8 = std::pow( tmp7, -0.71875 );
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp0[ 1 ] + tmp0[ 1 ];
      assert(tmp9 == tmp9);
    const auto tmp10 = 0.28125 * tmp9;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * tmp8;
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 * tmp4;
      assert(tmp12 == tmp12);
    const auto tmp13 = std::pow( tmp7, 0.28125 );
      assert(tmp13 == tmp13);
    const auto tmp14 = std::cos( tmp3 );
      assert(tmp14 == tmp14);
    const auto tmp15 = std::pow( tmp0[ 1 ], 2 );
      assert(tmp15 == tmp15);
    const auto tmp16 = std::pow( tmp0[ 0 ], 2 );
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp16 + tmp15;
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp0[ 0 ] / tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = 0.5625 * tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp19 * tmp14;
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp20 * tmp13;
      assert(tmp21 == tmp21);
    const auto tmp22 = tmp21 + tmp12;
      assert(tmp22 == tmp22);
    const auto tmp23 = -1 * tmp22;
      assert(tmp23 == tmp23);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp24 = jacobianCoefficient< 0 >( x );
      assert(tmp24 == tmp24);
    const auto tmp25 = (tmp24[ 0 ])[ 1 ] + tmp23;
      assert(tmp25 == tmp25);
    const auto tmp26 = tmp25 * tmp25;
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp0[ 0 ] + tmp0[ 0 ];
      assert(tmp27 == tmp27);
    const auto tmp28 = 0.28125 * tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp28 * tmp8;
      assert(tmp29 == tmp29);
    const auto tmp30 = tmp29 * tmp4;
      assert(tmp30 == tmp30);
    const auto tmp31 = -1 * tmp0[ 1 ];
      assert(tmp31 == tmp31);
    const auto tmp32 = tmp31 / tmp17;
      assert(tmp32 == tmp32);
    const auto tmp33 = 0.5625 * tmp32;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 * tmp14;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp34 * tmp13;
      assert(tmp35 == tmp35);
    const auto tmp36 = tmp35 + tmp30;
      assert(tmp36 == tmp36);
    const auto tmp37 = -1 * tmp36;
      assert(tmp37 == tmp37);
    const auto tmp38 = (tmp24[ 0 ])[ 0 ] + tmp37;
      assert(tmp38 == tmp38);
    const auto tmp39 = tmp38 * tmp38;
      assert(tmp39 == tmp39);
    const auto tmp40 = tmp39 + tmp26;
      assert(tmp40 == tmp40);
    result[ 0 ] = tmp40;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::atan2;
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = std::atan2( tmp0[ 1 ], tmp0[ 0 ] );
      assert(tmp1 == tmp1);
    const auto tmp2 = (tmp0[ 1 ] < 0.0 ? 6.283185307179586 : 0.0) + tmp1;
      assert(tmp2 == tmp2);
    const auto tmp3 = 0.5625 * tmp2;
      assert(tmp3 == tmp3);
    const auto tmp4 = std::sin( tmp3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp0[ 1 ] * tmp0[ 1 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp0[ 0 ] * tmp0[ 0 ];
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 + tmp5;
      assert(tmp7 == tmp7);
    const auto tmp8 = std::pow( tmp7, -0.71875 );
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp0[ 1 ] + tmp0[ 1 ];
      assert(tmp9 == tmp9);
    const auto tmp10 = 0.28125 * tmp9;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * tmp8;
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 * tmp4;
      assert(tmp12 == tmp12);
    const auto tmp13 = std::pow( tmp7, 0.28125 );
      assert(tmp13 == tmp13);
    const auto tmp14 = std::cos( tmp3 );
      assert(tmp14 == tmp14);
    const auto tmp15 = std::pow( tmp0[ 1 ], 2 );
      assert(tmp15 == tmp15);
    const auto tmp16 = std::pow( tmp0[ 0 ], 2 );
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp16 + tmp15;
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp0[ 0 ] / tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = 0.5625 * tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp19 * tmp14;
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp20 * tmp13;
      assert(tmp21 == tmp21);
    const auto tmp22 = tmp21 + tmp12;
      assert(tmp22 == tmp22);
    const auto tmp23 = -1 * tmp22;
      assert(tmp23 == tmp23);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp24 = jacobianCoefficient< 0 >( x );
      assert(tmp24 == tmp24);
    const auto tmp25 = (tmp24[ 0 ])[ 1 ] + tmp23;
      assert(tmp25 == tmp25);
    const auto tmp26 = std::pow( tmp7, -1.71875 );
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp0[ 0 ] + tmp0[ 0 ];
      assert(tmp27 == tmp27);
    const auto tmp28 = -0.71875 * tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp28 * tmp26;
      assert(tmp29 == tmp29);
    const auto tmp30 = tmp10 * tmp29;
      assert(tmp30 == tmp30);
    const auto tmp31 = tmp30 * tmp4;
      assert(tmp31 == tmp31);
    const auto tmp32 = -1 * tmp0[ 1 ];
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp32 / tmp17;
      assert(tmp33 == tmp33);
    const auto tmp34 = 0.5625 * tmp33;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp34 * tmp14;
      assert(tmp35 == tmp35);
    const auto tmp36 = tmp11 * tmp35;
      assert(tmp36 == tmp36);
    const auto tmp37 = tmp36 + tmp31;
      assert(tmp37 == tmp37);
    const auto tmp38 = 2 * tmp0[ 0 ];
      assert(tmp38 == tmp38);
    const auto tmp39 = tmp38 * tmp18;
      assert(tmp39 == tmp39);
    const auto tmp40 = -1 * tmp39;
      assert(tmp40 == tmp40);
    const auto tmp41 = 1 + tmp40;
      assert(tmp41 == tmp41);
    const auto tmp42 = tmp41 / tmp17;
      assert(tmp42 == tmp42);
    const auto tmp43 = 0.5625 * tmp42;
      assert(tmp43 == tmp43);
    const auto tmp44 = tmp43 * tmp14;
      assert(tmp44 == tmp44);
    const auto tmp45 = -1 * tmp4;
      assert(tmp45 == tmp45);
    const auto tmp46 = tmp34 * tmp45;
      assert(tmp46 == tmp46);
    const auto tmp47 = tmp46 * tmp19;
      assert(tmp47 == tmp47);
    const auto tmp48 = tmp47 + tmp44;
      assert(tmp48 == tmp48);
    const auto tmp49 = tmp48 * tmp13;
      assert(tmp49 == tmp49);
    const auto tmp50 = 0.28125 * tmp27;
      assert(tmp50 == tmp50);
    const auto tmp51 = tmp50 * tmp8;
      assert(tmp51 == tmp51);
    const auto tmp52 = tmp51 * tmp20;
      assert(tmp52 == tmp52);
    const auto tmp53 = tmp52 + tmp49;
      assert(tmp53 == tmp53);
    const auto tmp54 = tmp53 + tmp37;
      assert(tmp54 == tmp54);
    const auto tmp55 = -1 * tmp54;
      assert(tmp55 == tmp55);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp56 = hessianCoefficient< 0 >( x );
      assert(tmp56 == tmp56);
    const auto tmp57 = ((tmp56[ 0 ])[ 1 ])[ 0 ] + tmp55;
      assert(tmp57 == tmp57);
    const auto tmp58 = tmp57 * tmp25;
      assert(tmp58 == tmp58);
    const auto tmp59 = tmp58 + tmp58;
      assert(tmp59 == tmp59);
    const auto tmp60 = tmp51 * tmp4;
      assert(tmp60 == tmp60);
    const auto tmp61 = tmp35 * tmp13;
      assert(tmp61 == tmp61);
    const auto tmp62 = tmp61 + tmp60;
      assert(tmp62 == tmp62);
    const auto tmp63 = -1 * tmp62;
      assert(tmp63 == tmp63);
    const auto tmp64 = (tmp24[ 0 ])[ 0 ] + tmp63;
      assert(tmp64 == tmp64);
    const auto tmp65 = 0.5625 * tmp8;
      assert(tmp65 == tmp65);
    const auto tmp66 = tmp50 * tmp29;
      assert(tmp66 == tmp66);
    const auto tmp67 = tmp66 + tmp65;
      assert(tmp67 == tmp67);
    const auto tmp68 = tmp67 * tmp4;
      assert(tmp68 == tmp68);
    const auto tmp69 = tmp51 * tmp35;
      assert(tmp69 == tmp69);
    const auto tmp70 = tmp69 + tmp68;
      assert(tmp70 == tmp70);
    const auto tmp71 = tmp38 * tmp33;
      assert(tmp71 == tmp71);
    const auto tmp72 = -1 * tmp71;
      assert(tmp72 == tmp72);
    const auto tmp73 = tmp72 / tmp17;
      assert(tmp73 == tmp73);
    const auto tmp74 = 0.5625 * tmp73;
      assert(tmp74 == tmp74);
    const auto tmp75 = tmp74 * tmp14;
      assert(tmp75 == tmp75);
    const auto tmp76 = tmp46 * tmp34;
      assert(tmp76 == tmp76);
    const auto tmp77 = tmp76 + tmp75;
      assert(tmp77 == tmp77);
    const auto tmp78 = tmp77 * tmp13;
      assert(tmp78 == tmp78);
    const auto tmp79 = tmp69 + tmp78;
      assert(tmp79 == tmp79);
    const auto tmp80 = tmp79 + tmp70;
      assert(tmp80 == tmp80);
    const auto tmp81 = -1 * tmp80;
      assert(tmp81 == tmp81);
    const auto tmp82 = ((tmp56[ 0 ])[ 0 ])[ 0 ] + tmp81;
      assert(tmp82 == tmp82);
    const auto tmp83 = tmp82 * tmp64;
      assert(tmp83 == tmp83);
    const auto tmp84 = tmp83 + tmp83;
      assert(tmp84 == tmp84);
    const auto tmp85 = tmp84 + tmp59;
      assert(tmp85 == tmp85);
    const auto tmp86 = -0.71875 * tmp9;
      assert(tmp86 == tmp86);
    const auto tmp87 = tmp86 * tmp26;
      assert(tmp87 == tmp87);
    const auto tmp88 = tmp10 * tmp87;
      assert(tmp88 == tmp88);
    const auto tmp89 = tmp88 + tmp65;
      assert(tmp89 == tmp89);
    const auto tmp90 = tmp89 * tmp4;
      assert(tmp90 == tmp90);
    const auto tmp91 = tmp11 * tmp20;
      assert(tmp91 == tmp91);
    const auto tmp92 = tmp91 + tmp90;
      assert(tmp92 == tmp92);
    const auto tmp93 = 2 * tmp0[ 1 ];
      assert(tmp93 == tmp93);
    const auto tmp94 = tmp93 * tmp18;
      assert(tmp94 == tmp94);
    const auto tmp95 = -1 * tmp94;
      assert(tmp95 == tmp95);
    const auto tmp96 = tmp95 / tmp17;
      assert(tmp96 == tmp96);
    const auto tmp97 = 0.5625 * tmp96;
      assert(tmp97 == tmp97);
    const auto tmp98 = tmp97 * tmp14;
      assert(tmp98 == tmp98);
    const auto tmp99 = tmp19 * tmp45;
      assert(tmp99 == tmp99);
    const auto tmp100 = tmp99 * tmp19;
      assert(tmp100 == tmp100);
    const auto tmp101 = tmp100 + tmp98;
      assert(tmp101 == tmp101);
    const auto tmp102 = tmp101 * tmp13;
      assert(tmp102 == tmp102);
    const auto tmp103 = tmp91 + tmp102;
      assert(tmp103 == tmp103);
    const auto tmp104 = tmp103 + tmp92;
      assert(tmp104 == tmp104);
    const auto tmp105 = -1 * tmp104;
      assert(tmp105 == tmp105);
    const auto tmp106 = ((tmp56[ 0 ])[ 1 ])[ 1 ] + tmp105;
      assert(tmp106 == tmp106);
    const auto tmp107 = tmp106 * tmp25;
      assert(tmp107 == tmp107);
    const auto tmp108 = tmp107 + tmp107;
      assert(tmp108 == tmp108);
    const auto tmp109 = tmp50 * tmp87;
      assert(tmp109 == tmp109);
    const auto tmp110 = tmp109 * tmp4;
      assert(tmp110 == tmp110);
    const auto tmp111 = tmp52 + tmp110;
      assert(tmp111 == tmp111);
    const auto tmp112 = tmp93 * tmp33;
      assert(tmp112 == tmp112);
    const auto tmp113 = -1 * tmp112;
      assert(tmp113 == tmp113);
    const auto tmp114 = -1 + tmp113;
      assert(tmp114 == tmp114);
    const auto tmp115 = tmp114 / tmp17;
      assert(tmp115 == tmp115);
    const auto tmp116 = 0.5625 * tmp115;
      assert(tmp116 == tmp116);
    const auto tmp117 = tmp116 * tmp14;
      assert(tmp117 == tmp117);
    const auto tmp118 = tmp99 * tmp34;
      assert(tmp118 == tmp118);
    const auto tmp119 = tmp118 + tmp117;
      assert(tmp119 == tmp119);
    const auto tmp120 = tmp119 * tmp13;
      assert(tmp120 == tmp120);
    const auto tmp121 = tmp36 + tmp120;
      assert(tmp121 == tmp121);
    const auto tmp122 = tmp121 + tmp111;
      assert(tmp122 == tmp122);
    const auto tmp123 = -1 * tmp122;
      assert(tmp123 == tmp123);
    const auto tmp124 = ((tmp56[ 0 ])[ 0 ])[ 1 ] + tmp123;
      assert(tmp124 == tmp124);
    const auto tmp125 = tmp124 * tmp64;
      assert(tmp125 == tmp125);
    const auto tmp126 = tmp125 + tmp125;
      assert(tmp126 == tmp126);
    const auto tmp127 = tmp126 + tmp108;
      assert(tmp127 == tmp127);
    (result[ 0 ])[ 0 ] = tmp85;
    (result[ 0 ])[ 1 ] = tmp127;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution > > coefficients_;
};

} // namespace UFLLocalFunctions_5072eec47a5359ed121016c9adbad382

PYBIND11_MODULE( localfunction_5072eec47a5359ed121016c9adbad382_af122c1df944c95cd395ec0f91d0f970, module )
{
  typedef UFLLocalFunctions_5072eec47a5359ed121016c9adbad382::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_5072eec47a5359ed121016c9adbad382::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_5072eec47a5359ed121016c9adbad382_af122c1df944c95cd395ec0f91d0f970.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > &coeffsolution ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>(gridView),name,order, coeffsolution); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
