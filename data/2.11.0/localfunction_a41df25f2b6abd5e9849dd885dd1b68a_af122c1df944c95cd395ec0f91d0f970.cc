#ifndef GUARD_a41df25f2b6abd5e9849dd885dd1b68a
#define GUARD_a41df25f2b6abd5e9849dd885dd1b68a

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
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

namespace UFLLocalFunctions_a41df25f2b6abd5e9849dd885dd1b68a
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
    using std::abs;
    using std::atan2;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = std::atan2( tmp0[ 1 ], tmp0[ 0 ] );
    const auto tmp2 = (tmp0[ 1 ] < 0.0 ? 6.283185307179586 : 0.0) + tmp1;
    const auto tmp3 = 0.5625 * tmp2;
    const auto tmp4 = std::sin( tmp3 );
    const auto tmp5 = tmp0[ 1 ] * tmp0[ 1 ];
    const auto tmp6 = tmp0[ 0 ] * tmp0[ 0 ];
    const auto tmp7 = tmp6 + tmp5;
    const auto tmp8 = std::pow( tmp7, 0.28125 );
    const auto tmp9 = tmp8 * tmp4;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp10 = evaluateCoefficient< 0 >( x );
    const auto tmp11 = -1 * tmp10[ 0 ];
    const auto tmp12 = tmp11 + tmp9;
    const auto tmp13 = std::abs( tmp12 );
    result[ 0 ] = tmp13;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::atan2;
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = std::atan2( tmp0[ 1 ], tmp0[ 0 ] );
    const auto tmp2 = (tmp0[ 1 ] < 0.0 ? 6.283185307179586 : 0.0) + tmp1;
    const auto tmp3 = 0.5625 * tmp2;
    const auto tmp4 = std::sin( tmp3 );
    const auto tmp5 = tmp0[ 1 ] * tmp0[ 1 ];
    const auto tmp6 = tmp0[ 0 ] * tmp0[ 0 ];
    const auto tmp7 = tmp6 + tmp5;
    const auto tmp8 = std::pow( tmp7, 0.28125 );
    const auto tmp9 = tmp8 * tmp4;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp10 = evaluateCoefficient< 0 >( x );
    const auto tmp11 = -1 * tmp10[ 0 ];
    const auto tmp12 = tmp11 + tmp9;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp13 = jacobianCoefficient< 0 >( x );
    const auto tmp14 = -1 * (tmp13[ 0 ])[ 0 ];
    const auto tmp15 = std::pow( tmp7, -0.71875 );
    const auto tmp16 = tmp0[ 0 ] + tmp0[ 0 ];
    const auto tmp17 = 0.28125 * tmp16;
    const auto tmp18 = tmp17 * tmp15;
    const auto tmp19 = tmp18 * tmp4;
    const auto tmp20 = std::cos( tmp3 );
    const auto tmp21 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp22 = std::pow( tmp0[ 0 ], 2 );
    const auto tmp23 = tmp22 + tmp21;
    const auto tmp24 = -1 * tmp0[ 1 ];
    const auto tmp25 = tmp24 / tmp23;
    const auto tmp26 = 0.5625 * tmp25;
    const auto tmp27 = tmp26 * tmp20;
    const auto tmp28 = tmp27 * tmp8;
    const auto tmp29 = tmp28 + tmp19;
    const auto tmp30 = tmp29 + tmp14;
    const auto tmp31 = tmp30 * (tmp12 == 0.0 ? 0.0 : (tmp12 < 0.0 ? -1 : 1));
    const auto tmp32 = -1 * (tmp13[ 0 ])[ 1 ];
    const auto tmp33 = tmp0[ 1 ] + tmp0[ 1 ];
    const auto tmp34 = 0.28125 * tmp33;
    const auto tmp35 = tmp34 * tmp15;
    const auto tmp36 = tmp35 * tmp4;
    const auto tmp37 = tmp0[ 0 ] / tmp23;
    const auto tmp38 = 0.5625 * tmp37;
    const auto tmp39 = tmp38 * tmp20;
    const auto tmp40 = tmp39 * tmp8;
    const auto tmp41 = tmp40 + tmp36;
    const auto tmp42 = tmp41 + tmp32;
    const auto tmp43 = tmp42 * (tmp12 == 0.0 ? 0.0 : (tmp12 < 0.0 ? -1 : 1));
    (result[ 0 ])[ 0 ] = tmp31;
    (result[ 0 ])[ 1 ] = tmp43;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::atan2;
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = std::atan2( tmp0[ 1 ], tmp0[ 0 ] );
    const auto tmp2 = (tmp0[ 1 ] < 0.0 ? 6.283185307179586 : 0.0) + tmp1;
    const auto tmp3 = 0.5625 * tmp2;
    const auto tmp4 = std::sin( tmp3 );
    const auto tmp5 = tmp0[ 1 ] * tmp0[ 1 ];
    const auto tmp6 = tmp0[ 0 ] * tmp0[ 0 ];
    const auto tmp7 = tmp6 + tmp5;
    const auto tmp8 = std::pow( tmp7, 0.28125 );
    const auto tmp9 = tmp8 * tmp4;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp10 = evaluateCoefficient< 0 >( x );
    const auto tmp11 = -1 * tmp10[ 0 ];
    const auto tmp12 = tmp11 + tmp9;
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp13 = hessianCoefficient< 0 >( x );
    const auto tmp14 = -1 * ((tmp13[ 0 ])[ 0 ])[ 0 ];
    const auto tmp15 = std::pow( tmp7, -0.71875 );
    const auto tmp16 = 0.5625 * tmp15;
    const auto tmp17 = std::pow( tmp7, -1.71875 );
    const auto tmp18 = tmp0[ 0 ] + tmp0[ 0 ];
    const auto tmp19 = -0.71875 * tmp18;
    const auto tmp20 = tmp19 * tmp17;
    const auto tmp21 = 0.28125 * tmp18;
    const auto tmp22 = tmp21 * tmp20;
    const auto tmp23 = tmp22 + tmp16;
    const auto tmp24 = tmp23 * tmp4;
    const auto tmp25 = std::cos( tmp3 );
    const auto tmp26 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp27 = std::pow( tmp0[ 0 ], 2 );
    const auto tmp28 = tmp27 + tmp26;
    const auto tmp29 = -1 * tmp0[ 1 ];
    const auto tmp30 = tmp29 / tmp28;
    const auto tmp31 = 0.5625 * tmp30;
    const auto tmp32 = tmp31 * tmp25;
    const auto tmp33 = tmp21 * tmp15;
    const auto tmp34 = tmp33 * tmp32;
    const auto tmp35 = tmp34 + tmp24;
    const auto tmp36 = 2 * tmp0[ 0 ];
    const auto tmp37 = tmp36 * tmp30;
    const auto tmp38 = -1 * tmp37;
    const auto tmp39 = tmp38 / tmp28;
    const auto tmp40 = 0.5625 * tmp39;
    const auto tmp41 = tmp40 * tmp25;
    const auto tmp42 = -1 * tmp4;
    const auto tmp43 = tmp31 * tmp42;
    const auto tmp44 = tmp43 * tmp31;
    const auto tmp45 = tmp44 + tmp41;
    const auto tmp46 = tmp45 * tmp8;
    const auto tmp47 = tmp34 + tmp46;
    const auto tmp48 = tmp47 + tmp35;
    const auto tmp49 = tmp48 + tmp14;
    const auto tmp50 = tmp49 * (tmp12 == 0.0 ? 0.0 : (tmp12 < 0.0 ? -1 : 1));
    const auto tmp51 = -1 * ((tmp13[ 0 ])[ 0 ])[ 1 ];
    const auto tmp52 = tmp0[ 1 ] + tmp0[ 1 ];
    const auto tmp53 = -0.71875 * tmp52;
    const auto tmp54 = tmp53 * tmp17;
    const auto tmp55 = tmp21 * tmp54;
    const auto tmp56 = tmp55 * tmp4;
    const auto tmp57 = tmp0[ 0 ] / tmp28;
    const auto tmp58 = 0.5625 * tmp57;
    const auto tmp59 = tmp58 * tmp25;
    const auto tmp60 = tmp33 * tmp59;
    const auto tmp61 = tmp60 + tmp56;
    const auto tmp62 = 2 * tmp0[ 1 ];
    const auto tmp63 = tmp62 * tmp30;
    const auto tmp64 = -1 * tmp63;
    const auto tmp65 = -1 + tmp64;
    const auto tmp66 = tmp65 / tmp28;
    const auto tmp67 = 0.5625 * tmp66;
    const auto tmp68 = tmp67 * tmp25;
    const auto tmp69 = tmp58 * tmp42;
    const auto tmp70 = tmp69 * tmp31;
    const auto tmp71 = tmp70 + tmp68;
    const auto tmp72 = tmp71 * tmp8;
    const auto tmp73 = 0.28125 * tmp52;
    const auto tmp74 = tmp73 * tmp15;
    const auto tmp75 = tmp74 * tmp32;
    const auto tmp76 = tmp75 + tmp72;
    const auto tmp77 = tmp76 + tmp61;
    const auto tmp78 = tmp77 + tmp51;
    const auto tmp79 = tmp78 * (tmp12 == 0.0 ? 0.0 : (tmp12 < 0.0 ? -1 : 1));
    const auto tmp80 = -1 * ((tmp13[ 0 ])[ 1 ])[ 0 ];
    const auto tmp81 = tmp73 * tmp20;
    const auto tmp82 = tmp81 * tmp4;
    const auto tmp83 = tmp75 + tmp82;
    const auto tmp84 = tmp36 * tmp57;
    const auto tmp85 = -1 * tmp84;
    const auto tmp86 = 1 + tmp85;
    const auto tmp87 = tmp86 / tmp28;
    const auto tmp88 = 0.5625 * tmp87;
    const auto tmp89 = tmp88 * tmp25;
    const auto tmp90 = tmp43 * tmp58;
    const auto tmp91 = tmp90 + tmp89;
    const auto tmp92 = tmp91 * tmp8;
    const auto tmp93 = tmp60 + tmp92;
    const auto tmp94 = tmp93 + tmp83;
    const auto tmp95 = tmp94 + tmp80;
    const auto tmp96 = tmp95 * (tmp12 == 0.0 ? 0.0 : (tmp12 < 0.0 ? -1 : 1));
    const auto tmp97 = -1 * ((tmp13[ 0 ])[ 1 ])[ 1 ];
    const auto tmp98 = tmp73 * tmp54;
    const auto tmp99 = tmp98 + tmp16;
    const auto tmp100 = tmp99 * tmp4;
    const auto tmp101 = tmp74 * tmp59;
    const auto tmp102 = tmp101 + tmp100;
    const auto tmp103 = tmp62 * tmp57;
    const auto tmp104 = -1 * tmp103;
    const auto tmp105 = tmp104 / tmp28;
    const auto tmp106 = 0.5625 * tmp105;
    const auto tmp107 = tmp106 * tmp25;
    const auto tmp108 = tmp69 * tmp58;
    const auto tmp109 = tmp108 + tmp107;
    const auto tmp110 = tmp109 * tmp8;
    const auto tmp111 = tmp101 + tmp110;
    const auto tmp112 = tmp111 + tmp102;
    const auto tmp113 = tmp112 + tmp97;
    const auto tmp114 = tmp113 * (tmp12 == 0.0 ? 0.0 : (tmp12 < 0.0 ? -1 : 1));
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp50;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp79;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp96;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp114;
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
    std::get< i >( coefficients_ ).evaluate( x, result );;
    return result;
  }

  template< std::size_t i, class Point >
  typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
  {
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType result;
    std::get< i >( coefficients_ ).jacobian( x, result );;
    return result;
  }

  template< std::size_t i, class Point >
  typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
  {
    typename CoefficientFunctionSpaceType< i >::HessianRangeType result;
    std::get< i >( coefficients_ ).hessian( x, result );;
    return result;
  }
  ConstantTupleType constants_;
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution > > coefficients_;
};

} // namespace UFLLocalFunctions_a41df25f2b6abd5e9849dd885dd1b68a

PYBIND11_MODULE( localfunction_a41df25f2b6abd5e9849dd885dd1b68a_af122c1df944c95cd395ec0f91d0f970, module )
{
  typedef UFLLocalFunctions_a41df25f2b6abd5e9849dd885dd1b68a::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_a41df25f2b6abd5e9849dd885dd1b68a::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_a41df25f2b6abd5e9849dd885dd1b68a_af122c1df944c95cd395ec0f91d0f970.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > &coeffsolution ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>(gridView),name,order, coeffsolution); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
