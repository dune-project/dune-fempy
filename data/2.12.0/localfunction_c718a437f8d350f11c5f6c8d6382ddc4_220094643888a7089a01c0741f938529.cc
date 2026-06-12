#ifndef GUARD_c718a437f8d350f11c5f6c8d6382ddc4
#define GUARD_c718a437f8d350f11c5f6c8d6382ddc4

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

namespace UFLLocalFunctions_c718a437f8d350f11c5f6c8d6382ddc4
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffsolution, class Coeffiterate >
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
  typedef std::tuple<> ConstantTupleType;
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< Coeffsolution, Coeffiterate > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsolution>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffiterate>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffsolution &coeffsolution, const Coeffiterate &coeffiterate, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffsolution >( coeffsolution ), Dune::Fem::ConstLocalFunction< Coeffiterate >( coeffiterate ) )
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
    using std::pow;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = std::pow( tmp0[ 1 ], 2 );
      assert(tmp1 == tmp1);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp2 = evaluateCoefficient< 1 >( x );
      assert(tmp2 == tmp2);
    const auto tmp3 = -1 * tmp2[ 1 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = tmp0[ 1 ] + tmp3;
      assert(tmp4 == tmp4);
    const auto tmp5 = std::pow( tmp4, 2 );
      assert(tmp5 == tmp5);
    result[ 0 ] = tmp1;
    result[ 1 ] = tmp5;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp0 = jacobianCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = 2 * (tmp0[ 1 ])[ 0 ];
      assert(tmp1 == tmp1);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp2 = evaluateCoefficient< 0 >( x );
      assert(tmp2 == tmp2);
    const auto tmp3 = tmp2[ 1 ] * tmp1;
      assert(tmp3 == tmp3);
    const auto tmp4 = 2 * (tmp0[ 1 ])[ 1 ];
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp2[ 1 ] * tmp4;
      assert(tmp5 == tmp5);
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp6 = jacobianCoefficient< 1 >( x );
      assert(tmp6 == tmp6);
    const auto tmp7 = -1 * (tmp6[ 1 ])[ 0 ];
      assert(tmp7 == tmp7);
    const auto tmp8 = (tmp0[ 1 ])[ 0 ] + tmp7;
      assert(tmp8 == tmp8);
    const auto tmp9 = 2 * tmp8;
      assert(tmp9 == tmp9);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp10 = evaluateCoefficient< 1 >( x );
      assert(tmp10 == tmp10);
    const auto tmp11 = -1 * tmp10[ 1 ];
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp2[ 1 ] + tmp11;
      assert(tmp12 == tmp12);
    const auto tmp13 = tmp12 * tmp9;
      assert(tmp13 == tmp13);
    const auto tmp14 = -1 * (tmp6[ 1 ])[ 1 ];
      assert(tmp14 == tmp14);
    const auto tmp15 = (tmp0[ 1 ])[ 1 ] + tmp14;
      assert(tmp15 == tmp15);
    const auto tmp16 = 2 * tmp15;
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp12 * tmp16;
      assert(tmp17 == tmp17);
    (result[ 0 ])[ 0 ] = tmp3;
    (result[ 0 ])[ 1 ] = tmp5;
    (result[ 1 ])[ 0 ] = tmp13;
    (result[ 1 ])[ 1 ] = tmp17;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp0 = hessianCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = 2 * ((tmp0[ 1 ])[ 0 ])[ 0 ];
      assert(tmp1 == tmp1);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp2 = evaluateCoefficient< 0 >( x );
      assert(tmp2 == tmp2);
    const auto tmp3 = tmp2[ 1 ] * tmp1;
      assert(tmp3 == tmp3);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp4 = jacobianCoefficient< 0 >( x );
      assert(tmp4 == tmp4);
    const auto tmp5 = 2 * (tmp4[ 1 ])[ 0 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = (tmp4[ 1 ])[ 0 ] * tmp5;
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 + tmp3;
      assert(tmp7 == tmp7);
    const auto tmp8 = 2 * ((tmp0[ 1 ])[ 0 ])[ 1 ];
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp2[ 1 ] * tmp8;
      assert(tmp9 == tmp9);
    const auto tmp10 = (tmp4[ 1 ])[ 1 ] * tmp5;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 + tmp9;
      assert(tmp11 == tmp11);
    const auto tmp12 = 2 * ((tmp0[ 1 ])[ 1 ])[ 0 ];
      assert(tmp12 == tmp12);
    const auto tmp13 = tmp2[ 1 ] * tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = 2 * (tmp4[ 1 ])[ 1 ];
      assert(tmp14 == tmp14);
    const auto tmp15 = (tmp4[ 1 ])[ 0 ] * tmp14;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 + tmp13;
      assert(tmp16 == tmp16);
    const auto tmp17 = 2 * ((tmp0[ 1 ])[ 1 ])[ 1 ];
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp2[ 1 ] * tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = (tmp4[ 1 ])[ 1 ] * tmp14;
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp19 + tmp18;
      assert(tmp20 == tmp20);
    typename CoefficientFunctionSpaceType< 1 >::HessianRangeType tmp21 = hessianCoefficient< 1 >( x );
      assert(tmp21 == tmp21);
    const auto tmp22 = -1 * ((tmp21[ 1 ])[ 0 ])[ 0 ];
      assert(tmp22 == tmp22);
    const auto tmp23 = ((tmp0[ 1 ])[ 0 ])[ 0 ] + tmp22;
      assert(tmp23 == tmp23);
    const auto tmp24 = 2 * tmp23;
      assert(tmp24 == tmp24);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp25 = evaluateCoefficient< 1 >( x );
      assert(tmp25 == tmp25);
    const auto tmp26 = -1 * tmp25[ 1 ];
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp2[ 1 ] + tmp26;
      assert(tmp27 == tmp27);
    const auto tmp28 = tmp27 * tmp24;
      assert(tmp28 == tmp28);
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp29 = jacobianCoefficient< 1 >( x );
      assert(tmp29 == tmp29);
    const auto tmp30 = -1 * (tmp29[ 1 ])[ 0 ];
      assert(tmp30 == tmp30);
    const auto tmp31 = (tmp4[ 1 ])[ 0 ] + tmp30;
      assert(tmp31 == tmp31);
    const auto tmp32 = 2 * tmp31;
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp31 * tmp32;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 + tmp28;
      assert(tmp34 == tmp34);
    const auto tmp35 = -1 * ((tmp21[ 1 ])[ 0 ])[ 1 ];
      assert(tmp35 == tmp35);
    const auto tmp36 = ((tmp0[ 1 ])[ 0 ])[ 1 ] + tmp35;
      assert(tmp36 == tmp36);
    const auto tmp37 = 2 * tmp36;
      assert(tmp37 == tmp37);
    const auto tmp38 = tmp27 * tmp37;
      assert(tmp38 == tmp38);
    const auto tmp39 = -1 * (tmp29[ 1 ])[ 1 ];
      assert(tmp39 == tmp39);
    const auto tmp40 = (tmp4[ 1 ])[ 1 ] + tmp39;
      assert(tmp40 == tmp40);
    const auto tmp41 = tmp40 * tmp32;
      assert(tmp41 == tmp41);
    const auto tmp42 = tmp41 + tmp38;
      assert(tmp42 == tmp42);
    const auto tmp43 = -1 * ((tmp21[ 1 ])[ 1 ])[ 0 ];
      assert(tmp43 == tmp43);
    const auto tmp44 = ((tmp0[ 1 ])[ 1 ])[ 0 ] + tmp43;
      assert(tmp44 == tmp44);
    const auto tmp45 = 2 * tmp44;
      assert(tmp45 == tmp45);
    const auto tmp46 = tmp27 * tmp45;
      assert(tmp46 == tmp46);
    const auto tmp47 = 2 * tmp40;
      assert(tmp47 == tmp47);
    const auto tmp48 = tmp31 * tmp47;
      assert(tmp48 == tmp48);
    const auto tmp49 = tmp48 + tmp46;
      assert(tmp49 == tmp49);
    const auto tmp50 = -1 * ((tmp21[ 1 ])[ 1 ])[ 1 ];
      assert(tmp50 == tmp50);
    const auto tmp51 = ((tmp0[ 1 ])[ 1 ])[ 1 ] + tmp50;
      assert(tmp51 == tmp51);
    const auto tmp52 = 2 * tmp51;
      assert(tmp52 == tmp52);
    const auto tmp53 = tmp27 * tmp52;
      assert(tmp53 == tmp53);
    const auto tmp54 = tmp40 * tmp47;
      assert(tmp54 == tmp54);
    const auto tmp55 = tmp54 + tmp53;
      assert(tmp55 == tmp55);
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp7;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp11;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp16;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp20;
    ((result[ 1 ])[ 0 ])[ 0 ] = tmp34;
    ((result[ 1 ])[ 0 ])[ 1 ] = tmp42;
    ((result[ 1 ])[ 1 ])[ 0 ] = tmp49;
    ((result[ 1 ])[ 1 ])[ 1 ] = tmp55;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution >, Dune::Fem::ConstLocalFunction< Coeffiterate > > coefficients_;
};

} // namespace UFLLocalFunctions_c718a437f8d350f11c5f6c8d6382ddc4

PYBIND11_MODULE( localfunction_c718a437f8d350f11c5f6c8d6382ddc4_220094643888a7089a01c0741f938529, module )
{
  typedef UFLLocalFunctions_c718a437f8d350f11c5f6c8d6382ddc4::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_c718a437f8d350f11c5f6c8d6382ddc4::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_c718a437f8d350f11c5f6c8d6382ddc4_220094643888a7089a01c0741f938529.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffsolution, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>, Dune::FieldVector< double, 2 > > &coeffiterate ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::cube > >>(gridView),name,order, coeffsolution, coeffiterate); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
