#ifndef GUARD_a868408f9688e744b1aa50538fff55d6
#define GUARD_a868408f9688e744b1aa50538fff55d6

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

namespace UFLLocalFunctions_a868408f9688e744b1aa50538fff55d6
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffu_h >
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
  typedef std::tuple< Coeffu_h > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffu_h>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffu_h &coeffu_h, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffu_h >( coeffu_h ) )
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
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = std::pow( tmp0[ 0 ], 2 );
      assert(tmp1 == tmp1);
    const auto tmp2 = 3 * tmp1;
      assert(tmp2 == tmp2);
    const auto tmp3 = 0.3333333333333333 * tmp2;
      assert(tmp3 == tmp3);
    const auto tmp4 = -1 * tmp3;
      assert(tmp4 == tmp4);
    const auto tmp5 = 2 * tmp0[ 0 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = 0.5 * tmp5;
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 + tmp4;
      assert(tmp7 == tmp7);
    const auto tmp8 = -1 * tmp7;
      assert(tmp8 == tmp8);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp9 = jacobianCoefficient< 0 >( x );
      assert(tmp9 == tmp9);
    const auto tmp10 = (tmp9[ 0 ])[ 0 ] + tmp8;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * tmp10;
      assert(tmp11 == tmp11);
    const auto tmp12 = std::pow( tmp0[ 1 ], 2 );
      assert(tmp12 == tmp12);
    const auto tmp13 = 3 * tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = -1 * tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = 0.3333333333333333 * tmp14;
      assert(tmp15 == tmp15);
    const auto tmp16 = -1 * tmp15;
      assert(tmp16 == tmp16);
    const auto tmp17 = 2 * tmp0[ 1 ];
      assert(tmp17 == tmp17);
    const auto tmp18 = 0.5 * tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 + tmp16;
      assert(tmp19 == tmp19);
    const auto tmp20 = -1 * tmp19;
      assert(tmp20 == tmp20);
    const auto tmp21 = (tmp9[ 0 ])[ 1 ] + tmp20;
      assert(tmp21 == tmp21);
    const auto tmp22 = tmp21 * tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 + tmp11;
      assert(tmp23 == tmp23);
    result[ 0 ] = tmp23;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = std::pow( tmp0[ 0 ], 2 );
      assert(tmp1 == tmp1);
    const auto tmp2 = 3 * tmp1;
      assert(tmp2 == tmp2);
    const auto tmp3 = 0.3333333333333333 * tmp2;
      assert(tmp3 == tmp3);
    const auto tmp4 = -1 * tmp3;
      assert(tmp4 == tmp4);
    const auto tmp5 = 2 * tmp0[ 0 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = 0.5 * tmp5;
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 + tmp4;
      assert(tmp7 == tmp7);
    const auto tmp8 = -1 * tmp7;
      assert(tmp8 == tmp8);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp9 = jacobianCoefficient< 0 >( x );
      assert(tmp9 == tmp9);
    const auto tmp10 = (tmp9[ 0 ])[ 0 ] + tmp8;
      assert(tmp10 == tmp10);
    const auto tmp11 = 3 * tmp5;
      assert(tmp11 == tmp11);
    const auto tmp12 = 0.3333333333333333 * tmp11;
      assert(tmp12 == tmp12);
    const auto tmp13 = -1 * tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = 1.0 + tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = -1 * tmp14;
      assert(tmp15 == tmp15);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp16 = hessianCoefficient< 0 >( x );
      assert(tmp16 == tmp16);
    const auto tmp17 = ((tmp16[ 0 ])[ 0 ])[ 0 ] + tmp15;
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp17 * tmp10;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 + tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = std::pow( tmp0[ 1 ], 2 );
      assert(tmp20 == tmp20);
    const auto tmp21 = 3 * tmp20;
      assert(tmp21 == tmp21);
    const auto tmp22 = -1 * tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = 0.3333333333333333 * tmp22;
      assert(tmp23 == tmp23);
    const auto tmp24 = -1 * tmp23;
      assert(tmp24 == tmp24);
    const auto tmp25 = 2 * tmp0[ 1 ];
      assert(tmp25 == tmp25);
    const auto tmp26 = 0.5 * tmp25;
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp26 + tmp24;
      assert(tmp27 == tmp27);
    const auto tmp28 = -1 * tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = (tmp9[ 0 ])[ 1 ] + tmp28;
      assert(tmp29 == tmp29);
    const auto tmp30 = ((tmp16[ 0 ])[ 1 ])[ 0 ] * tmp29;
      assert(tmp30 == tmp30);
    const auto tmp31 = tmp30 + tmp30;
      assert(tmp31 == tmp31);
    const auto tmp32 = tmp31 + tmp19;
      assert(tmp32 == tmp32);
    const auto tmp33 = ((tmp16[ 0 ])[ 0 ])[ 1 ] * tmp10;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 + tmp33;
      assert(tmp34 == tmp34);
    const auto tmp35 = 3 * tmp25;
      assert(tmp35 == tmp35);
    const auto tmp36 = -1 * tmp35;
      assert(tmp36 == tmp36);
    const auto tmp37 = 0.3333333333333333 * tmp36;
      assert(tmp37 == tmp37);
    const auto tmp38 = -1 * tmp37;
      assert(tmp38 == tmp38);
    const auto tmp39 = 1.0 + tmp38;
      assert(tmp39 == tmp39);
    const auto tmp40 = -1 * tmp39;
      assert(tmp40 == tmp40);
    const auto tmp41 = ((tmp16[ 0 ])[ 1 ])[ 1 ] + tmp40;
      assert(tmp41 == tmp41);
    const auto tmp42 = tmp41 * tmp29;
      assert(tmp42 == tmp42);
    const auto tmp43 = tmp42 + tmp42;
      assert(tmp43 == tmp43);
    const auto tmp44 = tmp43 + tmp34;
      assert(tmp44 == tmp44);
    (result[ 0 ])[ 0 ] = tmp32;
    (result[ 0 ])[ 1 ] = tmp44;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffu_h > > coefficients_;
};

} // namespace UFLLocalFunctions_a868408f9688e744b1aa50538fff55d6

PYBIND11_MODULE( localfunction_a868408f9688e744b1aa50538fff55d6_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_a868408f9688e744b1aa50538fff55d6::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_a868408f9688e744b1aa50538fff55d6::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_a868408f9688e744b1aa50538fff55d6_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffu_h ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffu_h); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
