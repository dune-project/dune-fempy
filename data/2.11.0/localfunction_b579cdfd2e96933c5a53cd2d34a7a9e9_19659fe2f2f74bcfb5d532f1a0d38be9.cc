#ifndef GUARD_b579cdfd2e96933c5a53cd2d34a7a9e9
#define GUARD_b579cdfd2e96933c5a53cd2d34a7a9e9

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

namespace UFLLocalFunctions_b579cdfd2e96933c5a53cd2d34a7a9e9
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffu_h >
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
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 6.283185307179586 * tmp0[ 1 ];
    const auto tmp2 = std::cos( tmp1 );
    const auto tmp3 = 6.283185307179586 * tmp0[ 0 ];
    const auto tmp4 = std::cos( tmp3 );
    const auto tmp5 = tmp4 * tmp2;
    const auto tmp6 = -1 * tmp5;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp7 = evaluateCoefficient< 0 >( x );
    const auto tmp8 = tmp7[ 0 ] + tmp6;
    const auto tmp9 = std::pow( tmp8, 2 );
    const auto tmp10 = std::sin( tmp3 );
    const auto tmp11 = -1 * tmp10;
    const auto tmp12 = 6.283185307179586 * tmp11;
    const auto tmp13 = tmp12 * tmp2;
    const auto tmp14 = -1 * tmp13;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp15 = jacobianCoefficient< 0 >( x );
    const auto tmp16 = (tmp15[ 0 ])[ 0 ] + tmp14;
    const auto tmp17 = tmp16 * tmp16;
    const auto tmp18 = std::sin( tmp1 );
    const auto tmp19 = -1 * tmp18;
    const auto tmp20 = 6.283185307179586 * tmp19;
    const auto tmp21 = tmp20 * tmp4;
    const auto tmp22 = -1 * tmp21;
    const auto tmp23 = (tmp15[ 0 ])[ 1 ] + tmp22;
    const auto tmp24 = tmp23 * tmp23;
    const auto tmp25 = tmp24 + tmp17;
    result[ 0 ] = tmp9;
    result[ 1 ] = tmp25;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 6.283185307179586 * tmp0[ 1 ];
    const auto tmp2 = std::cos( tmp1 );
    const auto tmp3 = 6.283185307179586 * tmp0[ 0 ];
    const auto tmp4 = std::sin( tmp3 );
    const auto tmp5 = -1 * tmp4;
    const auto tmp6 = 6.283185307179586 * tmp5;
    const auto tmp7 = tmp6 * tmp2;
    const auto tmp8 = -1 * tmp7;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp9 = jacobianCoefficient< 0 >( x );
    const auto tmp10 = (tmp9[ 0 ])[ 0 ] + tmp8;
    const auto tmp11 = 2 * tmp10;
    const auto tmp12 = std::cos( tmp3 );
    const auto tmp13 = tmp12 * tmp2;
    const auto tmp14 = -1 * tmp13;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp15 = evaluateCoefficient< 0 >( x );
    const auto tmp16 = tmp15[ 0 ] + tmp14;
    const auto tmp17 = tmp16 * tmp11;
    const auto tmp18 = std::sin( tmp1 );
    const auto tmp19 = -1 * tmp18;
    const auto tmp20 = 6.283185307179586 * tmp19;
    const auto tmp21 = tmp20 * tmp12;
    const auto tmp22 = -1 * tmp21;
    const auto tmp23 = (tmp9[ 0 ])[ 1 ] + tmp22;
    const auto tmp24 = 2 * tmp23;
    const auto tmp25 = tmp16 * tmp24;
    const auto tmp26 = 6.283185307179586 * tmp12;
    const auto tmp27 = -1 * tmp26;
    const auto tmp28 = 6.283185307179586 * tmp27;
    const auto tmp29 = tmp28 * tmp2;
    const auto tmp30 = -1 * tmp29;
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp31 = hessianCoefficient< 0 >( x );
    const auto tmp32 = ((tmp31[ 0 ])[ 0 ])[ 0 ] + tmp30;
    const auto tmp33 = tmp32 * tmp10;
    const auto tmp34 = tmp10 * tmp32;
    const auto tmp35 = tmp34 + tmp33;
    const auto tmp36 = tmp6 * tmp20;
    const auto tmp37 = -1 * tmp36;
    const auto tmp38 = ((tmp31[ 0 ])[ 1 ])[ 0 ] + tmp37;
    const auto tmp39 = tmp38 * tmp23;
    const auto tmp40 = tmp23 * tmp38;
    const auto tmp41 = tmp40 + tmp39;
    const auto tmp42 = tmp41 + tmp35;
    const auto tmp43 = ((tmp31[ 0 ])[ 0 ])[ 1 ] + tmp37;
    const auto tmp44 = tmp43 * tmp10;
    const auto tmp45 = tmp10 * tmp43;
    const auto tmp46 = tmp45 + tmp44;
    const auto tmp47 = 6.283185307179586 * tmp2;
    const auto tmp48 = -1 * tmp47;
    const auto tmp49 = 6.283185307179586 * tmp48;
    const auto tmp50 = tmp49 * tmp12;
    const auto tmp51 = -1 * tmp50;
    const auto tmp52 = ((tmp31[ 0 ])[ 1 ])[ 1 ] + tmp51;
    const auto tmp53 = tmp52 * tmp23;
    const auto tmp54 = tmp23 * tmp52;
    const auto tmp55 = tmp54 + tmp53;
    const auto tmp56 = tmp55 + tmp46;
    (result[ 0 ])[ 0 ] = tmp17;
    (result[ 0 ])[ 1 ] = tmp25;
    (result[ 1 ])[ 0 ] = tmp42;
    (result[ 1 ])[ 1 ] = tmp56;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffu_h > > coefficients_;
};

} // namespace UFLLocalFunctions_b579cdfd2e96933c5a53cd2d34a7a9e9

PYBIND11_MODULE( localfunction_b579cdfd2e96933c5a53cd2d34a7a9e9_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_b579cdfd2e96933c5a53cd2d34a7a9e9::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_b579cdfd2e96933c5a53cd2d34a7a9e9::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_b579cdfd2e96933c5a53cd2d34a7a9e9_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffu_h ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffu_h); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
