#ifndef GUARD_c5753a4079884ed9a6778e98a54f18ad
#define GUARD_c5753a4079884ed9a6778e98a54f18ad

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

namespace UFLLocalFunctions_c5753a4079884ed9a6778e98a54f18ad
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
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = std::pow( tmp0[ 0 ], 3 );
    const auto tmp2 = std::pow( tmp0[ 1 ], 3 );
    const auto tmp3 = -1 * tmp2;
    const auto tmp4 = tmp3 + tmp1;
    const auto tmp5 = 0.3333333333333333 * tmp4;
    const auto tmp6 = -1 * tmp5;
    const auto tmp7 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp8 = std::pow( tmp0[ 0 ], 2 );
    const auto tmp9 = tmp8 + tmp7;
    const auto tmp10 = 0.5 * tmp9;
    const auto tmp11 = tmp10 + tmp6;
    const auto tmp12 = 1 + tmp11;
    const auto tmp13 = -1 * tmp12;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp14 = evaluateCoefficient< 0 >( x );
    const auto tmp15 = tmp14[ 0 ] + tmp13;
    const auto tmp16 = tmp15 * tmp15;
    const auto tmp17 = 3 * tmp8;
    const auto tmp18 = 0.3333333333333333 * tmp17;
    const auto tmp19 = -1 * tmp18;
    const auto tmp20 = 2 * tmp0[ 0 ];
    const auto tmp21 = 0.5 * tmp20;
    const auto tmp22 = tmp21 + tmp19;
    const auto tmp23 = -1 * tmp22;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp24 = jacobianCoefficient< 0 >( x );
    const auto tmp25 = (tmp24[ 0 ])[ 0 ] + tmp23;
    const auto tmp26 = tmp25 * tmp25;
    const auto tmp27 = 3 * tmp7;
    const auto tmp28 = -1 * tmp27;
    const auto tmp29 = 0.3333333333333333 * tmp28;
    const auto tmp30 = -1 * tmp29;
    const auto tmp31 = 2 * tmp0[ 1 ];
    const auto tmp32 = 0.5 * tmp31;
    const auto tmp33 = tmp32 + tmp30;
    const auto tmp34 = -1 * tmp33;
    const auto tmp35 = (tmp24[ 0 ])[ 1 ] + tmp34;
    const auto tmp36 = tmp35 * tmp35;
    const auto tmp37 = tmp36 + tmp26;
    result[ 0 ] = tmp16;
    result[ 1 ] = tmp37;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = std::pow( tmp0[ 0 ], 2 );
    const auto tmp2 = 3 * tmp1;
    const auto tmp3 = 0.3333333333333333 * tmp2;
    const auto tmp4 = -1 * tmp3;
    const auto tmp5 = 2 * tmp0[ 0 ];
    const auto tmp6 = 0.5 * tmp5;
    const auto tmp7 = tmp6 + tmp4;
    const auto tmp8 = -1 * tmp7;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp9 = jacobianCoefficient< 0 >( x );
    const auto tmp10 = (tmp9[ 0 ])[ 0 ] + tmp8;
    const auto tmp11 = std::pow( tmp0[ 0 ], 3 );
    const auto tmp12 = std::pow( tmp0[ 1 ], 3 );
    const auto tmp13 = -1 * tmp12;
    const auto tmp14 = tmp13 + tmp11;
    const auto tmp15 = 0.3333333333333333 * tmp14;
    const auto tmp16 = -1 * tmp15;
    const auto tmp17 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp18 = tmp1 + tmp17;
    const auto tmp19 = 0.5 * tmp18;
    const auto tmp20 = tmp19 + tmp16;
    const auto tmp21 = 1 + tmp20;
    const auto tmp22 = -1 * tmp21;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp23 = evaluateCoefficient< 0 >( x );
    const auto tmp24 = tmp23[ 0 ] + tmp22;
    const auto tmp25 = tmp24 * tmp10;
    const auto tmp26 = tmp25 + tmp25;
    const auto tmp27 = 3 * tmp17;
    const auto tmp28 = -1 * tmp27;
    const auto tmp29 = 0.3333333333333333 * tmp28;
    const auto tmp30 = -1 * tmp29;
    const auto tmp31 = 2 * tmp0[ 1 ];
    const auto tmp32 = 0.5 * tmp31;
    const auto tmp33 = tmp32 + tmp30;
    const auto tmp34 = -1 * tmp33;
    const auto tmp35 = (tmp9[ 0 ])[ 1 ] + tmp34;
    const auto tmp36 = tmp24 * tmp35;
    const auto tmp37 = tmp36 + tmp36;
    const auto tmp38 = 3 * tmp5;
    const auto tmp39 = 0.3333333333333333 * tmp38;
    const auto tmp40 = -1 * tmp39;
    const auto tmp41 = 1.0 + tmp40;
    const auto tmp42 = -1 * tmp41;
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp43 = hessianCoefficient< 0 >( x );
    const auto tmp44 = ((tmp43[ 0 ])[ 0 ])[ 0 ] + tmp42;
    const auto tmp45 = tmp44 * tmp10;
    const auto tmp46 = tmp45 + tmp45;
    const auto tmp47 = ((tmp43[ 0 ])[ 1 ])[ 0 ] * tmp35;
    const auto tmp48 = tmp47 + tmp47;
    const auto tmp49 = tmp48 + tmp46;
    const auto tmp50 = ((tmp43[ 0 ])[ 0 ])[ 1 ] * tmp10;
    const auto tmp51 = tmp50 + tmp50;
    const auto tmp52 = 3 * tmp31;
    const auto tmp53 = -1 * tmp52;
    const auto tmp54 = 0.3333333333333333 * tmp53;
    const auto tmp55 = -1 * tmp54;
    const auto tmp56 = 1.0 + tmp55;
    const auto tmp57 = -1 * tmp56;
    const auto tmp58 = ((tmp43[ 0 ])[ 1 ])[ 1 ] + tmp57;
    const auto tmp59 = tmp58 * tmp35;
    const auto tmp60 = tmp59 + tmp59;
    const auto tmp61 = tmp60 + tmp51;
    (result[ 0 ])[ 0 ] = tmp26;
    (result[ 0 ])[ 1 ] = tmp37;
    (result[ 1 ])[ 0 ] = tmp49;
    (result[ 1 ])[ 1 ] = tmp61;
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

} // namespace UFLLocalFunctions_c5753a4079884ed9a6778e98a54f18ad

PYBIND11_MODULE( localfunction_c5753a4079884ed9a6778e98a54f18ad_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_c5753a4079884ed9a6778e98a54f18ad::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_c5753a4079884ed9a6778e98a54f18ad::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_c5753a4079884ed9a6778e98a54f18ad_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffu_h ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffu_h); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
