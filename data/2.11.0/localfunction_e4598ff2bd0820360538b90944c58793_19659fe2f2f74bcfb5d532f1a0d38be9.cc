#ifndef GUARD_e4598ff2bd0820360538b90944c58793
#define GUARD_e4598ff2bd0820360538b90944c58793

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

namespace UFLLocalFunctions_e4598ff2bd0820360538b90944c58793
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
    const auto tmp13 = -1 + tmp12;
    const auto tmp14 = 0.6065306597126334 * tmp13;
    const auto tmp15 = 1 + tmp14;
    const auto tmp16 = -1 * tmp15;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp17 = evaluateCoefficient< 0 >( x );
    const auto tmp18 = tmp17[ 0 ] + tmp16;
    const auto tmp19 = tmp18 * tmp18;
    result[ 0 ] = tmp19;
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
    const auto tmp8 = 0.6065306597126334 * tmp7;
    const auto tmp9 = -1 * tmp8;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp10 = jacobianCoefficient< 0 >( x );
    const auto tmp11 = (tmp10[ 0 ])[ 0 ] + tmp9;
    const auto tmp12 = std::pow( tmp0[ 0 ], 3 );
    const auto tmp13 = std::pow( tmp0[ 1 ], 3 );
    const auto tmp14 = -1 * tmp13;
    const auto tmp15 = tmp14 + tmp12;
    const auto tmp16 = 0.3333333333333333 * tmp15;
    const auto tmp17 = -1 * tmp16;
    const auto tmp18 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp19 = tmp1 + tmp18;
    const auto tmp20 = 0.5 * tmp19;
    const auto tmp21 = tmp20 + tmp17;
    const auto tmp22 = 1 + tmp21;
    const auto tmp23 = -1 + tmp22;
    const auto tmp24 = 0.6065306597126334 * tmp23;
    const auto tmp25 = 1 + tmp24;
    const auto tmp26 = -1 * tmp25;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp27 = evaluateCoefficient< 0 >( x );
    const auto tmp28 = tmp27[ 0 ] + tmp26;
    const auto tmp29 = tmp28 * tmp11;
    const auto tmp30 = tmp29 + tmp29;
    const auto tmp31 = 3 * tmp18;
    const auto tmp32 = -1 * tmp31;
    const auto tmp33 = 0.3333333333333333 * tmp32;
    const auto tmp34 = -1 * tmp33;
    const auto tmp35 = 2 * tmp0[ 1 ];
    const auto tmp36 = 0.5 * tmp35;
    const auto tmp37 = tmp36 + tmp34;
    const auto tmp38 = 0.6065306597126334 * tmp37;
    const auto tmp39 = -1 * tmp38;
    const auto tmp40 = (tmp10[ 0 ])[ 1 ] + tmp39;
    const auto tmp41 = tmp28 * tmp40;
    const auto tmp42 = tmp41 + tmp41;
    (result[ 0 ])[ 0 ] = tmp30;
    (result[ 0 ])[ 1 ] = tmp42;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
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
    const auto tmp8 = 0.6065306597126334 * tmp7;
    const auto tmp9 = -1 * tmp8;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp10 = jacobianCoefficient< 0 >( x );
    const auto tmp11 = (tmp10[ 0 ])[ 0 ] + tmp9;
    const auto tmp12 = tmp11 * tmp11;
    const auto tmp13 = 3 * tmp5;
    const auto tmp14 = 0.3333333333333333 * tmp13;
    const auto tmp15 = -1 * tmp14;
    const auto tmp16 = 1.0 + tmp15;
    const auto tmp17 = 0.6065306597126334 * tmp16;
    const auto tmp18 = -1 * tmp17;
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp19 = hessianCoefficient< 0 >( x );
    const auto tmp20 = ((tmp19[ 0 ])[ 0 ])[ 0 ] + tmp18;
    const auto tmp21 = std::pow( tmp0[ 0 ], 3 );
    const auto tmp22 = std::pow( tmp0[ 1 ], 3 );
    const auto tmp23 = -1 * tmp22;
    const auto tmp24 = tmp23 + tmp21;
    const auto tmp25 = 0.3333333333333333 * tmp24;
    const auto tmp26 = -1 * tmp25;
    const auto tmp27 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp28 = tmp1 + tmp27;
    const auto tmp29 = 0.5 * tmp28;
    const auto tmp30 = tmp29 + tmp26;
    const auto tmp31 = 1 + tmp30;
    const auto tmp32 = -1 + tmp31;
    const auto tmp33 = 0.6065306597126334 * tmp32;
    const auto tmp34 = 1 + tmp33;
    const auto tmp35 = -1 * tmp34;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp36 = evaluateCoefficient< 0 >( x );
    const auto tmp37 = tmp36[ 0 ] + tmp35;
    const auto tmp38 = tmp37 * tmp20;
    const auto tmp39 = tmp38 + tmp12;
    const auto tmp40 = tmp39 + tmp39;
    const auto tmp41 = 3 * tmp27;
    const auto tmp42 = -1 * tmp41;
    const auto tmp43 = 0.3333333333333333 * tmp42;
    const auto tmp44 = -1 * tmp43;
    const auto tmp45 = 2 * tmp0[ 1 ];
    const auto tmp46 = 0.5 * tmp45;
    const auto tmp47 = tmp46 + tmp44;
    const auto tmp48 = 0.6065306597126334 * tmp47;
    const auto tmp49 = -1 * tmp48;
    const auto tmp50 = (tmp10[ 0 ])[ 1 ] + tmp49;
    const auto tmp51 = tmp50 * tmp11;
    const auto tmp52 = ((tmp19[ 0 ])[ 0 ])[ 1 ] * tmp37;
    const auto tmp53 = tmp52 + tmp51;
    const auto tmp54 = tmp53 + tmp53;
    const auto tmp55 = ((tmp19[ 0 ])[ 1 ])[ 0 ] * tmp37;
    const auto tmp56 = tmp55 + tmp51;
    const auto tmp57 = tmp56 + tmp56;
    const auto tmp58 = tmp50 * tmp50;
    const auto tmp59 = 3 * tmp45;
    const auto tmp60 = -1 * tmp59;
    const auto tmp61 = 0.3333333333333333 * tmp60;
    const auto tmp62 = -1 * tmp61;
    const auto tmp63 = 1.0 + tmp62;
    const auto tmp64 = 0.6065306597126334 * tmp63;
    const auto tmp65 = -1 * tmp64;
    const auto tmp66 = ((tmp19[ 0 ])[ 1 ])[ 1 ] + tmp65;
    const auto tmp67 = tmp37 * tmp66;
    const auto tmp68 = tmp67 + tmp58;
    const auto tmp69 = tmp68 + tmp68;
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp40;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp54;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp57;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp69;
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

} // namespace UFLLocalFunctions_e4598ff2bd0820360538b90944c58793

PYBIND11_MODULE( localfunction_e4598ff2bd0820360538b90944c58793_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_e4598ff2bd0820360538b90944c58793::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_e4598ff2bd0820360538b90944c58793::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_e4598ff2bd0820360538b90944c58793_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffu_h ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffu_h); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
