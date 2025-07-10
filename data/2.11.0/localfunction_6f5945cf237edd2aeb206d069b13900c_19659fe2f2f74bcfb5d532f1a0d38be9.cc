#ifndef GUARD_6f5945cf237edd2aeb206d069b13900c
#define GUARD_6f5945cf237edd2aeb206d069b13900c

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

namespace UFLLocalFunctions_6f5945cf237edd2aeb206d069b13900c
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffcopy_of_u, class Coeffcopy_of_v >
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
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< Coeffcopy_of_u, Coeffcopy_of_v > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffcopy_of_u>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffcopy_of_v>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffcopy_of_u &coeffcopy_of_u, const Coeffcopy_of_v &coeffcopy_of_v, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffcopy_of_u >( coeffcopy_of_u ), Dune::Fem::ConstLocalFunction< Coeffcopy_of_v >( coeffcopy_of_v ) )
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
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp0 = evaluateCoefficient< 1 >( x );
    const auto tmp1 = -1 * tmp0[ 0 ];
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp2 = evaluateCoefficient< 0 >( x );
    const auto tmp3 = tmp2[ 0 ] + tmp1;
    const auto tmp4 = 0.25 * tmp3;
    const auto tmp5 = tmp0[ 0 ] + tmp4;
    result[ 0 ] = tmp5;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp0 = jacobianCoefficient< 1 >( x );
    const auto tmp1 = -1 * (tmp0[ 0 ])[ 0 ];
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp2 = jacobianCoefficient< 0 >( x );
    const auto tmp3 = (tmp2[ 0 ])[ 0 ] + tmp1;
    const auto tmp4 = 0.25 * tmp3;
    const auto tmp5 = (tmp0[ 0 ])[ 0 ] + tmp4;
    const auto tmp6 = -1 * (tmp0[ 0 ])[ 1 ];
    const auto tmp7 = (tmp2[ 0 ])[ 1 ] + tmp6;
    const auto tmp8 = 0.25 * tmp7;
    const auto tmp9 = (tmp0[ 0 ])[ 1 ] + tmp8;
    (result[ 0 ])[ 0 ] = tmp5;
    (result[ 0 ])[ 1 ] = tmp9;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 1 >::HessianRangeType tmp0 = hessianCoefficient< 1 >( x );
    const auto tmp1 = -1 * ((tmp0[ 0 ])[ 0 ])[ 0 ];
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp2 = hessianCoefficient< 0 >( x );
    const auto tmp3 = ((tmp2[ 0 ])[ 0 ])[ 0 ] + tmp1;
    const auto tmp4 = 0.25 * tmp3;
    const auto tmp5 = ((tmp0[ 0 ])[ 0 ])[ 0 ] + tmp4;
    const auto tmp6 = -1 * ((tmp0[ 0 ])[ 0 ])[ 1 ];
    const auto tmp7 = ((tmp2[ 0 ])[ 0 ])[ 1 ] + tmp6;
    const auto tmp8 = 0.25 * tmp7;
    const auto tmp9 = ((tmp0[ 0 ])[ 0 ])[ 1 ] + tmp8;
    const auto tmp10 = -1 * ((tmp0[ 0 ])[ 1 ])[ 0 ];
    const auto tmp11 = ((tmp2[ 0 ])[ 1 ])[ 0 ] + tmp10;
    const auto tmp12 = 0.25 * tmp11;
    const auto tmp13 = ((tmp0[ 0 ])[ 1 ])[ 0 ] + tmp12;
    const auto tmp14 = -1 * ((tmp0[ 0 ])[ 1 ])[ 1 ];
    const auto tmp15 = ((tmp2[ 0 ])[ 1 ])[ 1 ] + tmp14;
    const auto tmp16 = 0.25 * tmp15;
    const auto tmp17 = ((tmp0[ 0 ])[ 1 ])[ 1 ] + tmp16;
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp5;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp9;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp13;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp17;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffcopy_of_u >, Dune::Fem::ConstLocalFunction< Coeffcopy_of_v > > coefficients_;
};

} // namespace UFLLocalFunctions_6f5945cf237edd2aeb206d069b13900c

PYBIND11_MODULE( localfunction_6f5945cf237edd2aeb206d069b13900c_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_6f5945cf237edd2aeb206d069b13900c::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_6f5945cf237edd2aeb206d069b13900c::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_6f5945cf237edd2aeb206d069b13900c_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffcopy_of_u, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffcopy_of_v ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffcopy_of_u, coeffcopy_of_v); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
