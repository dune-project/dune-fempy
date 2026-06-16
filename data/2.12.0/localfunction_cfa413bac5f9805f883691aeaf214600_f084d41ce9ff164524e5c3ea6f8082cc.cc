#ifndef GUARD_cfa413bac5f9805f883691aeaf214600
#define GUARD_cfa413bac5f9805f883691aeaf214600

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/common/exceptions.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_cfa413bac5f9805f883691aeaf214600
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffdisplacement >
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
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< Coeffdisplacement > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffdisplacement>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffdisplacement &coeffdisplacement, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffdisplacement >( coeffdisplacement ) )
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
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    GlobalCoordinateType tmp1 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp1 == tmp1);
    const auto tmp2 = tmp1[ 0 ] + tmp0[ 0 ];
      assert(tmp2 == tmp2);
    const auto tmp3 = tmp1[ 1 ] + tmp0[ 1 ];
      assert(tmp3 == tmp3);
    result[ 0 ] = tmp2;
    result[ 1 ] = tmp3;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp0 = jacobianCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = 1 + (tmp0[ 0 ])[ 0 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = 1 + (tmp0[ 1 ])[ 1 ];
      assert(tmp2 == tmp2);
    (result[ 0 ])[ 0 ] = tmp1;
    (result[ 0 ])[ 1 ] = (tmp0[ 0 ])[ 1 ];
    (result[ 1 ])[ 0 ] = (tmp0[ 1 ])[ 0 ];
    (result[ 1 ])[ 1 ] = tmp2;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp0 = hessianCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    ((result[ 0 ])[ 0 ])[ 0 ] = ((tmp0[ 0 ])[ 0 ])[ 0 ];
    ((result[ 0 ])[ 0 ])[ 1 ] = ((tmp0[ 0 ])[ 0 ])[ 1 ];
    ((result[ 0 ])[ 1 ])[ 0 ] = ((tmp0[ 0 ])[ 1 ])[ 0 ];
    ((result[ 0 ])[ 1 ])[ 1 ] = ((tmp0[ 0 ])[ 1 ])[ 1 ];
    ((result[ 1 ])[ 0 ])[ 0 ] = ((tmp0[ 1 ])[ 0 ])[ 0 ];
    ((result[ 1 ])[ 0 ])[ 1 ] = ((tmp0[ 1 ])[ 0 ])[ 1 ];
    ((result[ 1 ])[ 1 ])[ 0 ] = ((tmp0[ 1 ])[ 1 ])[ 0 ];
    ((result[ 1 ])[ 1 ])[ 1 ] = ((tmp0[ 1 ])[ 1 ])[ 1 ];
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffdisplacement > > coefficients_;
};

} // namespace UFLLocalFunctions_cfa413bac5f9805f883691aeaf214600

PYBIND11_MODULE( localfunction_cfa413bac5f9805f883691aeaf214600_f084d41ce9ff164524e5c3ea6f8082cc, module )
{
  typedef UFLLocalFunctions_cfa413bac5f9805f883691aeaf214600::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::Grid<2,2>::LeafGridView>, Dune::FieldVector< double, 2 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_cfa413bac5f9805f883691aeaf214600::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::Grid<2,2>::LeafGridView>, Dune::FieldVector< double, 2 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_cfa413bac5f9805f883691aeaf214600_f084d41ce9ff164524e5c3ea6f8082cc.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::Grid<2,2>::LeafGridView>, Dune::FieldVector< double, 2 > > &coeffdisplacement ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::Vem::Grid<2,2>::LeafGridView>(gridView),name,order, coeffdisplacement); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
