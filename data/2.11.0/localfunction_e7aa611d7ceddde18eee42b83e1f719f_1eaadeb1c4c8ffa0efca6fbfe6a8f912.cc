#ifndef GUARD_e7aa611d7ceddde18eee42b83e1f719f
#define GUARD_e7aa611d7ceddde18eee42b83e1f719f

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/polygongrid/dgf.hh>
#include <dune/polygongrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/common/exceptions.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_e7aa611d7ceddde18eee42b83e1f719f
{

  // UFLLocalFunction
// ----------------

template< class GridPart >
struct UFLLocalFunction
  : public Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<4>>
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridViewType GridView;
  typedef typename GridView::ctype ctype;
  typedef Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<4>> BaseType;
  typedef Dune::Fem::GridFunctionSpace<GridPartType,Dune::Dim<4>> FunctionSpaceType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef typename EntityType::Geometry Geometry;
  typedef typename Geometry::GlobalCoordinate GlobalCoordinateType;
  typedef Dune::Fem::IntersectionSide Side;
  typedef std::tuple<> ConstantTupleType;
  typedef std::tuple<> CoefficientTupleType;
  static constexpr bool gridPartValid = true;
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order)
  {}

  void bind ( const IntersectionType &intersection, Side side )
  {
    BaseType::bind(intersection, side);
  }

  void bind ( const EntityType &entity )
  {
    BaseType::bind(entity);
  }

  void unbind ()
  {
    BaseType::unbind();
  }

  template< class Point >
  void evaluate ( const Point &x, typename FunctionSpaceType::RangeType &result ) const
  {
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
    const auto tmp2 = tmp0[ 0 ] * tmp0[ 0 ];
    const auto tmp3 = tmp2 + tmp1;
    result[ 0 ] = (tmp3 < 0.1 ? 1 : 0.125);
    result[ 1 ] = (tmp3 < 0.1 ? 0.0 : 0.0);
    result[ 2 ] = (tmp3 < 0.1 ? 0.0 : 0.0);
    result[ 3 ] = (tmp3 < 0.1 ? 2.5 : 0.25);
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    (result[ 0 ])[ 0 ] = 0.0;
    (result[ 0 ])[ 1 ] = 0.0;
    (result[ 1 ])[ 0 ] = 0.0;
    (result[ 1 ])[ 1 ] = 0.0;
    (result[ 2 ])[ 0 ] = 0.0;
    (result[ 2 ])[ 1 ] = 0.0;
    (result[ 3 ])[ 0 ] = 0.0;
    (result[ 3 ])[ 1 ] = 0.0;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    ((result[ 0 ])[ 0 ])[ 0 ] = 0.0;
    ((result[ 0 ])[ 0 ])[ 1 ] = 0.0;
    ((result[ 0 ])[ 1 ])[ 0 ] = 0.0;
    ((result[ 0 ])[ 1 ])[ 1 ] = 0.0;
    ((result[ 1 ])[ 0 ])[ 0 ] = 0.0;
    ((result[ 1 ])[ 0 ])[ 1 ] = 0.0;
    ((result[ 1 ])[ 1 ])[ 0 ] = 0.0;
    ((result[ 1 ])[ 1 ])[ 1 ] = 0.0;
    ((result[ 2 ])[ 0 ])[ 0 ] = 0.0;
    ((result[ 2 ])[ 0 ])[ 1 ] = 0.0;
    ((result[ 2 ])[ 1 ])[ 0 ] = 0.0;
    ((result[ 2 ])[ 1 ])[ 1 ] = 0.0;
    ((result[ 3 ])[ 0 ])[ 0 ] = 0.0;
    ((result[ 3 ])[ 0 ])[ 1 ] = 0.0;
    ((result[ 3 ])[ 1 ])[ 0 ] = 0.0;
    ((result[ 3 ])[ 1 ])[ 1 ] = 0.0;
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
  ConstantTupleType constants_;
  std::tuple<  > coefficients_;
};

} // namespace UFLLocalFunctions_e7aa611d7ceddde18eee42b83e1f719f

PYBIND11_MODULE( localfunction_e7aa611d7ceddde18eee42b83e1f719f_1eaadeb1c4c8ffa0efca6fbfe6a8f912, module )
{
  typedef UFLLocalFunctions_e7aa611d7ceddde18eee42b83e1f719f::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::PolygonGrid< double >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_e7aa611d7ceddde18eee42b83e1f719f::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::PolygonGrid< double >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_e7aa611d7ceddde18eee42b83e1f719f_1eaadeb1c4c8ffa0efca6fbfe6a8f912.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::PolygonGrid< double >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
