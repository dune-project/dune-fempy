#ifndef GUARD_cb799f09cf3672c39d40b01f45e766d4
#define GUARD_cb799f09cf3672c39d40b01f45e766d4

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/common/exceptions.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_cb799f09cf3672c39d40b01f45e766d4
{

  // UFLLocalFunction
// ----------------

template< class GridPart >
struct UFLLocalFunction
  : public Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<3>>
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridViewType GridView;
  typedef typename GridView::ctype ctype;
  typedef Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<3>> BaseType;
  typedef Dune::Fem::GridFunctionSpace<GridPartType,Dune::Dim<3>> FunctionSpaceType;
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
    result[ 0 ] = tmp0[ 0 ];
    result[ 1 ] = tmp0[ 1 ];
    result[ 2 ] = tmp0[ 2 ];
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    (result[ 0 ])[ 0 ] = 1;
    (result[ 0 ])[ 1 ] = 0.0;
    (result[ 0 ])[ 2 ] = 0.0;
    (result[ 1 ])[ 0 ] = 0.0;
    (result[ 1 ])[ 1 ] = 1;
    (result[ 1 ])[ 2 ] = 0.0;
    (result[ 2 ])[ 0 ] = 0.0;
    (result[ 2 ])[ 1 ] = 0.0;
    (result[ 2 ])[ 2 ] = 1;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    ((result[ 0 ])[ 0 ])[ 0 ] = 0.0;
    ((result[ 0 ])[ 0 ])[ 1 ] = 0.0;
    ((result[ 0 ])[ 0 ])[ 2 ] = 0.0;
    ((result[ 0 ])[ 1 ])[ 0 ] = 0.0;
    ((result[ 0 ])[ 1 ])[ 1 ] = 0.0;
    ((result[ 0 ])[ 1 ])[ 2 ] = 0.0;
    ((result[ 0 ])[ 2 ])[ 0 ] = 0.0;
    ((result[ 0 ])[ 2 ])[ 1 ] = 0.0;
    ((result[ 0 ])[ 2 ])[ 2 ] = 0.0;
    ((result[ 1 ])[ 0 ])[ 0 ] = 0.0;
    ((result[ 1 ])[ 0 ])[ 1 ] = 0.0;
    ((result[ 1 ])[ 0 ])[ 2 ] = 0.0;
    ((result[ 1 ])[ 1 ])[ 0 ] = 0.0;
    ((result[ 1 ])[ 1 ])[ 1 ] = 0.0;
    ((result[ 1 ])[ 1 ])[ 2 ] = 0.0;
    ((result[ 1 ])[ 2 ])[ 0 ] = 0.0;
    ((result[ 1 ])[ 2 ])[ 1 ] = 0.0;
    ((result[ 1 ])[ 2 ])[ 2 ] = 0.0;
    ((result[ 2 ])[ 0 ])[ 0 ] = 0.0;
    ((result[ 2 ])[ 0 ])[ 1 ] = 0.0;
    ((result[ 2 ])[ 0 ])[ 2 ] = 0.0;
    ((result[ 2 ])[ 1 ])[ 0 ] = 0.0;
    ((result[ 2 ])[ 1 ])[ 1 ] = 0.0;
    ((result[ 2 ])[ 1 ])[ 2 ] = 0.0;
    ((result[ 2 ])[ 2 ])[ 0 ] = 0.0;
    ((result[ 2 ])[ 2 ])[ 1 ] = 0.0;
    ((result[ 2 ])[ 2 ])[ 2 ] = 0.0;
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

} // namespace UFLLocalFunctions_cb799f09cf3672c39d40b01f45e766d4

PYBIND11_MODULE( localfunction_cb799f09cf3672c39d40b01f45e766d4_d3813aca2f63526b9a4d82f7fa6d93d1, module )
{
  typedef UFLLocalFunctions_cb799f09cf3672c39d40b01f45e766d4::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_cb799f09cf3672c39d40b01f45e766d4::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_cb799f09cf3672c39d40b01f45e766d4_d3813aca2f63526b9a4d82f7fa6d93d1.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
