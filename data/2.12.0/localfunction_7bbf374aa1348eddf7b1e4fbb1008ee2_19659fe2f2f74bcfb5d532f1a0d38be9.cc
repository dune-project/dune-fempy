#ifndef GUARD_7bbf374aa1348eddf7b1e4fbb1008ee2
#define GUARD_7bbf374aa1348eddf7b1e4fbb1008ee2

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
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_7bbf374aa1348eddf7b1e4fbb1008ee2
{

  // UFLLocalFunction
// ----------------

template< class GridPart >
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

  void decouple ()
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
      assert(tmp0 == tmp0);
    const auto tmp1 = tmp0[ 0 ] + tmp0[ 1 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = tmp1 / 1.4142135623730951;
      assert(tmp2 == tmp2);
    const auto tmp3 = -1 * tmp0[ 0 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = tmp0[ 1 ] + tmp3;
      assert(tmp4 == tmp4);
    const auto tmp5 = 1.4142135623730951 * tmp4;
      assert(tmp5 == tmp5);
    result[ 0 ] = tmp2;
    result[ 1 ] = tmp5;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    (result[ 0 ])[ 0 ] = 0.7071067811865475;
    (result[ 0 ])[ 1 ] = 0.7071067811865475;
    (result[ 1 ])[ 0 ] = -1.4142135623730951;
    (result[ 1 ])[ 1 ] = 1.4142135623730951;
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

} // namespace UFLLocalFunctions_7bbf374aa1348eddf7b1e4fbb1008ee2

PYBIND11_MODULE( localfunction_7bbf374aa1348eddf7b1e4fbb1008ee2_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_7bbf374aa1348eddf7b1e4fbb1008ee2::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_7bbf374aa1348eddf7b1e4fbb1008ee2::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_7bbf374aa1348eddf7b1e4fbb1008ee2_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
