#ifndef GUARD_0bb034af7b5f344bba30482c71ee2d67
#define GUARD_0bb034af7b5f344bba30482c71ee2d67

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

namespace UFLLocalFunctions_0bb034af7b5f344bba30482c71ee2d67
{

  // UFLLocalFunction
// ----------------

template< class GridPart >
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
    result[ 0 ] = tmp12;
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
    const auto tmp8 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp9 = 3 * tmp8;
    const auto tmp10 = -1 * tmp9;
    const auto tmp11 = 0.3333333333333333 * tmp10;
    const auto tmp12 = -1 * tmp11;
    const auto tmp13 = 2 * tmp0[ 1 ];
    const auto tmp14 = 0.5 * tmp13;
    const auto tmp15 = tmp14 + tmp12;
    (result[ 0 ])[ 0 ] = tmp7;
    (result[ 0 ])[ 1 ] = tmp15;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 2 * tmp0[ 0 ];
    const auto tmp2 = 3 * tmp1;
    const auto tmp3 = 0.3333333333333333 * tmp2;
    const auto tmp4 = -1 * tmp3;
    const auto tmp5 = 1.0 + tmp4;
    const auto tmp6 = 2 * tmp0[ 1 ];
    const auto tmp7 = 3 * tmp6;
    const auto tmp8 = -1 * tmp7;
    const auto tmp9 = 0.3333333333333333 * tmp8;
    const auto tmp10 = -1 * tmp9;
    const auto tmp11 = 1.0 + tmp10;
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp5;
    ((result[ 0 ])[ 0 ])[ 1 ] = 0.0;
    ((result[ 0 ])[ 1 ])[ 0 ] = 0.0;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp11;
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

} // namespace UFLLocalFunctions_0bb034af7b5f344bba30482c71ee2d67

PYBIND11_MODULE( localfunction_0bb034af7b5f344bba30482c71ee2d67_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_0bb034af7b5f344bba30482c71ee2d67::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_0bb034af7b5f344bba30482c71ee2d67::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_0bb034af7b5f344bba30482c71ee2d67_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
