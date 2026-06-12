#ifndef GUARD_b10ff5627bd2f22e898386281281f7db
#define GUARD_b10ff5627bd2f22e898386281281f7db

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

namespace UFLLocalFunctions_b10ff5627bd2f22e898386281281f7db
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
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = std::cos( tmp0[ 0 ] );
      assert(tmp1 == tmp1);
    const auto tmp2 = -1 * tmp0[ 0 ];
      assert(tmp2 == tmp2);
    const auto tmp3 = 25 + tmp2;
      assert(tmp3 == tmp3);
    const auto tmp4 = tmp3 * tmp1;
      assert(tmp4 == tmp4);
    const auto tmp5 = std::sin( tmp0[ 0 ] );
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp3 * tmp5;
      assert(tmp6 == tmp6);
    result[ 0 ] = tmp4;
    result[ 1 ] = tmp6;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = std::cos( tmp0[ 0 ] );
      assert(tmp1 == tmp1);
    const auto tmp2 = -1 * tmp1;
      assert(tmp2 == tmp2);
    const auto tmp3 = std::sin( tmp0[ 0 ] );
      assert(tmp3 == tmp3);
    const auto tmp4 = -1 * tmp3;
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * tmp0[ 0 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = 25 + tmp5;
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 * tmp4;
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 + tmp2;
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp6 * tmp1;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 + tmp4;
      assert(tmp10 == tmp10);
    (result[ 0 ])[ 0 ] = tmp8;
    (result[ 1 ])[ 0 ] = tmp10;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = std::sin( tmp0[ 0 ] );
      assert(tmp1 == tmp1);
    const auto tmp2 = -1 * tmp1;
      assert(tmp2 == tmp2);
    const auto tmp3 = -1 * tmp2;
      assert(tmp3 == tmp3);
    const auto tmp4 = std::cos( tmp0[ 0 ] );
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = -1 * tmp0[ 0 ];
      assert(tmp6 == tmp6);
    const auto tmp7 = 25 + tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 * tmp5;
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp8 + tmp3;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 + tmp3;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp7 * tmp2;
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 + tmp5;
      assert(tmp12 == tmp12);
    const auto tmp13 = tmp12 + tmp5;
      assert(tmp13 == tmp13);
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp10;
    ((result[ 1 ])[ 0 ])[ 0 ] = tmp13;
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

} // namespace UFLLocalFunctions_b10ff5627bd2f22e898386281281f7db

PYBIND11_MODULE( localfunction_b10ff5627bd2f22e898386281281f7db_e14d829faa8914e3ff98a63559ccacdd, module )
{
  typedef UFLLocalFunctions_b10ff5627bd2f22e898386281281f7db::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_b10ff5627bd2f22e898386281281f7db::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_b10ff5627bd2f22e898386281281f7db_e14d829faa8914e3ff98a63559ccacdd.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
