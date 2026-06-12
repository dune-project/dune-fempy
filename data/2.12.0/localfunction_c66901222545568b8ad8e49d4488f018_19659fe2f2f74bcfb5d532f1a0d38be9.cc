#ifndef GUARD_c66901222545568b8ad8e49d4488f018
#define GUARD_c66901222545568b8ad8e49d4488f018

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

namespace UFLLocalFunctions_c66901222545568b8ad8e49d4488f018
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
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 6.283185307179586 * tmp0[ 1 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::cos( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = 6.283185307179586 * tmp0[ 0 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = std::cos( tmp3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp4 * tmp2;
      assert(tmp5 == tmp5);
    result[ 0 ] = tmp5;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 6.283185307179586 * tmp0[ 1 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::cos( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = 6.283185307179586 * tmp0[ 0 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = std::sin( tmp3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = 6.283185307179586 * tmp5;
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 * tmp2;
      assert(tmp7 == tmp7);
    const auto tmp8 = std::cos( tmp3 );
      assert(tmp8 == tmp8);
    const auto tmp9 = std::sin( tmp1 );
      assert(tmp9 == tmp9);
    const auto tmp10 = -1 * tmp9;
      assert(tmp10 == tmp10);
    const auto tmp11 = 6.283185307179586 * tmp10;
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 * tmp8;
      assert(tmp12 == tmp12);
    (result[ 0 ])[ 0 ] = tmp7;
    (result[ 0 ])[ 1 ] = tmp12;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 6.283185307179586 * tmp0[ 1 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::cos( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = 6.283185307179586 * tmp0[ 0 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = std::cos( tmp3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = 6.283185307179586 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = -1 * tmp5;
      assert(tmp6 == tmp6);
    const auto tmp7 = 6.283185307179586 * tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 * tmp2;
      assert(tmp8 == tmp8);
    const auto tmp9 = std::sin( tmp1 );
      assert(tmp9 == tmp9);
    const auto tmp10 = -1 * tmp9;
      assert(tmp10 == tmp10);
    const auto tmp11 = 6.283185307179586 * tmp10;
      assert(tmp11 == tmp11);
    const auto tmp12 = std::sin( tmp3 );
      assert(tmp12 == tmp12);
    const auto tmp13 = -1 * tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = 6.283185307179586 * tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 * tmp11;
      assert(tmp15 == tmp15);
    const auto tmp16 = 6.283185307179586 * tmp2;
      assert(tmp16 == tmp16);
    const auto tmp17 = -1 * tmp16;
      assert(tmp17 == tmp17);
    const auto tmp18 = 6.283185307179586 * tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 * tmp4;
      assert(tmp19 == tmp19);
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp8;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp15;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp15;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp19;
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

} // namespace UFLLocalFunctions_c66901222545568b8ad8e49d4488f018

PYBIND11_MODULE( localfunction_c66901222545568b8ad8e49d4488f018_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_c66901222545568b8ad8e49d4488f018::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_c66901222545568b8ad8e49d4488f018::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_c66901222545568b8ad8e49d4488f018_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
