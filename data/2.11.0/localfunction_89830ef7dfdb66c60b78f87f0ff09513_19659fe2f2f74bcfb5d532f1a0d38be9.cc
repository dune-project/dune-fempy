#ifndef GUARD_89830ef7dfdb66c60b78f87f0ff09513
#define GUARD_89830ef7dfdb66c60b78f87f0ff09513

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

namespace UFLLocalFunctions_89830ef7dfdb66c60b78f87f0ff09513
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
    const auto tmp1 = -1 + tmp0[ 1 ];
    const auto tmp2 = std::pow( tmp1, 2 );
    const auto tmp3 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp4 = 256 * tmp0[ 0 ];
    const auto tmp5 = -1 + tmp0[ 0 ];
    const auto tmp6 = tmp5 * tmp4;
    const auto tmp7 = 2 * tmp0[ 0 ];
    const auto tmp8 = -1 + tmp7;
    const auto tmp9 = tmp8 * tmp6;
    const auto tmp10 = tmp9 * tmp3;
    const auto tmp11 = tmp10 * tmp2;
    result[ 0 ] = tmp11;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = -1 + tmp0[ 1 ];
    const auto tmp2 = std::pow( tmp1, 2 );
    const auto tmp3 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp4 = 256 * tmp0[ 0 ];
    const auto tmp5 = -1 + tmp0[ 0 ];
    const auto tmp6 = tmp5 * tmp4;
    const auto tmp7 = 2 * tmp6;
    const auto tmp8 = 256 * tmp5;
    const auto tmp9 = tmp4 + tmp8;
    const auto tmp10 = 2 * tmp0[ 0 ];
    const auto tmp11 = -1 + tmp10;
    const auto tmp12 = tmp11 * tmp9;
    const auto tmp13 = tmp12 + tmp7;
    const auto tmp14 = tmp13 * tmp3;
    const auto tmp15 = tmp14 * tmp2;
    const auto tmp16 = tmp11 * tmp6;
    const auto tmp17 = 2 * tmp0[ 1 ];
    const auto tmp18 = tmp17 * tmp16;
    const auto tmp19 = tmp18 * tmp2;
    const auto tmp20 = tmp16 * tmp3;
    const auto tmp21 = 2 * tmp1;
    const auto tmp22 = tmp21 * tmp20;
    const auto tmp23 = tmp22 + tmp19;
    (result[ 0 ])[ 0 ] = tmp15;
    (result[ 0 ])[ 1 ] = tmp23;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = -1 + tmp0[ 1 ];
    const auto tmp2 = std::pow( tmp1, 2 );
    const auto tmp3 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp4 = -1 + tmp0[ 0 ];
    const auto tmp5 = 256 * tmp4;
    const auto tmp6 = 256 * tmp0[ 0 ];
    const auto tmp7 = tmp6 + tmp5;
    const auto tmp8 = 2 * tmp7;
    const auto tmp9 = 2 * tmp0[ 0 ];
    const auto tmp10 = -1 + tmp9;
    const auto tmp11 = 512 * tmp10;
    const auto tmp12 = tmp11 + tmp8;
    const auto tmp13 = tmp12 + tmp8;
    const auto tmp14 = tmp13 * tmp3;
    const auto tmp15 = tmp14 * tmp2;
    const auto tmp16 = 2 * tmp0[ 1 ];
    const auto tmp17 = tmp4 * tmp6;
    const auto tmp18 = 2 * tmp17;
    const auto tmp19 = tmp10 * tmp7;
    const auto tmp20 = tmp19 + tmp18;
    const auto tmp21 = tmp20 * tmp16;
    const auto tmp22 = tmp21 * tmp2;
    const auto tmp23 = tmp20 * tmp3;
    const auto tmp24 = 2 * tmp1;
    const auto tmp25 = tmp24 * tmp23;
    const auto tmp26 = tmp25 + tmp22;
    const auto tmp27 = tmp10 * tmp17;
    const auto tmp28 = 2 * tmp27;
    const auto tmp29 = tmp28 * tmp2;
    const auto tmp30 = tmp16 * tmp27;
    const auto tmp31 = tmp24 * tmp30;
    const auto tmp32 = tmp31 + tmp29;
    const auto tmp33 = tmp27 * tmp3;
    const auto tmp34 = 2 * tmp33;
    const auto tmp35 = tmp31 + tmp34;
    const auto tmp36 = tmp35 + tmp32;
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp15;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp26;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp26;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp36;
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

} // namespace UFLLocalFunctions_89830ef7dfdb66c60b78f87f0ff09513

PYBIND11_MODULE( localfunction_89830ef7dfdb66c60b78f87f0ff09513_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_89830ef7dfdb66c60b78f87f0ff09513::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_89830ef7dfdb66c60b78f87f0ff09513::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_89830ef7dfdb66c60b78f87f0ff09513_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
