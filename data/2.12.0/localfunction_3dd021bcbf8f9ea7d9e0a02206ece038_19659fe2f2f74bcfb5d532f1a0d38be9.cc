#ifndef GUARD_3dd021bcbf8f9ea7d9e0a02206ece038
#define GUARD_3dd021bcbf8f9ea7d9e0a02206ece038

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

namespace UFLLocalFunctions_3dd021bcbf8f9ea7d9e0a02206ece038
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
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 + tmp0[ 0 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::pow( tmp1, 2 );
      assert(tmp2 == tmp2);
    const auto tmp3 = std::pow( tmp0[ 0 ], 2 );
      assert(tmp3 == tmp3);
    const auto tmp4 = 256 * tmp0[ 1 ];
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 + tmp0[ 1 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp5 * tmp4;
      assert(tmp6 == tmp6);
    const auto tmp7 = 2 * tmp0[ 1 ];
      assert(tmp7 == tmp7);
    const auto tmp8 = -1 + tmp7;
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp8 * tmp6;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 * tmp3;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * tmp2;
      assert(tmp11 == tmp11);
    result[ 0 ] = tmp11;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 + tmp0[ 0 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::pow( tmp1, 2 );
      assert(tmp2 == tmp2);
    const auto tmp3 = 256 * tmp0[ 1 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = -1 + tmp0[ 1 ];
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp4 * tmp3;
      assert(tmp5 == tmp5);
    const auto tmp6 = 2 * tmp0[ 1 ];
      assert(tmp6 == tmp6);
    const auto tmp7 = -1 + tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 * tmp5;
      assert(tmp8 == tmp8);
    const auto tmp9 = 2 * tmp0[ 0 ];
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 * tmp8;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * tmp2;
      assert(tmp11 == tmp11);
    const auto tmp12 = std::pow( tmp0[ 0 ], 2 );
      assert(tmp12 == tmp12);
    const auto tmp13 = tmp8 * tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = 2 * tmp1;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 * tmp13;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 + tmp11;
      assert(tmp16 == tmp16);
    const auto tmp17 = 2 * tmp5;
      assert(tmp17 == tmp17);
    const auto tmp18 = 256 * tmp4;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp3 + tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp7 * tmp19;
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp20 + tmp17;
      assert(tmp21 == tmp21);
    const auto tmp22 = tmp21 * tmp12;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 * tmp2;
      assert(tmp23 == tmp23);
    (result[ 0 ])[ 0 ] = tmp16;
    (result[ 0 ])[ 1 ] = tmp23;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 + tmp0[ 0 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::pow( tmp1, 2 );
      assert(tmp2 == tmp2);
    const auto tmp3 = 256 * tmp0[ 1 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = -1 + tmp0[ 1 ];
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp4 * tmp3;
      assert(tmp5 == tmp5);
    const auto tmp6 = 2 * tmp0[ 1 ];
      assert(tmp6 == tmp6);
    const auto tmp7 = -1 + tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 * tmp5;
      assert(tmp8 == tmp8);
    const auto tmp9 = 2 * tmp8;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 * tmp2;
      assert(tmp10 == tmp10);
    const auto tmp11 = 2 * tmp0[ 0 ];
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 * tmp8;
      assert(tmp12 == tmp12);
    const auto tmp13 = 2 * tmp1;
      assert(tmp13 == tmp13);
    const auto tmp14 = tmp13 * tmp12;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 + tmp10;
      assert(tmp15 == tmp15);
    const auto tmp16 = std::pow( tmp0[ 0 ], 2 );
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp8 * tmp16;
      assert(tmp17 == tmp17);
    const auto tmp18 = 2 * tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp14 + tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp19 + tmp15;
      assert(tmp20 == tmp20);
    const auto tmp21 = 2 * tmp5;
      assert(tmp21 == tmp21);
    const auto tmp22 = 256 * tmp4;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp3 + tmp22;
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp7 * tmp23;
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp24 + tmp21;
      assert(tmp25 == tmp25);
    const auto tmp26 = tmp25 * tmp11;
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp26 * tmp2;
      assert(tmp27 == tmp27);
    const auto tmp28 = tmp25 * tmp16;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp13 * tmp28;
      assert(tmp29 == tmp29);
    const auto tmp30 = tmp29 + tmp27;
      assert(tmp30 == tmp30);
    const auto tmp31 = 2 * tmp23;
      assert(tmp31 == tmp31);
    const auto tmp32 = 512 * tmp7;
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp32 + tmp31;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 + tmp31;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp34 * tmp16;
      assert(tmp35 == tmp35);
    const auto tmp36 = tmp35 * tmp2;
      assert(tmp36 == tmp36);
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp20;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp30;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp30;
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

} // namespace UFLLocalFunctions_3dd021bcbf8f9ea7d9e0a02206ece038

PYBIND11_MODULE( localfunction_3dd021bcbf8f9ea7d9e0a02206ece038_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_3dd021bcbf8f9ea7d9e0a02206ece038::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_3dd021bcbf8f9ea7d9e0a02206ece038::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_3dd021bcbf8f9ea7d9e0a02206ece038_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
