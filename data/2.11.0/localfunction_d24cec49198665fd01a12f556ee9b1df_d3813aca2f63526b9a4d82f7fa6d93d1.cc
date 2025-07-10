#ifndef GUARD_d24cec49198665fd01a12f556ee9b1df
#define GUARD_d24cec49198665fd01a12f556ee9b1df

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

namespace UFLLocalFunctions_d24cec49198665fd01a12f556ee9b1df
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
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 0.7853981633974483 * tmp0[ 2 ];
    const auto tmp2 = std::cos( tmp1 );
    const auto tmp3 = tmp0[ 0 ] + tmp0[ 1 ];
    const auto tmp4 = 6.283185307179586 * tmp3;
    const auto tmp5 = std::sin( tmp4 );
    const auto tmp6 = 0.5 * tmp5;
    const auto tmp7 = tmp6 * tmp2;
    const auto tmp8 = 1 + tmp7;
    const auto tmp9 = tmp0[ 0 ] * tmp8;
    const auto tmp10 = tmp0[ 1 ] * tmp8;
    const auto tmp11 = tmp0[ 2 ] * tmp8;
    result[ 0 ] = tmp9;
    result[ 1 ] = tmp10;
    result[ 2 ] = tmp11;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 0.7853981633974483 * tmp0[ 2 ];
    const auto tmp2 = std::cos( tmp1 );
    const auto tmp3 = tmp0[ 0 ] + tmp0[ 1 ];
    const auto tmp4 = 6.283185307179586 * tmp3;
    const auto tmp5 = std::cos( tmp4 );
    const auto tmp6 = 6.283185307179586 * tmp5;
    const auto tmp7 = 0.5 * tmp6;
    const auto tmp8 = tmp7 * tmp2;
    const auto tmp9 = tmp0[ 0 ] * tmp8;
    const auto tmp10 = std::sin( tmp4 );
    const auto tmp11 = 0.5 * tmp10;
    const auto tmp12 = tmp11 * tmp2;
    const auto tmp13 = 1 + tmp12;
    const auto tmp14 = tmp13 + tmp9;
    const auto tmp15 = std::sin( tmp1 );
    const auto tmp16 = -1 * tmp15;
    const auto tmp17 = 0.7853981633974483 * tmp16;
    const auto tmp18 = tmp17 * tmp11;
    const auto tmp19 = tmp0[ 0 ] * tmp18;
    const auto tmp20 = tmp0[ 1 ] * tmp8;
    const auto tmp21 = tmp13 + tmp20;
    const auto tmp22 = tmp0[ 1 ] * tmp18;
    const auto tmp23 = tmp0[ 2 ] * tmp8;
    const auto tmp24 = tmp0[ 2 ] * tmp18;
    const auto tmp25 = tmp13 + tmp24;
    (result[ 0 ])[ 0 ] = tmp14;
    (result[ 0 ])[ 1 ] = tmp9;
    (result[ 0 ])[ 2 ] = tmp19;
    (result[ 1 ])[ 0 ] = tmp20;
    (result[ 1 ])[ 1 ] = tmp21;
    (result[ 1 ])[ 2 ] = tmp22;
    (result[ 2 ])[ 0 ] = tmp23;
    (result[ 2 ])[ 1 ] = tmp23;
    (result[ 2 ])[ 2 ] = tmp25;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 0.7853981633974483 * tmp0[ 2 ];
    const auto tmp2 = std::cos( tmp1 );
    const auto tmp3 = tmp0[ 0 ] + tmp0[ 1 ];
    const auto tmp4 = 6.283185307179586 * tmp3;
    const auto tmp5 = std::cos( tmp4 );
    const auto tmp6 = 6.283185307179586 * tmp5;
    const auto tmp7 = 0.5 * tmp6;
    const auto tmp8 = tmp7 * tmp2;
    const auto tmp9 = std::sin( tmp4 );
    const auto tmp10 = -1 * tmp9;
    const auto tmp11 = 6.283185307179586 * tmp10;
    const auto tmp12 = 6.283185307179586 * tmp11;
    const auto tmp13 = 0.5 * tmp12;
    const auto tmp14 = tmp13 * tmp2;
    const auto tmp15 = tmp0[ 0 ] * tmp14;
    const auto tmp16 = tmp15 + tmp8;
    const auto tmp17 = tmp16 + tmp8;
    const auto tmp18 = 0.5 * tmp9;
    const auto tmp19 = std::sin( tmp1 );
    const auto tmp20 = -1 * tmp19;
    const auto tmp21 = 0.7853981633974483 * tmp20;
    const auto tmp22 = tmp21 * tmp18;
    const auto tmp23 = tmp7 * tmp21;
    const auto tmp24 = tmp0[ 0 ] * tmp23;
    const auto tmp25 = tmp24 + tmp22;
    const auto tmp26 = 0.7853981633974483 * tmp2;
    const auto tmp27 = -1 * tmp26;
    const auto tmp28 = 0.7853981633974483 * tmp27;
    const auto tmp29 = tmp28 * tmp18;
    const auto tmp30 = tmp0[ 0 ] * tmp29;
    const auto tmp31 = tmp0[ 1 ] * tmp14;
    const auto tmp32 = tmp31 + tmp8;
    const auto tmp33 = tmp0[ 1 ] * tmp23;
    const auto tmp34 = tmp32 + tmp8;
    const auto tmp35 = tmp33 + tmp22;
    const auto tmp36 = tmp0[ 1 ] * tmp29;
    const auto tmp37 = tmp0[ 2 ] * tmp14;
    const auto tmp38 = tmp0[ 2 ] * tmp23;
    const auto tmp39 = tmp38 + tmp8;
    const auto tmp40 = tmp0[ 2 ] * tmp29;
    const auto tmp41 = tmp40 + tmp22;
    const auto tmp42 = tmp41 + tmp22;
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp17;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp16;
    ((result[ 0 ])[ 0 ])[ 2 ] = tmp25;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp16;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp15;
    ((result[ 0 ])[ 1 ])[ 2 ] = tmp24;
    ((result[ 0 ])[ 2 ])[ 0 ] = tmp25;
    ((result[ 0 ])[ 2 ])[ 1 ] = tmp24;
    ((result[ 0 ])[ 2 ])[ 2 ] = tmp30;
    ((result[ 1 ])[ 0 ])[ 0 ] = tmp31;
    ((result[ 1 ])[ 0 ])[ 1 ] = tmp32;
    ((result[ 1 ])[ 0 ])[ 2 ] = tmp33;
    ((result[ 1 ])[ 1 ])[ 0 ] = tmp32;
    ((result[ 1 ])[ 1 ])[ 1 ] = tmp34;
    ((result[ 1 ])[ 1 ])[ 2 ] = tmp35;
    ((result[ 1 ])[ 2 ])[ 0 ] = tmp33;
    ((result[ 1 ])[ 2 ])[ 1 ] = tmp35;
    ((result[ 1 ])[ 2 ])[ 2 ] = tmp36;
    ((result[ 2 ])[ 0 ])[ 0 ] = tmp37;
    ((result[ 2 ])[ 0 ])[ 1 ] = tmp37;
    ((result[ 2 ])[ 0 ])[ 2 ] = tmp39;
    ((result[ 2 ])[ 1 ])[ 0 ] = tmp37;
    ((result[ 2 ])[ 1 ])[ 1 ] = tmp37;
    ((result[ 2 ])[ 1 ])[ 2 ] = tmp39;
    ((result[ 2 ])[ 2 ])[ 0 ] = tmp39;
    ((result[ 2 ])[ 2 ])[ 1 ] = tmp39;
    ((result[ 2 ])[ 2 ])[ 2 ] = tmp42;
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

} // namespace UFLLocalFunctions_d24cec49198665fd01a12f556ee9b1df

PYBIND11_MODULE( localfunction_d24cec49198665fd01a12f556ee9b1df_d3813aca2f63526b9a4d82f7fa6d93d1, module )
{
  typedef UFLLocalFunctions_d24cec49198665fd01a12f556ee9b1df::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_d24cec49198665fd01a12f556ee9b1df::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_d24cec49198665fd01a12f556ee9b1df_d3813aca2f63526b9a4d82f7fa6d93d1.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
