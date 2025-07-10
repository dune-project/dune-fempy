#ifndef GUARD_5cc181267fddc3c844daa3bec40763e7
#define GUARD_5cc181267fddc3c844daa3bec40763e7

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
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_5cc181267fddc3c844daa3bec40763e7
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
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = -0.5 + tmp0[ 1 ];
    const auto tmp2 = 12.566370614359172 * tmp1;
    const auto tmp3 = std::sin( tmp2 );
    const auto tmp4 = std::pow( tmp3, 3 );
    const auto tmp5 = -0.5 + tmp0[ 0 ];
    const auto tmp6 = 12.566370614359172 * tmp5;
    const auto tmp7 = std::sin( tmp6 );
    const auto tmp8 = std::pow( tmp7, 3 );
    const auto tmp9 = (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : tmp8)) * (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : tmp4));
    const auto tmp10 = 0.5 * tmp9;
    result[ 0 ] = tmp10;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = -0.5 + tmp0[ 1 ];
    const auto tmp2 = 12.566370614359172 * tmp1;
    const auto tmp3 = std::sin( tmp2 );
    const auto tmp4 = std::pow( tmp3, 3 );
    const auto tmp5 = -0.5 + tmp0[ 0 ];
    const auto tmp6 = 12.566370614359172 * tmp5;
    const auto tmp7 = std::sin( tmp6 );
    const auto tmp8 = std::pow( tmp7, 2 );
    const auto tmp9 = std::cos( tmp6 );
    const auto tmp10 = 12.566370614359172 * tmp9;
    const auto tmp11 = 3 * tmp10;
    const auto tmp12 = tmp11 * tmp8;
    const auto tmp13 = (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : tmp12)) * (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : tmp4));
    const auto tmp14 = std::pow( tmp7, 3 );
    const auto tmp15 = (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : 0.0)) * (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : tmp14));
    const auto tmp16 = tmp15 + tmp13;
    const auto tmp17 = 0.5 * tmp16;
    const auto tmp18 = (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : 0.0)) * (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : tmp4));
    const auto tmp19 = std::pow( tmp3, 2 );
    const auto tmp20 = std::cos( tmp2 );
    const auto tmp21 = 12.566370614359172 * tmp20;
    const auto tmp22 = 3 * tmp21;
    const auto tmp23 = tmp22 * tmp19;
    const auto tmp24 = (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : tmp23)) * (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : tmp14));
    const auto tmp25 = tmp24 + tmp18;
    const auto tmp26 = 0.5 * tmp25;
    (result[ 0 ])[ 0 ] = tmp17;
    (result[ 0 ])[ 1 ] = tmp26;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = -0.5 + tmp0[ 1 ];
    const auto tmp2 = 12.566370614359172 * tmp1;
    const auto tmp3 = std::sin( tmp2 );
    const auto tmp4 = std::pow( tmp3, 3 );
    const auto tmp5 = -0.5 + tmp0[ 0 ];
    const auto tmp6 = 12.566370614359172 * tmp5;
    const auto tmp7 = std::sin( tmp6 );
    const auto tmp8 = std::pow( tmp7, 2 );
    const auto tmp9 = -1 * tmp7;
    const auto tmp10 = 12.566370614359172 * tmp9;
    const auto tmp11 = 12.566370614359172 * tmp10;
    const auto tmp12 = 3 * tmp11;
    const auto tmp13 = tmp12 * tmp8;
    const auto tmp14 = std::cos( tmp6 );
    const auto tmp15 = 12.566370614359172 * tmp14;
    const auto tmp16 = 2 * tmp15;
    const auto tmp17 = tmp16 * tmp7;
    const auto tmp18 = 3 * tmp15;
    const auto tmp19 = tmp18 * tmp17;
    const auto tmp20 = tmp19 + tmp13;
    const auto tmp21 = (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : tmp20)) * (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : tmp4));
    const auto tmp22 = tmp18 * tmp8;
    const auto tmp23 = (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : 0.0)) * (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : tmp22));
    const auto tmp24 = tmp23 + tmp21;
    const auto tmp25 = std::pow( tmp7, 3 );
    const auto tmp26 = (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : 0.0)) * (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : tmp25));
    const auto tmp27 = tmp23 + tmp26;
    const auto tmp28 = tmp27 + tmp24;
    const auto tmp29 = 0.5 * tmp28;
    const auto tmp30 = (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : 0.0)) * (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : tmp4));
    const auto tmp31 = std::pow( tmp3, 2 );
    const auto tmp32 = std::cos( tmp2 );
    const auto tmp33 = 12.566370614359172 * tmp32;
    const auto tmp34 = 3 * tmp33;
    const auto tmp35 = tmp34 * tmp31;
    const auto tmp36 = (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : tmp22)) * (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : tmp35));
    const auto tmp37 = tmp36 + tmp30;
    const auto tmp38 = (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : 0.0)) * (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : 0.0));
    const auto tmp39 = tmp38 + tmp26;
    const auto tmp40 = tmp39 + tmp37;
    const auto tmp41 = 0.5 * tmp40;
    const auto tmp42 = tmp38 + tmp30;
    const auto tmp43 = tmp36 + tmp26;
    const auto tmp44 = tmp43 + tmp42;
    const auto tmp45 = 0.5 * tmp44;
    const auto tmp46 = (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : 0.0)) * (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : tmp35));
    const auto tmp47 = tmp46 + tmp30;
    const auto tmp48 = -1 * tmp3;
    const auto tmp49 = 12.566370614359172 * tmp48;
    const auto tmp50 = 12.566370614359172 * tmp49;
    const auto tmp51 = 3 * tmp50;
    const auto tmp52 = tmp51 * tmp31;
    const auto tmp53 = 2 * tmp33;
    const auto tmp54 = tmp53 * tmp3;
    const auto tmp55 = tmp34 * tmp54;
    const auto tmp56 = tmp55 + tmp52;
    const auto tmp57 = (tmp1 < -0.25 ? 0.0 : (tmp1 > 0.25 ? 0.0 : tmp56)) * (tmp5 < -0.25 ? 0.0 : (tmp5 > 0.25 ? 0.0 : tmp25));
    const auto tmp58 = tmp46 + tmp57;
    const auto tmp59 = tmp58 + tmp47;
    const auto tmp60 = 0.5 * tmp59;
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp29;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp41;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp45;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp60;
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

} // namespace UFLLocalFunctions_5cc181267fddc3c844daa3bec40763e7

PYBIND11_MODULE( localfunction_5cc181267fddc3c844daa3bec40763e7_d5631e01dd7008453673985fc163d5bb, module )
{
  typedef UFLLocalFunctions_5cc181267fddc3c844daa3bec40763e7::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_5cc181267fddc3c844daa3bec40763e7::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_5cc181267fddc3c844daa3bec40763e7_d5631e01dd7008453673985fc163d5bb.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
