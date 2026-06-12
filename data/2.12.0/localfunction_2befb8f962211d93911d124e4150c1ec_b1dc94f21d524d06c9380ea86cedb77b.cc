#ifndef GUARD_2befb8f962211d93911d124e4150c1ec
#define GUARD_2befb8f962211d93911d124e4150c1ec

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/common/exceptions.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_2befb8f962211d93911d124e4150c1ec
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffdiscrete_solution >
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
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< Coeffdiscrete_solution > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffdiscrete_solution>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffdiscrete_solution &coeffdiscrete_solution, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffdiscrete_solution >( coeffdiscrete_solution ) )
  {}

  void decouple ()
  {}

  void bind ( const IntersectionType &intersection, Side side )
  {
    BaseType::bind(intersection, side);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
  }

  void bind ( const EntityType &entity )
  {
    BaseType::bind(entity);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
  }

  void unbind ()
  {
    BaseType::unbind();
    std::get< 0 >( coefficients_ ).unbind();
  }

  template< class Point >
  void evaluate ( const Point &x, typename FunctionSpaceType::RangeType &result ) const
  {
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::sin( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = std::pow( tmp2, 2 );
      assert(tmp3 == tmp3);
    const auto tmp4 = 3.141592653589793 * tmp0[ 0 ];
      assert(tmp4 == tmp4);
    const auto tmp5 = std::sin( tmp4 );
      assert(tmp5 == tmp5);
    const auto tmp6 = std::pow( tmp5, 2 );
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 * tmp3;
      assert(tmp7 == tmp7);
    const auto tmp8 = -1 * tmp7;
      assert(tmp8 == tmp8);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp9 = evaluateCoefficient< 0 >( x );
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9[ 0 ] + tmp8;
      assert(tmp10 == tmp10);
    const auto tmp11 = std::pow( tmp10, 2 );
      assert(tmp11 == tmp11);
    result[ 0 ] = tmp11;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::sin( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = std::pow( tmp2, 2 );
      assert(tmp3 == tmp3);
    const auto tmp4 = 3.141592653589793 * tmp0[ 0 ];
      assert(tmp4 == tmp4);
    const auto tmp5 = std::sin( tmp4 );
      assert(tmp5 == tmp5);
    const auto tmp6 = std::cos( tmp4 );
      assert(tmp6 == tmp6);
    const auto tmp7 = 3.141592653589793 * tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = 2 * tmp7;
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp8 * tmp5;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 * tmp3;
      assert(tmp10 == tmp10);
    const auto tmp11 = -1 * tmp10;
      assert(tmp11 == tmp11);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp12 = jacobianCoefficient< 0 >( x );
      assert(tmp12 == tmp12);
    const auto tmp13 = (tmp12[ 0 ])[ 0 ] + tmp11;
      assert(tmp13 == tmp13);
    const auto tmp14 = 2 * tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = std::pow( tmp5, 2 );
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 * tmp3;
      assert(tmp16 == tmp16);
    const auto tmp17 = -1 * tmp16;
      assert(tmp17 == tmp17);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp18 = evaluateCoefficient< 0 >( x );
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18[ 0 ] + tmp17;
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp19 * tmp14;
      assert(tmp20 == tmp20);
    const auto tmp21 = std::cos( tmp1 );
      assert(tmp21 == tmp21);
    const auto tmp22 = 3.141592653589793 * tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = 2 * tmp22;
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp23 * tmp2;
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp24 * tmp15;
      assert(tmp25 == tmp25);
    const auto tmp26 = -1 * tmp25;
      assert(tmp26 == tmp26);
    const auto tmp27 = (tmp12[ 0 ])[ 1 ] + tmp26;
      assert(tmp27 == tmp27);
    const auto tmp28 = 2 * tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp19 * tmp28;
      assert(tmp29 == tmp29);
    (result[ 0 ])[ 0 ] = tmp20;
    (result[ 0 ])[ 1 ] = tmp29;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::sin( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = std::pow( tmp2, 2 );
      assert(tmp3 == tmp3);
    const auto tmp4 = 3.141592653589793 * tmp0[ 0 ];
      assert(tmp4 == tmp4);
    const auto tmp5 = std::sin( tmp4 );
      assert(tmp5 == tmp5);
    const auto tmp6 = std::cos( tmp4 );
      assert(tmp6 == tmp6);
    const auto tmp7 = 3.141592653589793 * tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = 2 * tmp7;
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp8 * tmp5;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 * tmp3;
      assert(tmp10 == tmp10);
    const auto tmp11 = -1 * tmp10;
      assert(tmp11 == tmp11);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp12 = jacobianCoefficient< 0 >( x );
      assert(tmp12 == tmp12);
    const auto tmp13 = (tmp12[ 0 ])[ 0 ] + tmp11;
      assert(tmp13 == tmp13);
    const auto tmp14 = 2 * tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp13 * tmp14;
      assert(tmp15 == tmp15);
    const auto tmp16 = -1 * tmp5;
      assert(tmp16 == tmp16);
    const auto tmp17 = 3.141592653589793 * tmp16;
      assert(tmp17 == tmp17);
    const auto tmp18 = 3.141592653589793 * tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = 2 * tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp19 * tmp5;
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp8 * tmp7;
      assert(tmp21 == tmp21);
    const auto tmp22 = tmp21 + tmp20;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 * tmp3;
      assert(tmp23 == tmp23);
    const auto tmp24 = -1 * tmp23;
      assert(tmp24 == tmp24);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp25 = hessianCoefficient< 0 >( x );
      assert(tmp25 == tmp25);
    const auto tmp26 = ((tmp25[ 0 ])[ 0 ])[ 0 ] + tmp24;
      assert(tmp26 == tmp26);
    const auto tmp27 = 2 * tmp26;
      assert(tmp27 == tmp27);
    const auto tmp28 = std::pow( tmp5, 2 );
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp28 * tmp3;
      assert(tmp29 == tmp29);
    const auto tmp30 = -1 * tmp29;
      assert(tmp30 == tmp30);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp31 = evaluateCoefficient< 0 >( x );
      assert(tmp31 == tmp31);
    const auto tmp32 = tmp31[ 0 ] + tmp30;
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp32 * tmp27;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 + tmp15;
      assert(tmp34 == tmp34);
    const auto tmp35 = std::cos( tmp1 );
      assert(tmp35 == tmp35);
    const auto tmp36 = 3.141592653589793 * tmp35;
      assert(tmp36 == tmp36);
    const auto tmp37 = 2 * tmp36;
      assert(tmp37 == tmp37);
    const auto tmp38 = tmp37 * tmp2;
      assert(tmp38 == tmp38);
    const auto tmp39 = tmp38 * tmp28;
      assert(tmp39 == tmp39);
    const auto tmp40 = -1 * tmp39;
      assert(tmp40 == tmp40);
    const auto tmp41 = (tmp12[ 0 ])[ 1 ] + tmp40;
      assert(tmp41 == tmp41);
    const auto tmp42 = tmp41 * tmp14;
      assert(tmp42 == tmp42);
    const auto tmp43 = tmp9 * tmp38;
      assert(tmp43 == tmp43);
    const auto tmp44 = -1 * tmp43;
      assert(tmp44 == tmp44);
    const auto tmp45 = ((tmp25[ 0 ])[ 0 ])[ 1 ] + tmp44;
      assert(tmp45 == tmp45);
    const auto tmp46 = 2 * tmp45;
      assert(tmp46 == tmp46);
    const auto tmp47 = tmp32 * tmp46;
      assert(tmp47 == tmp47);
    const auto tmp48 = tmp47 + tmp42;
      assert(tmp48 == tmp48);
    const auto tmp49 = 2 * tmp41;
      assert(tmp49 == tmp49);
    const auto tmp50 = tmp13 * tmp49;
      assert(tmp50 == tmp50);
    const auto tmp51 = ((tmp25[ 0 ])[ 1 ])[ 0 ] + tmp44;
      assert(tmp51 == tmp51);
    const auto tmp52 = 2 * tmp51;
      assert(tmp52 == tmp52);
    const auto tmp53 = tmp32 * tmp52;
      assert(tmp53 == tmp53);
    const auto tmp54 = tmp53 + tmp50;
      assert(tmp54 == tmp54);
    const auto tmp55 = tmp41 * tmp49;
      assert(tmp55 == tmp55);
    const auto tmp56 = -1 * tmp2;
      assert(tmp56 == tmp56);
    const auto tmp57 = 3.141592653589793 * tmp56;
      assert(tmp57 == tmp57);
    const auto tmp58 = 3.141592653589793 * tmp57;
      assert(tmp58 == tmp58);
    const auto tmp59 = 2 * tmp58;
      assert(tmp59 == tmp59);
    const auto tmp60 = tmp59 * tmp2;
      assert(tmp60 == tmp60);
    const auto tmp61 = tmp37 * tmp36;
      assert(tmp61 == tmp61);
    const auto tmp62 = tmp61 + tmp60;
      assert(tmp62 == tmp62);
    const auto tmp63 = tmp62 * tmp28;
      assert(tmp63 == tmp63);
    const auto tmp64 = -1 * tmp63;
      assert(tmp64 == tmp64);
    const auto tmp65 = ((tmp25[ 0 ])[ 1 ])[ 1 ] + tmp64;
      assert(tmp65 == tmp65);
    const auto tmp66 = 2 * tmp65;
      assert(tmp66 == tmp66);
    const auto tmp67 = tmp32 * tmp66;
      assert(tmp67 == tmp67);
    const auto tmp68 = tmp67 + tmp55;
      assert(tmp68 == tmp68);
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp34;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp48;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp54;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp68;
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

  template< std::size_t i, class Point >
  typename CoefficientFunctionSpaceType< i >::RangeType evaluateCoefficient ( const Point &x ) const
  {
    typename CoefficientFunctionSpaceType< i >::RangeType result;
    std::get< i >( coefficients_ ).evaluate( x, result );
    return result;
  }

  template< std::size_t i, class Point >
  typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
  {
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType result;
    std::get< i >( coefficients_ ).jacobian( x, result );
    return result;
  }

  template< std::size_t i, class Point >
  typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
  {
    typename CoefficientFunctionSpaceType< i >::HessianRangeType result;
    std::get< i >( coefficients_ ).hessian( x, result );
    return result;
  }
  ConstantTupleType constants_;
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffdiscrete_solution > > coefficients_;
};

} // namespace UFLLocalFunctions_2befb8f962211d93911d124e4150c1ec

PYBIND11_MODULE( localfunction_2befb8f962211d93911d124e4150c1ec_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef UFLLocalFunctions_2befb8f962211d93911d124e4150c1ec::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_2befb8f962211d93911d124e4150c1ec::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_2befb8f962211d93911d124e4150c1ec_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffdiscrete_solution ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>(gridView),name,order, coeffdiscrete_solution); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
