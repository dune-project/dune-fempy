#ifndef GUARD_aad4e337e18f41b1a1c5da88cf44b28e
#define GUARD_aad4e337e18f41b1a1c5da88cf44b28e

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
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_aad4e337e18f41b1a1c5da88cf44b28e
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class CoeffgradSol, class Coeffsol >
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
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< CoeffgradSol, Coeffsol > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<CoeffgradSol>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsol>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const CoeffgradSol &coeffgradSol, const Coeffsol &coeffsol, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< CoeffgradSol >( coeffgradSol ), Dune::Fem::ConstLocalFunction< Coeffsol >( coeffsol ) )
  {}

  void bind ( const IntersectionType &intersection, Side side )
  {
    BaseType::bind(intersection, side);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
    std::get< 1 >( coefficients_ ).bind( this->entity() );
  }

  void bind ( const EntityType &entity )
  {
    BaseType::bind(entity);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
    std::get< 1 >( coefficients_ ).bind( this->entity() );
  }

  void unbind ()
  {
    BaseType::unbind();
    std::get< 0 >( coefficients_ ).unbind();
    std::get< 1 >( coefficients_ ).unbind();
  }

  template< class Point >
  void evaluate ( const Point &x, typename FunctionSpaceType::RangeType &result ) const
  {
    using std::cos;
    using std::pow;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
    const auto tmp2 = std::cos( tmp1 );
    const auto tmp3 = 3.141592653589793 * tmp0[ 0 ];
    const auto tmp4 = std::cos( tmp3 );
    const auto tmp5 = tmp4 * tmp2;
    const auto tmp6 = -1 * tmp5;
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp7 = evaluateCoefficient< 1 >( x );
    const auto tmp8 = tmp7[ 0 ] + tmp6;
    const auto tmp9 = std::pow( tmp8, 2 );
    const auto tmp10 = std::sin( tmp3 );
    const auto tmp11 = -1 * tmp10;
    const auto tmp12 = 3.141592653589793 * tmp11;
    const auto tmp13 = tmp12 * tmp2;
    const auto tmp14 = -1 * tmp13;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp15 = evaluateCoefficient< 0 >( x );
    const auto tmp16 = tmp15[ 0 ] + tmp14;
    const auto tmp17 = tmp16 * tmp16;
    const auto tmp18 = std::sin( tmp1 );
    const auto tmp19 = -1 * tmp18;
    const auto tmp20 = 3.141592653589793 * tmp19;
    const auto tmp21 = tmp20 * tmp4;
    const auto tmp22 = -1 * tmp21;
    const auto tmp23 = tmp15[ 1 ] + tmp22;
    const auto tmp24 = tmp23 * tmp23;
    const auto tmp25 = tmp24 + tmp17;
    const auto tmp26 = 3.141592653589793 * tmp4;
    const auto tmp27 = -1 * tmp26;
    const auto tmp28 = 3.141592653589793 * tmp27;
    const auto tmp29 = tmp28 * tmp2;
    const auto tmp30 = -1 * tmp29;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp31 = jacobianCoefficient< 0 >( x );
    const auto tmp32 = (tmp31[ 0 ])[ 0 ] + tmp30;
    const auto tmp33 = 3.141592653589793 * tmp2;
    const auto tmp34 = -1 * tmp33;
    const auto tmp35 = 3.141592653589793 * tmp34;
    const auto tmp36 = tmp35 * tmp4;
    const auto tmp37 = -1 * tmp36;
    const auto tmp38 = (tmp31[ 1 ])[ 1 ] + tmp37;
    const auto tmp39 = tmp38 + tmp32;
    const auto tmp40 = std::pow( tmp39, 2 );
    result[ 0 ] = tmp9;
    result[ 1 ] = tmp25;
    result[ 2 ] = tmp40;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
    const auto tmp2 = std::cos( tmp1 );
    const auto tmp3 = 3.141592653589793 * tmp0[ 0 ];
    const auto tmp4 = std::sin( tmp3 );
    const auto tmp5 = -1 * tmp4;
    const auto tmp6 = 3.141592653589793 * tmp5;
    const auto tmp7 = tmp6 * tmp2;
    const auto tmp8 = -1 * tmp7;
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp9 = jacobianCoefficient< 1 >( x );
    const auto tmp10 = (tmp9[ 0 ])[ 0 ] + tmp8;
    const auto tmp11 = 2 * tmp10;
    const auto tmp12 = std::cos( tmp3 );
    const auto tmp13 = tmp12 * tmp2;
    const auto tmp14 = -1 * tmp13;
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp15 = evaluateCoefficient< 1 >( x );
    const auto tmp16 = tmp15[ 0 ] + tmp14;
    const auto tmp17 = tmp16 * tmp11;
    const auto tmp18 = std::sin( tmp1 );
    const auto tmp19 = -1 * tmp18;
    const auto tmp20 = 3.141592653589793 * tmp19;
    const auto tmp21 = tmp20 * tmp12;
    const auto tmp22 = -1 * tmp21;
    const auto tmp23 = (tmp9[ 0 ])[ 1 ] + tmp22;
    const auto tmp24 = 2 * tmp23;
    const auto tmp25 = tmp16 * tmp24;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp26 = evaluateCoefficient< 0 >( x );
    const auto tmp27 = tmp26[ 0 ] + tmp8;
    const auto tmp28 = 3.141592653589793 * tmp12;
    const auto tmp29 = -1 * tmp28;
    const auto tmp30 = 3.141592653589793 * tmp29;
    const auto tmp31 = tmp30 * tmp2;
    const auto tmp32 = -1 * tmp31;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp33 = jacobianCoefficient< 0 >( x );
    const auto tmp34 = (tmp33[ 0 ])[ 0 ] + tmp32;
    const auto tmp35 = tmp34 * tmp27;
    const auto tmp36 = tmp27 * tmp34;
    const auto tmp37 = tmp36 + tmp35;
    const auto tmp38 = tmp26[ 1 ] + tmp22;
    const auto tmp39 = tmp6 * tmp20;
    const auto tmp40 = -1 * tmp39;
    const auto tmp41 = (tmp33[ 1 ])[ 0 ] + tmp40;
    const auto tmp42 = tmp41 * tmp38;
    const auto tmp43 = tmp38 * tmp41;
    const auto tmp44 = tmp43 + tmp42;
    const auto tmp45 = tmp44 + tmp37;
    const auto tmp46 = (tmp33[ 0 ])[ 1 ] + tmp40;
    const auto tmp47 = tmp46 * tmp27;
    const auto tmp48 = tmp27 * tmp46;
    const auto tmp49 = tmp48 + tmp47;
    const auto tmp50 = 3.141592653589793 * tmp2;
    const auto tmp51 = -1 * tmp50;
    const auto tmp52 = 3.141592653589793 * tmp51;
    const auto tmp53 = tmp52 * tmp12;
    const auto tmp54 = -1 * tmp53;
    const auto tmp55 = (tmp33[ 1 ])[ 1 ] + tmp54;
    const auto tmp56 = tmp55 * tmp38;
    const auto tmp57 = tmp38 * tmp55;
    const auto tmp58 = tmp57 + tmp56;
    const auto tmp59 = tmp58 + tmp49;
    const auto tmp60 = 3.141592653589793 * tmp6;
    const auto tmp61 = -1 * tmp60;
    const auto tmp62 = 3.141592653589793 * tmp61;
    const auto tmp63 = tmp62 * tmp2;
    const auto tmp64 = -1 * tmp63;
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp65 = hessianCoefficient< 0 >( x );
    const auto tmp66 = ((tmp65[ 0 ])[ 0 ])[ 0 ] + tmp64;
    const auto tmp67 = tmp52 * tmp6;
    const auto tmp68 = -1 * tmp67;
    const auto tmp69 = ((tmp65[ 1 ])[ 1 ])[ 0 ] + tmp68;
    const auto tmp70 = tmp69 + tmp66;
    const auto tmp71 = 2 * tmp70;
    const auto tmp72 = tmp55 + tmp34;
    const auto tmp73 = tmp72 * tmp71;
    const auto tmp74 = 3.141592653589793 * tmp20;
    const auto tmp75 = -1 * tmp74;
    const auto tmp76 = 3.141592653589793 * tmp75;
    const auto tmp77 = tmp76 * tmp12;
    const auto tmp78 = -1 * tmp77;
    const auto tmp79 = ((tmp65[ 1 ])[ 1 ])[ 1 ] + tmp78;
    const auto tmp80 = tmp30 * tmp20;
    const auto tmp81 = -1 * tmp80;
    const auto tmp82 = ((tmp65[ 0 ])[ 0 ])[ 1 ] + tmp81;
    const auto tmp83 = tmp82 + tmp79;
    const auto tmp84 = 2 * tmp83;
    const auto tmp85 = tmp72 * tmp84;
    (result[ 0 ])[ 0 ] = tmp17;
    (result[ 0 ])[ 1 ] = tmp25;
    (result[ 1 ])[ 0 ] = tmp45;
    (result[ 1 ])[ 1 ] = tmp59;
    (result[ 2 ])[ 0 ] = tmp73;
    (result[ 2 ])[ 1 ] = tmp85;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    DUNE_THROW(Dune::NotImplemented,"hessian method could not be generated for local function (TooHighDerivative('CodeGenerator does not allow for third order derivatives, yet.'))");
    result=typename FunctionSpaceType::HessianRangeType(0);
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
    std::get< i >( coefficients_ ).evaluate( x, result );;
    return result;
  }

  template< std::size_t i, class Point >
  typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
  {
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType result;
    std::get< i >( coefficients_ ).jacobian( x, result );;
    return result;
  }

  template< std::size_t i, class Point >
  typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
  {
    typename CoefficientFunctionSpaceType< i >::HessianRangeType result;
    std::get< i >( coefficients_ ).hessian( x, result );;
    return result;
  }
  ConstantTupleType constants_;
  std::tuple< Dune::Fem::ConstLocalFunction< CoeffgradSol >, Dune::Fem::ConstLocalFunction< Coeffsol > > coefficients_;
};

} // namespace UFLLocalFunctions_aad4e337e18f41b1a1c5da88cf44b28e

PYBIND11_MODULE( localfunction_aad4e337e18f41b1a1c5da88cf44b28e_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef UFLLocalFunctions_aad4e337e18f41b1a1c5da88cf44b28e::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_aad4e337e18f41b1a1c5da88cf44b28e::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_aad4e337e18f41b1a1c5da88cf44b28e_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > > &coeffgradSol, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffsol ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>(gridView),name,order, coeffgradSol, coeffsol); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
