#ifndef GUARD_1d8980fcdff5f6568291e9d2f08eb3c8
#define GUARD_1d8980fcdff5f6568291e9d2f08eb3c8

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

namespace UFLLocalFunctions_1d8980fcdff5f6568291e9d2f08eb3c8
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
    const auto tmp6 = -1 * tmp5;
      assert(tmp6 == tmp6);
    const auto tmp7 = 3.141592653589793 * tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = 3.141592653589793 * tmp7;
      assert(tmp8 == tmp8);
    const auto tmp9 = 2 * tmp8;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 * tmp5;
      assert(tmp10 == tmp10);
    const auto tmp11 = std::cos( tmp4 );
      assert(tmp11 == tmp11);
    const auto tmp12 = 3.141592653589793 * tmp11;
      assert(tmp12 == tmp12);
    const auto tmp13 = 2 * tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = tmp13 * tmp12;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 + tmp10;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 * tmp3;
      assert(tmp16 == tmp16);
    const auto tmp17 = -1 * tmp16;
      assert(tmp17 == tmp17);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp18 = hessianCoefficient< 0 >( x );
      assert(tmp18 == tmp18);
    const auto tmp19 = ((tmp18[ 0 ])[ 0 ])[ 0 ] + tmp17;
      assert(tmp19 == tmp19);
    const auto tmp20 = std::pow( tmp5, 2 );
      assert(tmp20 == tmp20);
    const auto tmp21 = -1 * tmp2;
      assert(tmp21 == tmp21);
    const auto tmp22 = 3.141592653589793 * tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = 3.141592653589793 * tmp22;
      assert(tmp23 == tmp23);
    const auto tmp24 = 2 * tmp23;
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp24 * tmp2;
      assert(tmp25 == tmp25);
    const auto tmp26 = std::cos( tmp1 );
      assert(tmp26 == tmp26);
    const auto tmp27 = 3.141592653589793 * tmp26;
      assert(tmp27 == tmp27);
    const auto tmp28 = 2 * tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp28 * tmp27;
      assert(tmp29 == tmp29);
    const auto tmp30 = tmp29 + tmp25;
      assert(tmp30 == tmp30);
    const auto tmp31 = tmp30 * tmp20;
      assert(tmp31 == tmp31);
    const auto tmp32 = -1 * tmp31;
      assert(tmp32 == tmp32);
    const auto tmp33 = ((tmp18[ 0 ])[ 1 ])[ 1 ] + tmp32;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 + tmp19;
      assert(tmp34 == tmp34);
    const auto tmp35 = std::pow( tmp34, 2 );
      assert(tmp35 == tmp35);
    result[ 0 ] = tmp35;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    DUNE_THROW(Dune::NotImplemented,"jacobian method could not be generated for local function (TooHighDerivative('CodeGenerator does not allow for third order derivatives, yet.')).");
    result=typename FunctionSpaceType::JacobianRangeType(0);
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

} // namespace UFLLocalFunctions_1d8980fcdff5f6568291e9d2f08eb3c8

PYBIND11_MODULE( localfunction_1d8980fcdff5f6568291e9d2f08eb3c8_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef UFLLocalFunctions_1d8980fcdff5f6568291e9d2f08eb3c8::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_1d8980fcdff5f6568291e9d2f08eb3c8::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_1d8980fcdff5f6568291e9d2f08eb3c8_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffdiscrete_solution ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>(gridView),name,order, coeffdiscrete_solution); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
