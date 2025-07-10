#ifndef GUARD_9286a322b9092ac995ff923770a7c06b
#define GUARD_9286a322b9092ac995ff923770a7c06b

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

namespace UFLLocalFunctions_9286a322b9092ac995ff923770a7c06b
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
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
    const auto tmp2 = std::sin( tmp1 );
    const auto tmp3 = 3.141592653589793 * tmp0[ 0 ];
    const auto tmp4 = std::cos( tmp3 );
    const auto tmp5 = 3.141592653589793 * tmp4;
    const auto tmp6 = tmp5 * tmp2;
    const auto tmp7 = -1 * tmp6;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp8 = jacobianCoefficient< 0 >( x );
    const auto tmp9 = (tmp8[ 0 ])[ 0 ] + tmp7;
    const auto tmp10 = tmp9 * tmp9;
    const auto tmp11 = std::sin( tmp3 );
    const auto tmp12 = std::cos( tmp1 );
    const auto tmp13 = 3.141592653589793 * tmp12;
    const auto tmp14 = tmp13 * tmp11;
    const auto tmp15 = -1 * tmp14;
    const auto tmp16 = (tmp8[ 0 ])[ 1 ] + tmp15;
    const auto tmp17 = tmp16 * tmp16;
    const auto tmp18 = tmp17 + tmp10;
    result[ 0 ] = tmp18;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
    const auto tmp2 = std::sin( tmp1 );
    const auto tmp3 = 3.141592653589793 * tmp0[ 0 ];
    const auto tmp4 = std::cos( tmp3 );
    const auto tmp5 = 3.141592653589793 * tmp4;
    const auto tmp6 = tmp5 * tmp2;
    const auto tmp7 = -1 * tmp6;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp8 = jacobianCoefficient< 0 >( x );
    const auto tmp9 = (tmp8[ 0 ])[ 0 ] + tmp7;
    const auto tmp10 = std::sin( tmp3 );
    const auto tmp11 = -1 * tmp10;
    const auto tmp12 = 3.141592653589793 * tmp11;
    const auto tmp13 = 3.141592653589793 * tmp12;
    const auto tmp14 = tmp13 * tmp2;
    const auto tmp15 = -1 * tmp14;
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp16 = hessianCoefficient< 0 >( x );
    const auto tmp17 = ((tmp16[ 0 ])[ 0 ])[ 0 ] + tmp15;
    const auto tmp18 = tmp17 * tmp9;
    const auto tmp19 = tmp9 * tmp17;
    const auto tmp20 = tmp19 + tmp18;
    const auto tmp21 = std::cos( tmp1 );
    const auto tmp22 = 3.141592653589793 * tmp21;
    const auto tmp23 = tmp22 * tmp10;
    const auto tmp24 = -1 * tmp23;
    const auto tmp25 = (tmp8[ 0 ])[ 1 ] + tmp24;
    const auto tmp26 = tmp5 * tmp22;
    const auto tmp27 = -1 * tmp26;
    const auto tmp28 = ((tmp16[ 0 ])[ 1 ])[ 0 ] + tmp27;
    const auto tmp29 = tmp28 * tmp25;
    const auto tmp30 = tmp25 * tmp28;
    const auto tmp31 = tmp30 + tmp29;
    const auto tmp32 = tmp31 + tmp20;
    const auto tmp33 = ((tmp16[ 0 ])[ 0 ])[ 1 ] + tmp27;
    const auto tmp34 = tmp33 * tmp9;
    const auto tmp35 = tmp9 * tmp33;
    const auto tmp36 = tmp35 + tmp34;
    const auto tmp37 = -1 * tmp2;
    const auto tmp38 = 3.141592653589793 * tmp37;
    const auto tmp39 = 3.141592653589793 * tmp38;
    const auto tmp40 = tmp39 * tmp10;
    const auto tmp41 = -1 * tmp40;
    const auto tmp42 = ((tmp16[ 0 ])[ 1 ])[ 1 ] + tmp41;
    const auto tmp43 = tmp42 * tmp25;
    const auto tmp44 = tmp25 * tmp42;
    const auto tmp45 = tmp44 + tmp43;
    const auto tmp46 = tmp45 + tmp36;
    (result[ 0 ])[ 0 ] = tmp32;
    (result[ 0 ])[ 1 ] = tmp46;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffdiscrete_solution > > coefficients_;
};

} // namespace UFLLocalFunctions_9286a322b9092ac995ff923770a7c06b

PYBIND11_MODULE( localfunction_9286a322b9092ac995ff923770a7c06b_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef UFLLocalFunctions_9286a322b9092ac995ff923770a7c06b::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_9286a322b9092ac995ff923770a7c06b::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_9286a322b9092ac995ff923770a7c06b_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffdiscrete_solution ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>(gridView),name,order, coeffdiscrete_solution); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
