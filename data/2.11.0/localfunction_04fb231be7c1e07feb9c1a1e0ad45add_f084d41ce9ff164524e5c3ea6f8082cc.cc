#ifndef GUARD_04fb231be7c1e07feb9c1a1e0ad45add
#define GUARD_04fb231be7c1e07feb9c1a1e0ad45add

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
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_04fb231be7c1e07feb9c1a1e0ad45add
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffsolution >
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
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< Coeffsolution > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsolution>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffsolution &coeffsolution, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffsolution >( coeffsolution ) )
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
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
    const auto tmp2 = -1 * tmp1;
    const auto tmp3 = tmp0[ 1 ] + tmp2;
    const auto tmp4 = tmp0[ 0 ] * tmp0[ 0 ];
    const auto tmp5 = -1 * tmp4;
    const auto tmp6 = tmp0[ 0 ] + tmp5;
    const auto tmp7 = tmp6 * tmp3;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp8 = evaluateCoefficient< 0 >( x );
    const auto tmp9 = -1 * tmp8[ 0 ];
    const auto tmp10 = tmp9 + tmp7;
    const auto tmp11 = tmp10 * tmp10;
    const auto tmp12 = tmp0[ 0 ] + tmp0[ 0 ];
    const auto tmp13 = -1 * tmp12;
    const auto tmp14 = 1 + tmp13;
    const auto tmp15 = tmp14 * tmp3;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp16 = jacobianCoefficient< 0 >( x );
    const auto tmp17 = -1 * (tmp16[ 0 ])[ 0 ];
    const auto tmp18 = tmp17 + tmp15;
    const auto tmp19 = tmp18 * tmp18;
    const auto tmp20 = tmp0[ 1 ] + tmp0[ 1 ];
    const auto tmp21 = -1 * tmp20;
    const auto tmp22 = 1 + tmp21;
    const auto tmp23 = tmp22 * tmp6;
    const auto tmp24 = -1 * (tmp16[ 0 ])[ 1 ];
    const auto tmp25 = tmp24 + tmp23;
    const auto tmp26 = tmp25 * tmp25;
    const auto tmp27 = tmp26 + tmp19;
    result[ 0 ] = tmp11;
    result[ 1 ] = tmp27;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
    const auto tmp2 = -1 * tmp1;
    const auto tmp3 = tmp0[ 1 ] + tmp2;
    const auto tmp4 = tmp0[ 0 ] * tmp0[ 0 ];
    const auto tmp5 = -1 * tmp4;
    const auto tmp6 = tmp0[ 0 ] + tmp5;
    const auto tmp7 = tmp6 * tmp3;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp8 = evaluateCoefficient< 0 >( x );
    const auto tmp9 = -1 * tmp8[ 0 ];
    const auto tmp10 = tmp9 + tmp7;
    const auto tmp11 = tmp0[ 0 ] + tmp0[ 0 ];
    const auto tmp12 = -1 * tmp11;
    const auto tmp13 = 1 + tmp12;
    const auto tmp14 = tmp13 * tmp3;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp15 = jacobianCoefficient< 0 >( x );
    const auto tmp16 = -1 * (tmp15[ 0 ])[ 0 ];
    const auto tmp17 = tmp16 + tmp14;
    const auto tmp18 = tmp17 * tmp10;
    const auto tmp19 = tmp10 * tmp17;
    const auto tmp20 = tmp19 + tmp18;
    const auto tmp21 = tmp0[ 1 ] + tmp0[ 1 ];
    const auto tmp22 = -1 * tmp21;
    const auto tmp23 = 1 + tmp22;
    const auto tmp24 = tmp23 * tmp6;
    const auto tmp25 = -1 * (tmp15[ 0 ])[ 1 ];
    const auto tmp26 = tmp25 + tmp24;
    const auto tmp27 = tmp26 * tmp10;
    const auto tmp28 = tmp10 * tmp26;
    const auto tmp29 = tmp28 + tmp27;
    const auto tmp30 = -2 * tmp3;
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp31 = hessianCoefficient< 0 >( x );
    const auto tmp32 = -1 * ((tmp31[ 0 ])[ 0 ])[ 0 ];
    const auto tmp33 = tmp32 + tmp30;
    const auto tmp34 = tmp33 * tmp17;
    const auto tmp35 = tmp17 * tmp33;
    const auto tmp36 = tmp35 + tmp34;
    const auto tmp37 = tmp13 * tmp23;
    const auto tmp38 = -1 * ((tmp31[ 0 ])[ 1 ])[ 0 ];
    const auto tmp39 = tmp38 + tmp37;
    const auto tmp40 = tmp39 * tmp26;
    const auto tmp41 = tmp26 * tmp39;
    const auto tmp42 = tmp41 + tmp40;
    const auto tmp43 = tmp42 + tmp36;
    const auto tmp44 = -1 * ((tmp31[ 0 ])[ 0 ])[ 1 ];
    const auto tmp45 = tmp44 + tmp37;
    const auto tmp46 = tmp45 * tmp17;
    const auto tmp47 = tmp17 * tmp45;
    const auto tmp48 = tmp47 + tmp46;
    const auto tmp49 = -2 * tmp6;
    const auto tmp50 = -1 * ((tmp31[ 0 ])[ 1 ])[ 1 ];
    const auto tmp51 = tmp50 + tmp49;
    const auto tmp52 = tmp51 * tmp26;
    const auto tmp53 = tmp26 * tmp51;
    const auto tmp54 = tmp53 + tmp52;
    const auto tmp55 = tmp54 + tmp48;
    (result[ 0 ])[ 0 ] = tmp20;
    (result[ 0 ])[ 1 ] = tmp29;
    (result[ 1 ])[ 0 ] = tmp43;
    (result[ 1 ])[ 1 ] = tmp55;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution > > coefficients_;
};

} // namespace UFLLocalFunctions_04fb231be7c1e07feb9c1a1e0ad45add

PYBIND11_MODULE( localfunction_04fb231be7c1e07feb9c1a1e0ad45add_f084d41ce9ff164524e5c3ea6f8082cc, module )
{
  typedef UFLLocalFunctions_04fb231be7c1e07feb9c1a1e0ad45add::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::Grid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_04fb231be7c1e07feb9c1a1e0ad45add::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::Grid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_04fb231be7c1e07feb9c1a1e0ad45add_f084d41ce9ff164524e5c3ea6f8082cc.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::Grid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffsolution ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::Vem::Grid<2,2>::LeafGridView>(gridView),name,order, coeffsolution); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
