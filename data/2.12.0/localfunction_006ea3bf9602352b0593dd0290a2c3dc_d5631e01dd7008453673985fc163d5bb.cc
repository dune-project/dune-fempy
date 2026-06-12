#ifndef GUARD_006ea3bf9602352b0593dd0290a2c3dc
#define GUARD_006ea3bf9602352b0593dd0290a2c3dc

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

namespace UFLLocalFunctions_006ea3bf9602352b0593dd0290a2c3dc
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffsolution >
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
  typedef double Coneps;
  typedef std::tuple< std::shared_ptr< Coneps > > ConstantTupleType;
  template< std::size_t i >
  using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
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
  {
    std::get< 0 >( constants_ ) = std::make_shared< Coneps >( (Coneps(0)) );
  }

  void decouple ()
  {
    std::get< 0 >( constants_ ) = std::make_shared< Coneps >( (constant<0>()) );
  }

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
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = std::pow( tmp0[ 0 ], 2 );
      assert(tmp1 == tmp1);
    const auto tmp2 = -1 + tmp1;
      assert(tmp2 == tmp2);
    const auto tmp3 = std::pow( tmp2, 2 );
      assert(tmp3 == tmp3);
    const auto tmp4 = 0.25 * tmp3;
      assert(tmp4 == tmp4);
    double tmp5 = constant< 0 >();
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp5 * tmp5;
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 / 2;
      assert(tmp7 == tmp7);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp8 = jacobianCoefficient< 0 >( x );
      assert(tmp8 == tmp8);
    const auto tmp9 = (tmp8[ 0 ])[ 1 ] * (tmp8[ 0 ])[ 1 ];
      assert(tmp9 == tmp9);
    const auto tmp10 = (tmp8[ 0 ])[ 0 ] * (tmp8[ 0 ])[ 0 ];
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 + tmp9;
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 * tmp7;
      assert(tmp12 == tmp12);
    const auto tmp13 = tmp12 + tmp4;
      assert(tmp13 == tmp13);
    result[ 0 ] = tmp13;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::pow;
    double tmp0 = constant< 0 >();
      assert(tmp0 == tmp0);
    const auto tmp1 = tmp0 * tmp0;
      assert(tmp1 == tmp1);
    const auto tmp2 = tmp1 / 2;
      assert(tmp2 == tmp2);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp3 = hessianCoefficient< 0 >( x );
      assert(tmp3 == tmp3);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp4 = jacobianCoefficient< 0 >( x );
      assert(tmp4 == tmp4);
    const auto tmp5 = (tmp4[ 0 ])[ 1 ] * ((tmp3[ 0 ])[ 1 ])[ 0 ];
      assert(tmp5 == tmp5);
    const auto tmp6 = ((tmp3[ 0 ])[ 1 ])[ 0 ] * (tmp4[ 0 ])[ 1 ];
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp6 + tmp5;
      assert(tmp7 == tmp7);
    const auto tmp8 = (tmp4[ 0 ])[ 0 ] * ((tmp3[ 0 ])[ 0 ])[ 0 ];
      assert(tmp8 == tmp8);
    const auto tmp9 = ((tmp3[ 0 ])[ 0 ])[ 0 ] * (tmp4[ 0 ])[ 0 ];
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 + tmp8;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 + tmp7;
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 * tmp2;
      assert(tmp12 == tmp12);
    const auto tmp13 = 2 * (tmp4[ 0 ])[ 0 ];
      assert(tmp13 == tmp13);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp14 = evaluateCoefficient< 0 >( x );
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14[ 0 ] * tmp13;
      assert(tmp15 == tmp15);
    const auto tmp16 = 2 * tmp15;
      assert(tmp16 == tmp16);
    const auto tmp17 = std::pow( tmp14[ 0 ], 2 );
      assert(tmp17 == tmp17);
    const auto tmp18 = -1 + tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 * tmp16;
      assert(tmp19 == tmp19);
    const auto tmp20 = 0.25 * tmp19;
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp20 + tmp12;
      assert(tmp21 == tmp21);
    const auto tmp22 = (tmp4[ 0 ])[ 1 ] * ((tmp3[ 0 ])[ 1 ])[ 1 ];
      assert(tmp22 == tmp22);
    const auto tmp23 = ((tmp3[ 0 ])[ 1 ])[ 1 ] * (tmp4[ 0 ])[ 1 ];
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp23 + tmp22;
      assert(tmp24 == tmp24);
    const auto tmp25 = (tmp4[ 0 ])[ 0 ] * ((tmp3[ 0 ])[ 0 ])[ 1 ];
      assert(tmp25 == tmp25);
    const auto tmp26 = ((tmp3[ 0 ])[ 0 ])[ 1 ] * (tmp4[ 0 ])[ 0 ];
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp26 + tmp25;
      assert(tmp27 == tmp27);
    const auto tmp28 = tmp27 + tmp24;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp28 * tmp2;
      assert(tmp29 == tmp29);
    const auto tmp30 = 2 * (tmp4[ 0 ])[ 1 ];
      assert(tmp30 == tmp30);
    const auto tmp31 = tmp14[ 0 ] * tmp30;
      assert(tmp31 == tmp31);
    const auto tmp32 = 2 * tmp31;
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp18 * tmp32;
      assert(tmp33 == tmp33);
    const auto tmp34 = 0.25 * tmp33;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp34 + tmp29;
      assert(tmp35 == tmp35);
    (result[ 0 ])[ 0 ] = tmp21;
    (result[ 0 ])[ 1 ] = tmp35;
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

  const Coneps &coneps () const
  {
    return *std::get< 0 >( constants_ );
  }

  Coneps &coneps ()
  {
    return *std::get< 0 >( constants_ );
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution > > coefficients_;
};

} // namespace UFLLocalFunctions_006ea3bf9602352b0593dd0290a2c3dc

PYBIND11_MODULE( localfunction_006ea3bf9602352b0593dd0290a2c3dc_d5631e01dd7008453673985fc163d5bb, module )
{
  typedef UFLLocalFunctions_006ea3bf9602352b0593dd0290a2c3dc::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_006ea3bf9602352b0593dd0290a2c3dc::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_006ea3bf9602352b0593dd0290a2c3dc_d5631e01dd7008453673985fc163d5bb.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffsolution ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>(gridView),name,order, coeffsolution); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property( "eps", [] ( LocalFunctionType &self ) -> UFLLocalFunctions_006ea3bf9602352b0593dd0290a2c3dc::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > >::Coneps { return self.coneps(); }, [] ( LocalFunctionType &self, const UFLLocalFunctions_006ea3bf9602352b0593dd0290a2c3dc::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > >::Coneps &v ) { self.coneps() = v; } );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
