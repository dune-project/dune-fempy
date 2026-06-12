#ifndef GUARD_130c7f5da9879f68143b57057f95bb49
#define GUARD_130c7f5da9879f68143b57057f95bb49

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

namespace UFLLocalFunctions_130c7f5da9879f68143b57057f95bb49
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
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = tmp0[ 1 ] * tmp1;
      assert(tmp2 == tmp2);
    const auto tmp3 = std::cos( tmp2 );
      assert(tmp3 == tmp3);
    const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp4 * tmp3;
      assert(tmp5 == tmp5);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp6 = evaluateCoefficient< 0 >( x );
      assert(tmp6 == tmp6);
    const auto tmp7 = -1 * tmp6[ 0 ];
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 + tmp5;
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp8 * tmp8;
      assert(tmp9 == tmp9);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp10 = jacobianCoefficient< 0 >( x );
      assert(tmp10 == tmp10);
    const auto tmp11 = -1 * (tmp10[ 0 ])[ 1 ];
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp0[ 0 ] * tmp3;
      assert(tmp12 == tmp12);
    const auto tmp13 = std::sin( tmp2 );
      assert(tmp13 == tmp13);
    const auto tmp14 = -1 * tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp1 * tmp14;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp4 * tmp15;
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp16 + tmp12;
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp17 + tmp11;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 * tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = -1 * (tmp10[ 0 ])[ 0 ];
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp0[ 1 ] * tmp3;
      assert(tmp21 == tmp21);
    const auto tmp22 = 3.141592653589793 * tmp0[ 1 ];
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 * tmp14;
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp4 * tmp23;
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp24 + tmp21;
      assert(tmp25 == tmp25);
    const auto tmp26 = tmp25 + tmp20;
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp26 * tmp26;
      assert(tmp27 == tmp27);
    const auto tmp28 = tmp27 + tmp19;
      assert(tmp28 == tmp28);
    result[ 0 ] = tmp9;
    result[ 1 ] = tmp28;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = tmp0[ 1 ] * tmp1;
      assert(tmp2 == tmp2);
    const auto tmp3 = std::cos( tmp2 );
      assert(tmp3 == tmp3);
    const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp4 * tmp3;
      assert(tmp5 == tmp5);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp6 = evaluateCoefficient< 0 >( x );
      assert(tmp6 == tmp6);
    const auto tmp7 = -1 * tmp6[ 0 ];
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 + tmp5;
      assert(tmp8 == tmp8);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp9 = jacobianCoefficient< 0 >( x );
      assert(tmp9 == tmp9);
    const auto tmp10 = -1 * (tmp9[ 0 ])[ 0 ];
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp0[ 1 ] * tmp3;
      assert(tmp11 == tmp11);
    const auto tmp12 = std::sin( tmp2 );
      assert(tmp12 == tmp12);
    const auto tmp13 = -1 * tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = 3.141592653589793 * tmp0[ 1 ];
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 * tmp13;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp4 * tmp15;
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp16 + tmp11;
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp17 + tmp10;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 * tmp8;
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp8 * tmp18;
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp20 + tmp19;
      assert(tmp21 == tmp21);
    const auto tmp22 = -1 * (tmp9[ 0 ])[ 1 ];
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp0[ 0 ] * tmp3;
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp1 * tmp13;
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp4 * tmp24;
      assert(tmp25 == tmp25);
    const auto tmp26 = tmp25 + tmp23;
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp26 + tmp22;
      assert(tmp27 == tmp27);
    const auto tmp28 = tmp27 * tmp8;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp8 * tmp27;
      assert(tmp29 == tmp29);
    const auto tmp30 = tmp29 + tmp28;
      assert(tmp30 == tmp30);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp31 = hessianCoefficient< 0 >( x );
      assert(tmp31 == tmp31);
    const auto tmp32 = -1 * ((tmp31[ 0 ])[ 1 ])[ 0 ];
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp0[ 0 ] * tmp15;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 + tmp3;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp0[ 1 ] * tmp24;
      assert(tmp35 == tmp35);
    const auto tmp36 = 3.141592653589793 * tmp13;
      assert(tmp36 == tmp36);
    const auto tmp37 = tmp14 * tmp3;
      assert(tmp37 == tmp37);
    const auto tmp38 = -1 * tmp37;
      assert(tmp38 == tmp38);
    const auto tmp39 = tmp1 * tmp38;
      assert(tmp39 == tmp39);
    const auto tmp40 = tmp39 + tmp36;
      assert(tmp40 == tmp40);
    const auto tmp41 = tmp40 * tmp4;
      assert(tmp41 == tmp41);
    const auto tmp42 = tmp41 + tmp35;
      assert(tmp42 == tmp42);
    const auto tmp43 = tmp42 + tmp34;
      assert(tmp43 == tmp43);
    const auto tmp44 = tmp43 + tmp32;
      assert(tmp44 == tmp44);
    const auto tmp45 = tmp27 * tmp44;
      assert(tmp45 == tmp45);
    const auto tmp46 = tmp44 * tmp27;
      assert(tmp46 == tmp46);
    const auto tmp47 = tmp46 + tmp45;
      assert(tmp47 == tmp47);
    const auto tmp48 = -1 * ((tmp31[ 0 ])[ 0 ])[ 0 ];
      assert(tmp48 == tmp48);
    const auto tmp49 = tmp0[ 1 ] * tmp15;
      assert(tmp49 == tmp49);
    const auto tmp50 = tmp14 * tmp38;
      assert(tmp50 == tmp50);
    const auto tmp51 = tmp4 * tmp50;
      assert(tmp51 == tmp51);
    const auto tmp52 = tmp51 + tmp49;
      assert(tmp52 == tmp52);
    const auto tmp53 = tmp52 + tmp49;
      assert(tmp53 == tmp53);
    const auto tmp54 = tmp53 + tmp48;
      assert(tmp54 == tmp54);
    const auto tmp55 = tmp18 * tmp54;
      assert(tmp55 == tmp55);
    const auto tmp56 = tmp54 * tmp18;
      assert(tmp56 == tmp56);
    const auto tmp57 = tmp56 + tmp55;
      assert(tmp57 == tmp57);
    const auto tmp58 = tmp57 + tmp47;
      assert(tmp58 == tmp58);
    const auto tmp59 = -1 * ((tmp31[ 0 ])[ 1 ])[ 1 ];
      assert(tmp59 == tmp59);
    const auto tmp60 = tmp0[ 0 ] * tmp24;
      assert(tmp60 == tmp60);
    const auto tmp61 = tmp1 * tmp3;
      assert(tmp61 == tmp61);
    const auto tmp62 = -1 * tmp61;
      assert(tmp62 == tmp62);
    const auto tmp63 = tmp1 * tmp62;
      assert(tmp63 == tmp63);
    const auto tmp64 = tmp4 * tmp63;
      assert(tmp64 == tmp64);
    const auto tmp65 = tmp64 + tmp60;
      assert(tmp65 == tmp65);
    const auto tmp66 = tmp65 + tmp60;
      assert(tmp66 == tmp66);
    const auto tmp67 = tmp66 + tmp59;
      assert(tmp67 == tmp67);
    const auto tmp68 = tmp27 * tmp67;
      assert(tmp68 == tmp68);
    const auto tmp69 = tmp67 * tmp27;
      assert(tmp69 == tmp69);
    const auto tmp70 = tmp69 + tmp68;
      assert(tmp70 == tmp70);
    const auto tmp71 = -1 * ((tmp31[ 0 ])[ 0 ])[ 1 ];
      assert(tmp71 == tmp71);
    const auto tmp72 = tmp35 + tmp3;
      assert(tmp72 == tmp72);
    const auto tmp73 = tmp14 * tmp62;
      assert(tmp73 == tmp73);
    const auto tmp74 = tmp73 + tmp36;
      assert(tmp74 == tmp74);
    const auto tmp75 = tmp74 * tmp4;
      assert(tmp75 == tmp75);
    const auto tmp76 = tmp75 + tmp33;
      assert(tmp76 == tmp76);
    const auto tmp77 = tmp76 + tmp72;
      assert(tmp77 == tmp77);
    const auto tmp78 = tmp77 + tmp71;
      assert(tmp78 == tmp78);
    const auto tmp79 = tmp18 * tmp78;
      assert(tmp79 == tmp79);
    const auto tmp80 = tmp78 * tmp18;
      assert(tmp80 == tmp80);
    const auto tmp81 = tmp80 + tmp79;
      assert(tmp81 == tmp81);
    const auto tmp82 = tmp81 + tmp70;
      assert(tmp82 == tmp82);
    (result[ 0 ])[ 0 ] = tmp21;
    (result[ 0 ])[ 1 ] = tmp30;
    (result[ 1 ])[ 0 ] = tmp58;
    (result[ 1 ])[ 1 ] = tmp82;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution > > coefficients_;
};

} // namespace UFLLocalFunctions_130c7f5da9879f68143b57057f95bb49

PYBIND11_MODULE( localfunction_130c7f5da9879f68143b57057f95bb49_f084d41ce9ff164524e5c3ea6f8082cc, module )
{
  typedef UFLLocalFunctions_130c7f5da9879f68143b57057f95bb49::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::Grid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_130c7f5da9879f68143b57057f95bb49::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::Grid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_130c7f5da9879f68143b57057f95bb49_f084d41ce9ff164524e5c3ea6f8082cc.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::Grid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffsolution ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::Vem::Grid<2,2>::LeafGridView>(gridView),name,order, coeffsolution); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
