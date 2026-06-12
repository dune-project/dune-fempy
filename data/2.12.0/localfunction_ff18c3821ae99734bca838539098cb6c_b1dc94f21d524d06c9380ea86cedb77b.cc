#ifndef GUARD_ff18c3821ae99734bca838539098cb6c
#define GUARD_ff18c3821ae99734bca838539098cb6c

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

namespace UFLLocalFunctions_ff18c3821ae99734bca838539098cb6c
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

  void decouple ()
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
      assert(tmp0 == tmp0);
    const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::sin( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = 3.141592653589793 * tmp0[ 0 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = std::sin( tmp3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = tmp4 * tmp2;
      assert(tmp5 == tmp5);
    const auto tmp6 = -1 * tmp5;
      assert(tmp6 == tmp6);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp7 = evaluateCoefficient< 1 >( x );
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7[ 0 ] + tmp6;
      assert(tmp8 == tmp8);
    const auto tmp9 = std::pow( tmp8, 2 );
      assert(tmp9 == tmp9);
    const auto tmp10 = std::cos( tmp3 );
      assert(tmp10 == tmp10);
    const auto tmp11 = 3.141592653589793 * tmp10;
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 * tmp2;
      assert(tmp12 == tmp12);
    const auto tmp13 = -1 * tmp12;
      assert(tmp13 == tmp13);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp14 = evaluateCoefficient< 0 >( x );
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14[ 0 ] + tmp13;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 * tmp15;
      assert(tmp16 == tmp16);
    const auto tmp17 = std::cos( tmp1 );
      assert(tmp17 == tmp17);
    const auto tmp18 = 3.141592653589793 * tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 * tmp4;
      assert(tmp19 == tmp19);
    const auto tmp20 = -1 * tmp19;
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp14[ 1 ] + tmp20;
      assert(tmp21 == tmp21);
    const auto tmp22 = tmp21 * tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 + tmp16;
      assert(tmp23 == tmp23);
    const auto tmp24 = -1 * tmp4;
      assert(tmp24 == tmp24);
    const auto tmp25 = 3.141592653589793 * tmp24;
      assert(tmp25 == tmp25);
    const auto tmp26 = 3.141592653589793 * tmp25;
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp26 * tmp2;
      assert(tmp27 == tmp27);
    const auto tmp28 = -1 * tmp27;
      assert(tmp28 == tmp28);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp29 = jacobianCoefficient< 0 >( x );
      assert(tmp29 == tmp29);
    const auto tmp30 = (tmp29[ 0 ])[ 0 ] + tmp28;
      assert(tmp30 == tmp30);
    const auto tmp31 = -1 * tmp2;
      assert(tmp31 == tmp31);
    const auto tmp32 = 3.141592653589793 * tmp31;
      assert(tmp32 == tmp32);
    const auto tmp33 = 3.141592653589793 * tmp32;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 * tmp4;
      assert(tmp34 == tmp34);
    const auto tmp35 = -1 * tmp34;
      assert(tmp35 == tmp35);
    const auto tmp36 = (tmp29[ 1 ])[ 1 ] + tmp35;
      assert(tmp36 == tmp36);
    const auto tmp37 = tmp36 + tmp30;
      assert(tmp37 == tmp37);
    const auto tmp38 = std::pow( tmp37, 2 );
      assert(tmp38 == tmp38);
    result[ 0 ] = tmp9;
    result[ 1 ] = tmp23;
    result[ 2 ] = tmp38;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cos;
    using std::sin;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp0 == tmp0);
    const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
      assert(tmp1 == tmp1);
    const auto tmp2 = std::sin( tmp1 );
      assert(tmp2 == tmp2);
    const auto tmp3 = 3.141592653589793 * tmp0[ 0 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = std::cos( tmp3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = 3.141592653589793 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp5 * tmp2;
      assert(tmp6 == tmp6);
    const auto tmp7 = -1 * tmp6;
      assert(tmp7 == tmp7);
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp8 = jacobianCoefficient< 1 >( x );
      assert(tmp8 == tmp8);
    const auto tmp9 = (tmp8[ 0 ])[ 0 ] + tmp7;
      assert(tmp9 == tmp9);
    const auto tmp10 = 2 * tmp9;
      assert(tmp10 == tmp10);
    const auto tmp11 = std::sin( tmp3 );
      assert(tmp11 == tmp11);
    const auto tmp12 = tmp11 * tmp2;
      assert(tmp12 == tmp12);
    const auto tmp13 = -1 * tmp12;
      assert(tmp13 == tmp13);
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp14 = evaluateCoefficient< 1 >( x );
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14[ 0 ] + tmp13;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 * tmp10;
      assert(tmp16 == tmp16);
    const auto tmp17 = std::cos( tmp1 );
      assert(tmp17 == tmp17);
    const auto tmp18 = 3.141592653589793 * tmp17;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 * tmp11;
      assert(tmp19 == tmp19);
    const auto tmp20 = -1 * tmp19;
      assert(tmp20 == tmp20);
    const auto tmp21 = (tmp8[ 0 ])[ 1 ] + tmp20;
      assert(tmp21 == tmp21);
    const auto tmp22 = 2 * tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp15 * tmp22;
      assert(tmp23 == tmp23);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp24 = evaluateCoefficient< 0 >( x );
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp24[ 0 ] + tmp7;
      assert(tmp25 == tmp25);
    const auto tmp26 = -1 * tmp11;
      assert(tmp26 == tmp26);
    const auto tmp27 = 3.141592653589793 * tmp26;
      assert(tmp27 == tmp27);
    const auto tmp28 = 3.141592653589793 * tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp28 * tmp2;
      assert(tmp29 == tmp29);
    const auto tmp30 = -1 * tmp29;
      assert(tmp30 == tmp30);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp31 = jacobianCoefficient< 0 >( x );
      assert(tmp31 == tmp31);
    const auto tmp32 = (tmp31[ 0 ])[ 0 ] + tmp30;
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp32 * tmp25;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp25 * tmp32;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp34 + tmp33;
      assert(tmp35 == tmp35);
    const auto tmp36 = tmp24[ 1 ] + tmp20;
      assert(tmp36 == tmp36);
    const auto tmp37 = tmp5 * tmp18;
      assert(tmp37 == tmp37);
    const auto tmp38 = -1 * tmp37;
      assert(tmp38 == tmp38);
    const auto tmp39 = (tmp31[ 1 ])[ 0 ] + tmp38;
      assert(tmp39 == tmp39);
    const auto tmp40 = tmp39 * tmp36;
      assert(tmp40 == tmp40);
    const auto tmp41 = tmp36 * tmp39;
      assert(tmp41 == tmp41);
    const auto tmp42 = tmp41 + tmp40;
      assert(tmp42 == tmp42);
    const auto tmp43 = tmp42 + tmp35;
      assert(tmp43 == tmp43);
    const auto tmp44 = (tmp31[ 0 ])[ 1 ] + tmp38;
      assert(tmp44 == tmp44);
    const auto tmp45 = tmp44 * tmp25;
      assert(tmp45 == tmp45);
    const auto tmp46 = tmp25 * tmp44;
      assert(tmp46 == tmp46);
    const auto tmp47 = tmp46 + tmp45;
      assert(tmp47 == tmp47);
    const auto tmp48 = -1 * tmp2;
      assert(tmp48 == tmp48);
    const auto tmp49 = 3.141592653589793 * tmp48;
      assert(tmp49 == tmp49);
    const auto tmp50 = 3.141592653589793 * tmp49;
      assert(tmp50 == tmp50);
    const auto tmp51 = tmp50 * tmp11;
      assert(tmp51 == tmp51);
    const auto tmp52 = -1 * tmp51;
      assert(tmp52 == tmp52);
    const auto tmp53 = (tmp31[ 1 ])[ 1 ] + tmp52;
      assert(tmp53 == tmp53);
    const auto tmp54 = tmp53 * tmp36;
      assert(tmp54 == tmp54);
    const auto tmp55 = tmp36 * tmp53;
      assert(tmp55 == tmp55);
    const auto tmp56 = tmp55 + tmp54;
      assert(tmp56 == tmp56);
    const auto tmp57 = tmp56 + tmp47;
      assert(tmp57 == tmp57);
    const auto tmp58 = -1 * tmp5;
      assert(tmp58 == tmp58);
    const auto tmp59 = 3.141592653589793 * tmp58;
      assert(tmp59 == tmp59);
    const auto tmp60 = 3.141592653589793 * tmp59;
      assert(tmp60 == tmp60);
    const auto tmp61 = tmp60 * tmp2;
      assert(tmp61 == tmp61);
    const auto tmp62 = -1 * tmp61;
      assert(tmp62 == tmp62);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp63 = hessianCoefficient< 0 >( x );
      assert(tmp63 == tmp63);
    const auto tmp64 = ((tmp63[ 0 ])[ 0 ])[ 0 ] + tmp62;
      assert(tmp64 == tmp64);
    const auto tmp65 = tmp50 * tmp5;
      assert(tmp65 == tmp65);
    const auto tmp66 = -1 * tmp65;
      assert(tmp66 == tmp66);
    const auto tmp67 = ((tmp63[ 1 ])[ 1 ])[ 0 ] + tmp66;
      assert(tmp67 == tmp67);
    const auto tmp68 = tmp67 + tmp64;
      assert(tmp68 == tmp68);
    const auto tmp69 = 2 * tmp68;
      assert(tmp69 == tmp69);
    const auto tmp70 = tmp53 + tmp32;
      assert(tmp70 == tmp70);
    const auto tmp71 = tmp70 * tmp69;
      assert(tmp71 == tmp71);
    const auto tmp72 = -1 * tmp18;
      assert(tmp72 == tmp72);
    const auto tmp73 = 3.141592653589793 * tmp72;
      assert(tmp73 == tmp73);
    const auto tmp74 = 3.141592653589793 * tmp73;
      assert(tmp74 == tmp74);
    const auto tmp75 = tmp74 * tmp11;
      assert(tmp75 == tmp75);
    const auto tmp76 = -1 * tmp75;
      assert(tmp76 == tmp76);
    const auto tmp77 = ((tmp63[ 1 ])[ 1 ])[ 1 ] + tmp76;
      assert(tmp77 == tmp77);
    const auto tmp78 = tmp28 * tmp18;
      assert(tmp78 == tmp78);
    const auto tmp79 = -1 * tmp78;
      assert(tmp79 == tmp79);
    const auto tmp80 = ((tmp63[ 0 ])[ 0 ])[ 1 ] + tmp79;
      assert(tmp80 == tmp80);
    const auto tmp81 = tmp80 + tmp77;
      assert(tmp81 == tmp81);
    const auto tmp82 = 2 * tmp81;
      assert(tmp82 == tmp82);
    const auto tmp83 = tmp70 * tmp82;
      assert(tmp83 == tmp83);
    (result[ 0 ])[ 0 ] = tmp16;
    (result[ 0 ])[ 1 ] = tmp23;
    (result[ 1 ])[ 0 ] = tmp43;
    (result[ 1 ])[ 1 ] = tmp57;
    (result[ 2 ])[ 0 ] = tmp71;
    (result[ 2 ])[ 1 ] = tmp83;
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
  std::tuple< Dune::Fem::ConstLocalFunction< CoeffgradSol >, Dune::Fem::ConstLocalFunction< Coeffsol > > coefficients_;
};

} // namespace UFLLocalFunctions_ff18c3821ae99734bca838539098cb6c

PYBIND11_MODULE( localfunction_ff18c3821ae99734bca838539098cb6c_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef UFLLocalFunctions_ff18c3821ae99734bca838539098cb6c::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_ff18c3821ae99734bca838539098cb6c::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_ff18c3821ae99734bca838539098cb6c_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > > &coeffgradSol, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffsol ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>(gridView),name,order, coeffgradSol, coeffsol); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
