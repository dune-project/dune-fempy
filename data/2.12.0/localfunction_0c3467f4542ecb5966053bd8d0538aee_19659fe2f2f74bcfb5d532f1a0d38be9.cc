#ifndef GUARD_0c3467f4542ecb5966053bd8d0538aee
#define GUARD_0c3467f4542ecb5966053bd8d0538aee

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/common/exceptions.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_0c3467f4542ecb5966053bd8d0538aee
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffcallback >
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
  typedef std::tuple< Coeffcallback > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffcallback>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffcallback &coeffcallback, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffcallback >( coeffcallback ) )
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
    using std::abs;
    using std::pow;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 * tmp0[ 0 ];
      assert(tmp1 == tmp1);
    GlobalCoordinateType tmp2 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp2 == tmp2);
    const auto tmp3 = std::pow( tmp2[ 0 ], 3 );
      assert(tmp3 == tmp3);
    const auto tmp4 = std::pow( tmp2[ 1 ], 3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp5 + tmp3;
      assert(tmp6 == tmp6);
    const auto tmp7 = 0.3333333333333333 * tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = -1 * tmp7;
      assert(tmp8 == tmp8);
    const auto tmp9 = std::pow( tmp2[ 1 ], 2 );
      assert(tmp9 == tmp9);
    const auto tmp10 = std::pow( tmp2[ 0 ], 2 );
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 + tmp9;
      assert(tmp11 == tmp11);
    const auto tmp12 = 0.5 * tmp11;
      assert(tmp12 == tmp12);
    const auto tmp13 = tmp12 + tmp8;
      assert(tmp13 == tmp13);
    const auto tmp14 = 1 + tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 + tmp1;
      assert(tmp15 == tmp15);
    const auto tmp16 = std::abs( tmp15 );
      assert(tmp16 == tmp16);
    result[ 0 ] = tmp16;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::pow;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 * tmp0[ 0 ];
      assert(tmp1 == tmp1);
    GlobalCoordinateType tmp2 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp2 == tmp2);
    const auto tmp3 = std::pow( tmp2[ 0 ], 3 );
      assert(tmp3 == tmp3);
    const auto tmp4 = std::pow( tmp2[ 1 ], 3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp5 + tmp3;
      assert(tmp6 == tmp6);
    const auto tmp7 = 0.3333333333333333 * tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = -1 * tmp7;
      assert(tmp8 == tmp8);
    const auto tmp9 = std::pow( tmp2[ 1 ], 2 );
      assert(tmp9 == tmp9);
    const auto tmp10 = std::pow( tmp2[ 0 ], 2 );
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 + tmp9;
      assert(tmp11 == tmp11);
    const auto tmp12 = 0.5 * tmp11;
      assert(tmp12 == tmp12);
    const auto tmp13 = tmp12 + tmp8;
      assert(tmp13 == tmp13);
    const auto tmp14 = 1 + tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 + tmp1;
      assert(tmp15 == tmp15);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp16 = jacobianCoefficient< 0 >( x );
      assert(tmp16 == tmp16);
    const auto tmp17 = -1 * (tmp16[ 0 ])[ 0 ];
      assert(tmp17 == tmp17);
    const auto tmp18 = 3 * tmp10;
      assert(tmp18 == tmp18);
    const auto tmp19 = 0.3333333333333333 * tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = -1 * tmp19;
      assert(tmp20 == tmp20);
    const auto tmp21 = 2 * tmp2[ 0 ];
      assert(tmp21 == tmp21);
    const auto tmp22 = 0.5 * tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 + tmp20;
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp23 + tmp17;
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp24 * (tmp15 == 0.0 ? 0.0 : (tmp15 < 0.0 ? -1 : 1));
      assert(tmp25 == tmp25);
    const auto tmp26 = -1 * (tmp16[ 0 ])[ 1 ];
      assert(tmp26 == tmp26);
    const auto tmp27 = 3 * tmp9;
      assert(tmp27 == tmp27);
    const auto tmp28 = -1 * tmp27;
      assert(tmp28 == tmp28);
    const auto tmp29 = 0.3333333333333333 * tmp28;
      assert(tmp29 == tmp29);
    const auto tmp30 = -1 * tmp29;
      assert(tmp30 == tmp30);
    const auto tmp31 = 2 * tmp2[ 1 ];
      assert(tmp31 == tmp31);
    const auto tmp32 = 0.5 * tmp31;
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp32 + tmp30;
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp33 + tmp26;
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp34 * (tmp15 == 0.0 ? 0.0 : (tmp15 < 0.0 ? -1 : 1));
      assert(tmp35 == tmp35);
    (result[ 0 ])[ 0 ] = tmp25;
    (result[ 0 ])[ 1 ] = tmp35;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::pow;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 * tmp0[ 0 ];
      assert(tmp1 == tmp1);
    GlobalCoordinateType tmp2 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp2 == tmp2);
    const auto tmp3 = std::pow( tmp2[ 0 ], 3 );
      assert(tmp3 == tmp3);
    const auto tmp4 = std::pow( tmp2[ 1 ], 3 );
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = tmp5 + tmp3;
      assert(tmp6 == tmp6);
    const auto tmp7 = 0.3333333333333333 * tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = -1 * tmp7;
      assert(tmp8 == tmp8);
    const auto tmp9 = std::pow( tmp2[ 1 ], 2 );
      assert(tmp9 == tmp9);
    const auto tmp10 = std::pow( tmp2[ 0 ], 2 );
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 + tmp9;
      assert(tmp11 == tmp11);
    const auto tmp12 = 0.5 * tmp11;
      assert(tmp12 == tmp12);
    const auto tmp13 = tmp12 + tmp8;
      assert(tmp13 == tmp13);
    const auto tmp14 = 1 + tmp13;
      assert(tmp14 == tmp14);
    const auto tmp15 = tmp14 + tmp1;
      assert(tmp15 == tmp15);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp16 = hessianCoefficient< 0 >( x );
      assert(tmp16 == tmp16);
    const auto tmp17 = -1 * ((tmp16[ 0 ])[ 0 ])[ 0 ];
      assert(tmp17 == tmp17);
    const auto tmp18 = 2 * tmp2[ 0 ];
      assert(tmp18 == tmp18);
    const auto tmp19 = 3 * tmp18;
      assert(tmp19 == tmp19);
    const auto tmp20 = 0.3333333333333333 * tmp19;
      assert(tmp20 == tmp20);
    const auto tmp21 = -1 * tmp20;
      assert(tmp21 == tmp21);
    const auto tmp22 = 1.0 + tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp22 + tmp17;
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp23 * (tmp15 == 0.0 ? 0.0 : (tmp15 < 0.0 ? -1 : 1));
      assert(tmp24 == tmp24);
    const auto tmp25 = -1 * ((tmp16[ 0 ])[ 0 ])[ 1 ];
      assert(tmp25 == tmp25);
    const auto tmp26 = tmp25 * (tmp15 == 0.0 ? 0.0 : (tmp15 < 0.0 ? -1 : 1));
      assert(tmp26 == tmp26);
    const auto tmp27 = -1 * ((tmp16[ 0 ])[ 1 ])[ 0 ];
      assert(tmp27 == tmp27);
    const auto tmp28 = tmp27 * (tmp15 == 0.0 ? 0.0 : (tmp15 < 0.0 ? -1 : 1));
      assert(tmp28 == tmp28);
    const auto tmp29 = -1 * ((tmp16[ 0 ])[ 1 ])[ 1 ];
      assert(tmp29 == tmp29);
    const auto tmp30 = 2 * tmp2[ 1 ];
      assert(tmp30 == tmp30);
    const auto tmp31 = 3 * tmp30;
      assert(tmp31 == tmp31);
    const auto tmp32 = -1 * tmp31;
      assert(tmp32 == tmp32);
    const auto tmp33 = 0.3333333333333333 * tmp32;
      assert(tmp33 == tmp33);
    const auto tmp34 = -1 * tmp33;
      assert(tmp34 == tmp34);
    const auto tmp35 = 1.0 + tmp34;
      assert(tmp35 == tmp35);
    const auto tmp36 = tmp35 + tmp29;
      assert(tmp36 == tmp36);
    const auto tmp37 = tmp36 * (tmp15 == 0.0 ? 0.0 : (tmp15 < 0.0 ? -1 : 1));
      assert(tmp37 == tmp37);
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp24;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp26;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp28;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp37;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffcallback > > coefficients_;
};

} // namespace UFLLocalFunctions_0c3467f4542ecb5966053bd8d0538aee

PYBIND11_MODULE( localfunction_0c3467f4542ecb5966053bd8d0538aee_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_0c3467f4542ecb5966053bd8d0538aee::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, pybind11::function > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_0c3467f4542ecb5966053bd8d0538aee::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, pybind11::function > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_0c3467f4542ecb5966053bd8d0538aee_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, pybind11::function > > &coeffcallback ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffcallback); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
