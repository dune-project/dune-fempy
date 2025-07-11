#ifndef GUARD_7d2cbec2fa1ee20560f6b1692c14ba44
#define GUARD_7d2cbec2fa1ee20560f6b1692c14ba44

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

namespace UFLLocalFunctions_7d2cbec2fa1ee20560f6b1692c14ba44
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class CoeffaTimesExact >
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
  typedef std::tuple< CoeffaTimesExact > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<CoeffaTimesExact>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const CoeffaTimesExact &coeffaTimesExact, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< CoeffaTimesExact >( coeffaTimesExact ) )
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
    const auto tmp1 = -1 * tmp0[ 0 ];
    GlobalCoordinateType tmp2 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp3 = std::pow( tmp2[ 0 ], 3 );
    const auto tmp4 = std::pow( tmp2[ 1 ], 3 );
    const auto tmp5 = -1 * tmp4;
    const auto tmp6 = tmp5 + tmp3;
    const auto tmp7 = 0.3333333333333333 * tmp6;
    const auto tmp8 = -1 * tmp7;
    const auto tmp9 = std::pow( tmp2[ 1 ], 2 );
    const auto tmp10 = std::pow( tmp2[ 0 ], 2 );
    const auto tmp11 = tmp10 + tmp9;
    const auto tmp12 = tmp11 + tmp8;
    const auto tmp13 = 1 + tmp12;
    const auto tmp14 = tmp13 + tmp1;
    const auto tmp15 = std::abs( tmp14 );
    result[ 0 ] = tmp15;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::pow;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
    const auto tmp1 = -1 * tmp0[ 0 ];
    GlobalCoordinateType tmp2 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp3 = std::pow( tmp2[ 0 ], 3 );
    const auto tmp4 = std::pow( tmp2[ 1 ], 3 );
    const auto tmp5 = -1 * tmp4;
    const auto tmp6 = tmp5 + tmp3;
    const auto tmp7 = 0.3333333333333333 * tmp6;
    const auto tmp8 = -1 * tmp7;
    const auto tmp9 = std::pow( tmp2[ 1 ], 2 );
    const auto tmp10 = std::pow( tmp2[ 0 ], 2 );
    const auto tmp11 = tmp10 + tmp9;
    const auto tmp12 = tmp11 + tmp8;
    const auto tmp13 = 1 + tmp12;
    const auto tmp14 = tmp13 + tmp1;
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp15 = jacobianCoefficient< 0 >( x );
    const auto tmp16 = -1 * (tmp15[ 0 ])[ 0 ];
    const auto tmp17 = 3 * tmp10;
    const auto tmp18 = 0.3333333333333333 * tmp17;
    const auto tmp19 = -1 * tmp18;
    const auto tmp20 = 2 * tmp2[ 0 ];
    const auto tmp21 = tmp20 + tmp19;
    const auto tmp22 = tmp21 + tmp16;
    const auto tmp23 = tmp22 * (tmp14 == 0.0 ? 0.0 : (tmp14 < 0.0 ? -1 : 1));
    const auto tmp24 = -1 * (tmp15[ 0 ])[ 1 ];
    const auto tmp25 = 3 * tmp9;
    const auto tmp26 = -1 * tmp25;
    const auto tmp27 = 0.3333333333333333 * tmp26;
    const auto tmp28 = -1 * tmp27;
    const auto tmp29 = 2 * tmp2[ 1 ];
    const auto tmp30 = tmp29 + tmp28;
    const auto tmp31 = tmp30 + tmp24;
    const auto tmp32 = tmp31 * (tmp14 == 0.0 ? 0.0 : (tmp14 < 0.0 ? -1 : 1));
    (result[ 0 ])[ 0 ] = tmp23;
    (result[ 0 ])[ 1 ] = tmp32;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::pow;
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
    const auto tmp1 = -1 * tmp0[ 0 ];
    GlobalCoordinateType tmp2 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp3 = std::pow( tmp2[ 0 ], 3 );
    const auto tmp4 = std::pow( tmp2[ 1 ], 3 );
    const auto tmp5 = -1 * tmp4;
    const auto tmp6 = tmp5 + tmp3;
    const auto tmp7 = 0.3333333333333333 * tmp6;
    const auto tmp8 = -1 * tmp7;
    const auto tmp9 = std::pow( tmp2[ 1 ], 2 );
    const auto tmp10 = std::pow( tmp2[ 0 ], 2 );
    const auto tmp11 = tmp10 + tmp9;
    const auto tmp12 = tmp11 + tmp8;
    const auto tmp13 = 1 + tmp12;
    const auto tmp14 = tmp13 + tmp1;
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp15 = hessianCoefficient< 0 >( x );
    const auto tmp16 = -1 * ((tmp15[ 0 ])[ 0 ])[ 0 ];
    const auto tmp17 = 2 * tmp2[ 0 ];
    const auto tmp18 = 3 * tmp17;
    const auto tmp19 = 0.3333333333333333 * tmp18;
    const auto tmp20 = -1 * tmp19;
    const auto tmp21 = 2 + tmp20;
    const auto tmp22 = tmp21 + tmp16;
    const auto tmp23 = tmp22 * (tmp14 == 0.0 ? 0.0 : (tmp14 < 0.0 ? -1 : 1));
    const auto tmp24 = -1 * ((tmp15[ 0 ])[ 0 ])[ 1 ];
    const auto tmp25 = tmp24 * (tmp14 == 0.0 ? 0.0 : (tmp14 < 0.0 ? -1 : 1));
    const auto tmp26 = -1 * ((tmp15[ 0 ])[ 1 ])[ 0 ];
    const auto tmp27 = tmp26 * (tmp14 == 0.0 ? 0.0 : (tmp14 < 0.0 ? -1 : 1));
    const auto tmp28 = -1 * ((tmp15[ 0 ])[ 1 ])[ 1 ];
    const auto tmp29 = 2 * tmp2[ 1 ];
    const auto tmp30 = 3 * tmp29;
    const auto tmp31 = -1 * tmp30;
    const auto tmp32 = 0.3333333333333333 * tmp31;
    const auto tmp33 = -1 * tmp32;
    const auto tmp34 = 2 + tmp33;
    const auto tmp35 = tmp34 + tmp28;
    const auto tmp36 = tmp35 * (tmp14 == 0.0 ? 0.0 : (tmp14 < 0.0 ? -1 : 1));
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp23;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp25;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp27;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp36;
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
  std::tuple< Dune::Fem::ConstLocalFunction< CoeffaTimesExact > > coefficients_;
};

} // namespace UFLLocalFunctions_7d2cbec2fa1ee20560f6b1692c14ba44

PYBIND11_MODULE( localfunction_7d2cbec2fa1ee20560f6b1692c14ba44_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_7d2cbec2fa1ee20560f6b1692c14ba44::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_7d2cbec2fa1ee20560f6b1692c14ba44::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_7d2cbec2fa1ee20560f6b1692c14ba44_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > &coeffaTimesExact ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffaTimesExact); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
