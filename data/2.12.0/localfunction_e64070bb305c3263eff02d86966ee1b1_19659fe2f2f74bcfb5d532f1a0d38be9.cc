#ifndef GUARD_e64070bb305c3263eff02d86966ee1b1
#define GUARD_e64070bb305c3263eff02d86966ee1b1

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

namespace UFLLocalFunctions_e64070bb305c3263eff02d86966ee1b1
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class CoeffaTimesExact, class CoeffaTimesExactA >
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
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< CoeffaTimesExact, CoeffaTimesExactA > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<CoeffaTimesExact>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<CoeffaTimesExactA>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const CoeffaTimesExact &coeffaTimesExact, const CoeffaTimesExactA &coeffaTimesExactA, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< CoeffaTimesExact >( coeffaTimesExact ), Dune::Fem::ConstLocalFunction< CoeffaTimesExactA >( coeffaTimesExactA ) )
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
    using std::abs;
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp0 = evaluateCoefficient< 1 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 * tmp0[ 0 ];
      assert(tmp1 == tmp1);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp2 = evaluateCoefficient< 0 >( x );
      assert(tmp2 == tmp2);
    const auto tmp3 = tmp2[ 0 ] + tmp1;
      assert(tmp3 == tmp3);
    const auto tmp4 = std::abs( tmp3 );
      assert(tmp4 == tmp4);
    result[ 0 ] = tmp4;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp0 = evaluateCoefficient< 1 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 * tmp0[ 0 ];
      assert(tmp1 == tmp1);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp2 = evaluateCoefficient< 0 >( x );
      assert(tmp2 == tmp2);
    const auto tmp3 = tmp2[ 0 ] + tmp1;
      assert(tmp3 == tmp3);
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp4 = jacobianCoefficient< 1 >( x );
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * (tmp4[ 0 ])[ 0 ];
      assert(tmp5 == tmp5);
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp6 = jacobianCoefficient< 0 >( x );
      assert(tmp6 == tmp6);
    const auto tmp7 = (tmp6[ 0 ])[ 0 ] + tmp5;
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 * (tmp3 == 0.0 ? 0.0 : (tmp3 < 0.0 ? -1 : 1));
      assert(tmp8 == tmp8);
    const auto tmp9 = -1 * (tmp4[ 0 ])[ 1 ];
      assert(tmp9 == tmp9);
    const auto tmp10 = (tmp6[ 0 ])[ 1 ] + tmp9;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * (tmp3 == 0.0 ? 0.0 : (tmp3 < 0.0 ? -1 : 1));
      assert(tmp11 == tmp11);
    (result[ 0 ])[ 0 ] = tmp8;
    (result[ 0 ])[ 1 ] = tmp11;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 1 >::RangeType tmp0 = evaluateCoefficient< 1 >( x );
      assert(tmp0 == tmp0);
    const auto tmp1 = -1 * tmp0[ 0 ];
      assert(tmp1 == tmp1);
    typename CoefficientFunctionSpaceType< 0 >::RangeType tmp2 = evaluateCoefficient< 0 >( x );
      assert(tmp2 == tmp2);
    const auto tmp3 = tmp2[ 0 ] + tmp1;
      assert(tmp3 == tmp3);
    typename CoefficientFunctionSpaceType< 1 >::HessianRangeType tmp4 = hessianCoefficient< 1 >( x );
      assert(tmp4 == tmp4);
    const auto tmp5 = -1 * ((tmp4[ 0 ])[ 0 ])[ 0 ];
      assert(tmp5 == tmp5);
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp6 = hessianCoefficient< 0 >( x );
      assert(tmp6 == tmp6);
    const auto tmp7 = ((tmp6[ 0 ])[ 0 ])[ 0 ] + tmp5;
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 * (tmp3 == 0.0 ? 0.0 : (tmp3 < 0.0 ? -1 : 1));
      assert(tmp8 == tmp8);
    const auto tmp9 = -1 * ((tmp4[ 0 ])[ 0 ])[ 1 ];
      assert(tmp9 == tmp9);
    const auto tmp10 = ((tmp6[ 0 ])[ 0 ])[ 1 ] + tmp9;
      assert(tmp10 == tmp10);
    const auto tmp11 = tmp10 * (tmp3 == 0.0 ? 0.0 : (tmp3 < 0.0 ? -1 : 1));
      assert(tmp11 == tmp11);
    const auto tmp12 = -1 * ((tmp4[ 0 ])[ 1 ])[ 0 ];
      assert(tmp12 == tmp12);
    const auto tmp13 = ((tmp6[ 0 ])[ 1 ])[ 0 ] + tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = tmp13 * (tmp3 == 0.0 ? 0.0 : (tmp3 < 0.0 ? -1 : 1));
      assert(tmp14 == tmp14);
    const auto tmp15 = -1 * ((tmp4[ 0 ])[ 1 ])[ 1 ];
      assert(tmp15 == tmp15);
    const auto tmp16 = ((tmp6[ 0 ])[ 1 ])[ 1 ] + tmp15;
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp16 * (tmp3 == 0.0 ? 0.0 : (tmp3 < 0.0 ? -1 : 1));
      assert(tmp17 == tmp17);
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp8;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp11;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp14;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp17;
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
  std::tuple< Dune::Fem::ConstLocalFunction< CoeffaTimesExact >, Dune::Fem::ConstLocalFunction< CoeffaTimesExactA > > coefficients_;
};

} // namespace UFLLocalFunctions_e64070bb305c3263eff02d86966ee1b1

PYBIND11_MODULE( localfunction_e64070bb305c3263eff02d86966ee1b1_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_e64070bb305c3263eff02d86966ee1b1::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > >, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_e64070bb305c3263eff02d86966ee1b1::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > >, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_e64070bb305c3263eff02d86966ee1b1_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > &coeffaTimesExact, const Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > &coeffaTimesExactA ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffaTimesExact, coeffaTimesExactA); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
