#ifndef GUARD_18e83e4b92a168a2aa116e94fc4541af
#define GUARD_18e83e4b92a168a2aa116e94fc4541af

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

namespace UFLLocalFunctions_18e83e4b92a168a2aa116e94fc4541af
{

  // UFLLocalFunction
// ----------------

template< class GridPart, class Coeffsolution_vector_cg, class Coeffanalytical_solution_x, class Coeffanalytical_solution_y, class Coeffanalytical_solution_pressure >
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
  typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 3 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
  typedef std::tuple< Coeffsolution_vector_cg, Coeffanalytical_solution_x, Coeffanalytical_solution_y, Coeffanalytical_solution_pressure > CoefficientTupleType;
  template< std::size_t i >
  using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
  template< std::size_t i >
  using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
  template< std::size_t i >
  using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
  static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffsolution_vector_cg>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_x>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_y>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_pressure>>();
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Coeffsolution_vector_cg &coeffsolution_vector_cg, const Coeffanalytical_solution_x &coeffanalytical_solution_x, const Coeffanalytical_solution_y &coeffanalytical_solution_y, const Coeffanalytical_solution_pressure &coeffanalytical_solution_pressure, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order),
      coefficients_( Dune::Fem::ConstLocalFunction< Coeffsolution_vector_cg >( coeffsolution_vector_cg ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >( coeffanalytical_solution_x ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y >( coeffanalytical_solution_y ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_pressure >( coeffanalytical_solution_pressure ) )
  {}

  void bind ( const IntersectionType &intersection, Side side )
  {
    BaseType::bind(intersection, side);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
    std::get< 1 >( coefficients_ ).bind( this->entity() );
    std::get< 2 >( coefficients_ ).bind( this->entity() );
    std::get< 3 >( coefficients_ ).bind( this->entity() );
  }

  void bind ( const EntityType &entity )
  {
    BaseType::bind(entity);
    std::get< 0 >( coefficients_ ).bind( this->entity() );
    std::get< 1 >( coefficients_ ).bind( this->entity() );
    std::get< 2 >( coefficients_ ).bind( this->entity() );
    std::get< 3 >( coefficients_ ).bind( this->entity() );
  }

  void unbind ()
  {
    BaseType::unbind();
    std::get< 0 >( coefficients_ ).unbind();
    std::get< 1 >( coefficients_ ).unbind();
    std::get< 2 >( coefficients_ ).unbind();
    std::get< 3 >( coefficients_ ).unbind();
  }

  template< class Point >
  void evaluate ( const Point &x, typename FunctionSpaceType::RangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 3 >::JacobianRangeType tmp0 = jacobianCoefficient< 3 >( x );
    const auto tmp1 = -1 * (tmp0[ 0 ])[ 1 ];
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp2 = jacobianCoefficient< 0 >( x );
    const auto tmp3 = (tmp2[ 2 ])[ 1 ] + tmp1;
    const auto tmp4 = tmp3 * tmp3;
    const auto tmp5 = -1 * (tmp0[ 0 ])[ 0 ];
    const auto tmp6 = (tmp2[ 2 ])[ 0 ] + tmp5;
    const auto tmp7 = tmp6 * tmp6;
    const auto tmp8 = tmp7 + tmp4;
    typename CoefficientFunctionSpaceType< 2 >::JacobianRangeType tmp9 = jacobianCoefficient< 2 >( x );
    const auto tmp10 = -1 * (tmp9[ 0 ])[ 1 ];
    const auto tmp11 = (tmp2[ 1 ])[ 1 ] + tmp10;
    const auto tmp12 = tmp11 * tmp11;
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp13 = jacobianCoefficient< 1 >( x );
    const auto tmp14 = -1 * (tmp13[ 0 ])[ 1 ];
    const auto tmp15 = (tmp2[ 0 ])[ 1 ] + tmp14;
    const auto tmp16 = tmp15 * tmp15;
    const auto tmp17 = tmp16 + tmp12;
    const auto tmp18 = -1 * (tmp9[ 0 ])[ 0 ];
    const auto tmp19 = (tmp2[ 1 ])[ 0 ] + tmp18;
    const auto tmp20 = tmp19 * tmp19;
    const auto tmp21 = -1 * (tmp13[ 0 ])[ 0 ];
    const auto tmp22 = (tmp2[ 0 ])[ 0 ] + tmp21;
    const auto tmp23 = tmp22 * tmp22;
    const auto tmp24 = tmp23 + tmp20;
    const auto tmp25 = tmp24 + tmp17;
    result[ 0 ] = tmp8;
    result[ 1 ] = tmp25;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    typename CoefficientFunctionSpaceType< 3 >::HessianRangeType tmp0 = hessianCoefficient< 3 >( x );
    const auto tmp1 = -1 * ((tmp0[ 0 ])[ 1 ])[ 0 ];
    typename CoefficientFunctionSpaceType< 0 >::HessianRangeType tmp2 = hessianCoefficient< 0 >( x );
    const auto tmp3 = ((tmp2[ 2 ])[ 1 ])[ 0 ] + tmp1;
    typename CoefficientFunctionSpaceType< 3 >::JacobianRangeType tmp4 = jacobianCoefficient< 3 >( x );
    const auto tmp5 = -1 * (tmp4[ 0 ])[ 1 ];
    typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp6 = jacobianCoefficient< 0 >( x );
    const auto tmp7 = (tmp6[ 2 ])[ 1 ] + tmp5;
    const auto tmp8 = tmp7 * tmp3;
    const auto tmp9 = tmp8 + tmp8;
    const auto tmp10 = -1 * ((tmp0[ 0 ])[ 0 ])[ 0 ];
    const auto tmp11 = ((tmp2[ 2 ])[ 0 ])[ 0 ] + tmp10;
    const auto tmp12 = -1 * (tmp4[ 0 ])[ 0 ];
    const auto tmp13 = (tmp6[ 2 ])[ 0 ] + tmp12;
    const auto tmp14 = tmp13 * tmp11;
    const auto tmp15 = tmp14 + tmp14;
    const auto tmp16 = tmp15 + tmp9;
    const auto tmp17 = -1 * ((tmp0[ 0 ])[ 1 ])[ 1 ];
    const auto tmp18 = ((tmp2[ 2 ])[ 1 ])[ 1 ] + tmp17;
    const auto tmp19 = tmp7 * tmp18;
    const auto tmp20 = tmp19 + tmp19;
    const auto tmp21 = -1 * ((tmp0[ 0 ])[ 0 ])[ 1 ];
    const auto tmp22 = ((tmp2[ 2 ])[ 0 ])[ 1 ] + tmp21;
    const auto tmp23 = tmp13 * tmp22;
    const auto tmp24 = tmp23 + tmp23;
    const auto tmp25 = tmp24 + tmp20;
    typename CoefficientFunctionSpaceType< 2 >::HessianRangeType tmp26 = hessianCoefficient< 2 >( x );
    const auto tmp27 = -1 * ((tmp26[ 0 ])[ 1 ])[ 0 ];
    const auto tmp28 = ((tmp2[ 1 ])[ 1 ])[ 0 ] + tmp27;
    typename CoefficientFunctionSpaceType< 2 >::JacobianRangeType tmp29 = jacobianCoefficient< 2 >( x );
    const auto tmp30 = -1 * (tmp29[ 0 ])[ 1 ];
    const auto tmp31 = (tmp6[ 1 ])[ 1 ] + tmp30;
    const auto tmp32 = tmp31 * tmp28;
    const auto tmp33 = tmp28 * tmp31;
    const auto tmp34 = tmp33 + tmp32;
    typename CoefficientFunctionSpaceType< 1 >::HessianRangeType tmp35 = hessianCoefficient< 1 >( x );
    const auto tmp36 = -1 * ((tmp35[ 0 ])[ 1 ])[ 0 ];
    const auto tmp37 = ((tmp2[ 0 ])[ 1 ])[ 0 ] + tmp36;
    typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp38 = jacobianCoefficient< 1 >( x );
    const auto tmp39 = -1 * (tmp38[ 0 ])[ 1 ];
    const auto tmp40 = (tmp6[ 0 ])[ 1 ] + tmp39;
    const auto tmp41 = tmp40 * tmp37;
    const auto tmp42 = tmp37 * tmp40;
    const auto tmp43 = tmp42 + tmp41;
    const auto tmp44 = tmp43 + tmp34;
    const auto tmp45 = -1 * ((tmp26[ 0 ])[ 0 ])[ 0 ];
    const auto tmp46 = ((tmp2[ 1 ])[ 0 ])[ 0 ] + tmp45;
    const auto tmp47 = -1 * (tmp29[ 0 ])[ 0 ];
    const auto tmp48 = (tmp6[ 1 ])[ 0 ] + tmp47;
    const auto tmp49 = tmp48 * tmp46;
    const auto tmp50 = tmp46 * tmp48;
    const auto tmp51 = tmp50 + tmp49;
    const auto tmp52 = -1 * ((tmp35[ 0 ])[ 0 ])[ 0 ];
    const auto tmp53 = ((tmp2[ 0 ])[ 0 ])[ 0 ] + tmp52;
    const auto tmp54 = -1 * (tmp38[ 0 ])[ 0 ];
    const auto tmp55 = (tmp6[ 0 ])[ 0 ] + tmp54;
    const auto tmp56 = tmp55 * tmp53;
    const auto tmp57 = tmp53 * tmp55;
    const auto tmp58 = tmp57 + tmp56;
    const auto tmp59 = tmp58 + tmp51;
    const auto tmp60 = tmp59 + tmp44;
    const auto tmp61 = -1 * ((tmp26[ 0 ])[ 1 ])[ 1 ];
    const auto tmp62 = ((tmp2[ 1 ])[ 1 ])[ 1 ] + tmp61;
    const auto tmp63 = tmp31 * tmp62;
    const auto tmp64 = tmp62 * tmp31;
    const auto tmp65 = tmp64 + tmp63;
    const auto tmp66 = -1 * ((tmp35[ 0 ])[ 1 ])[ 1 ];
    const auto tmp67 = ((tmp2[ 0 ])[ 1 ])[ 1 ] + tmp66;
    const auto tmp68 = tmp40 * tmp67;
    const auto tmp69 = tmp67 * tmp40;
    const auto tmp70 = tmp69 + tmp68;
    const auto tmp71 = tmp70 + tmp65;
    const auto tmp72 = -1 * ((tmp26[ 0 ])[ 0 ])[ 1 ];
    const auto tmp73 = ((tmp2[ 1 ])[ 0 ])[ 1 ] + tmp72;
    const auto tmp74 = tmp48 * tmp73;
    const auto tmp75 = tmp73 * tmp48;
    const auto tmp76 = tmp75 + tmp74;
    const auto tmp77 = -1 * ((tmp35[ 0 ])[ 0 ])[ 1 ];
    const auto tmp78 = ((tmp2[ 0 ])[ 0 ])[ 1 ] + tmp77;
    const auto tmp79 = tmp55 * tmp78;
    const auto tmp80 = tmp78 * tmp55;
    const auto tmp81 = tmp80 + tmp79;
    const auto tmp82 = tmp81 + tmp76;
    const auto tmp83 = tmp82 + tmp71;
    (result[ 0 ])[ 0 ] = tmp16;
    (result[ 0 ])[ 1 ] = tmp25;
    (result[ 1 ])[ 0 ] = tmp60;
    (result[ 1 ])[ 1 ] = tmp83;
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
  std::tuple< Dune::Fem::ConstLocalFunction< Coeffsolution_vector_cg >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_pressure > > coefficients_;
};

} // namespace UFLLocalFunctions_18e83e4b92a168a2aa116e94fc4541af

PYBIND11_MODULE( localfunction_18e83e4b92a168a2aa116e94fc4541af_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_18e83e4b92a168a2aa116e94fc4541af::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 3 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_18e83e4b92a168a2aa116e94fc4541af::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 3 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_18e83e4b92a168a2aa116e94fc4541af_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 3 > > &coeffsolution_vector_cg, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_x, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_y, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_pressure ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order, coeffsolution_vector_cg, coeffanalytical_solution_x, coeffanalytical_solution_y, coeffanalytical_solution_pressure); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >(), pybind11::keep_alive< 1, 5 >(), pybind11::keep_alive< 1, 6 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
