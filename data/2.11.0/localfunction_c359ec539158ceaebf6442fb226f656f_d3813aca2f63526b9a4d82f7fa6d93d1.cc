#ifndef GUARD_c359ec539158ceaebf6442fb226f656f
#define GUARD_c359ec539158ceaebf6442fb226f656f

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
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_c359ec539158ceaebf6442fb226f656f
{

  // UFLLocalFunction
// ----------------

template< class GridPart >
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
  typedef std::tuple<> CoefficientTupleType;
  static constexpr bool gridPartValid = true;
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order)
  {}

  void bind ( const IntersectionType &intersection, Side side )
  {
    BaseType::bind(intersection, side);
  }

  void bind ( const EntityType &entity )
  {
    BaseType::bind(entity);
  }

  void unbind ()
  {
    BaseType::unbind();
  }

  template< class Point >
  void evaluate ( const Point &x, typename FunctionSpaceType::RangeType &result ) const
  {
    using std::sqrt;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = tmp0[ 2 ] * tmp0[ 2 ];
    const auto tmp2 = tmp0[ 1 ] * tmp0[ 1 ];
    const auto tmp3 = tmp0[ 0 ] * tmp0[ 0 ];
    const auto tmp4 = tmp3 + tmp2;
    const auto tmp5 = tmp4 + tmp1;
    const auto tmp6 = std::sqrt( tmp5 );
    const auto tmp7 = 2.0 / tmp6;
    const auto tmp8 = tmp0[ 0 ] * tmp7;
    const auto tmp9 = tmp0[ 1 ] * tmp7;
    const auto tmp10 = tmp0[ 2 ] * tmp7;
    result[ 0 ] = tmp8;
    result[ 1 ] = tmp9;
    result[ 2 ] = tmp10;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::sqrt;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = tmp0[ 2 ] * tmp0[ 2 ];
    const auto tmp2 = tmp0[ 1 ] * tmp0[ 1 ];
    const auto tmp3 = tmp0[ 0 ] * tmp0[ 0 ];
    const auto tmp4 = tmp3 + tmp2;
    const auto tmp5 = tmp4 + tmp1;
    const auto tmp6 = std::sqrt( tmp5 );
    const auto tmp7 = 2.0 / tmp6;
    const auto tmp8 = 2 * tmp6;
    const auto tmp9 = tmp0[ 0 ] + tmp0[ 0 ];
    const auto tmp10 = tmp9 / tmp8;
    const auto tmp11 = tmp10 * tmp7;
    const auto tmp12 = -1 * tmp11;
    const auto tmp13 = tmp12 / tmp6;
    const auto tmp14 = tmp0[ 0 ] * tmp13;
    const auto tmp15 = tmp14 + tmp7;
    const auto tmp16 = tmp0[ 1 ] + tmp0[ 1 ];
    const auto tmp17 = tmp16 / tmp8;
    const auto tmp18 = tmp17 * tmp7;
    const auto tmp19 = -1 * tmp18;
    const auto tmp20 = tmp19 / tmp6;
    const auto tmp21 = tmp0[ 0 ] * tmp20;
    const auto tmp22 = tmp0[ 2 ] + tmp0[ 2 ];
    const auto tmp23 = tmp22 / tmp8;
    const auto tmp24 = tmp23 * tmp7;
    const auto tmp25 = -1 * tmp24;
    const auto tmp26 = tmp25 / tmp6;
    const auto tmp27 = tmp0[ 0 ] * tmp26;
    const auto tmp28 = tmp0[ 1 ] * tmp13;
    const auto tmp29 = tmp0[ 1 ] * tmp20;
    const auto tmp30 = tmp29 + tmp7;
    const auto tmp31 = tmp0[ 1 ] * tmp26;
    const auto tmp32 = tmp0[ 2 ] * tmp13;
    const auto tmp33 = tmp0[ 2 ] * tmp20;
    const auto tmp34 = tmp0[ 2 ] * tmp26;
    const auto tmp35 = tmp34 + tmp7;
    (result[ 0 ])[ 0 ] = tmp15;
    (result[ 0 ])[ 1 ] = tmp21;
    (result[ 0 ])[ 2 ] = tmp27;
    (result[ 1 ])[ 0 ] = tmp28;
    (result[ 1 ])[ 1 ] = tmp30;
    (result[ 1 ])[ 2 ] = tmp31;
    (result[ 2 ])[ 0 ] = tmp32;
    (result[ 2 ])[ 1 ] = tmp33;
    (result[ 2 ])[ 2 ] = tmp35;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::sqrt;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = tmp0[ 2 ] * tmp0[ 2 ];
    const auto tmp2 = tmp0[ 1 ] * tmp0[ 1 ];
    const auto tmp3 = tmp0[ 0 ] * tmp0[ 0 ];
    const auto tmp4 = tmp3 + tmp2;
    const auto tmp5 = tmp4 + tmp1;
    const auto tmp6 = std::sqrt( tmp5 );
    const auto tmp7 = 2.0 / tmp6;
    const auto tmp8 = 2 * tmp6;
    const auto tmp9 = tmp0[ 0 ] + tmp0[ 0 ];
    const auto tmp10 = tmp9 / tmp8;
    const auto tmp11 = tmp10 * tmp7;
    const auto tmp12 = -1 * tmp11;
    const auto tmp13 = tmp12 / tmp6;
    const auto tmp14 = tmp10 * tmp13;
    const auto tmp15 = -1 * tmp14;
    const auto tmp16 = 2 * tmp10;
    const auto tmp17 = tmp16 * tmp10;
    const auto tmp18 = -1 * tmp17;
    const auto tmp19 = 2 + tmp18;
    const auto tmp20 = tmp19 / tmp8;
    const auto tmp21 = tmp20 * tmp7;
    const auto tmp22 = tmp21 + tmp14;
    const auto tmp23 = -1 * tmp22;
    const auto tmp24 = tmp23 + tmp15;
    const auto tmp25 = tmp24 / tmp6;
    const auto tmp26 = tmp0[ 0 ] * tmp25;
    const auto tmp27 = tmp26 + tmp13;
    const auto tmp28 = tmp27 + tmp13;
    const auto tmp29 = tmp0[ 1 ] + tmp0[ 1 ];
    const auto tmp30 = tmp29 / tmp8;
    const auto tmp31 = tmp30 * tmp7;
    const auto tmp32 = -1 * tmp31;
    const auto tmp33 = tmp32 / tmp6;
    const auto tmp34 = tmp30 * tmp13;
    const auto tmp35 = -1 * tmp34;
    const auto tmp36 = tmp10 * tmp33;
    const auto tmp37 = 2 * tmp30;
    const auto tmp38 = tmp37 * tmp10;
    const auto tmp39 = -1 * tmp38;
    const auto tmp40 = tmp39 / tmp8;
    const auto tmp41 = tmp40 * tmp7;
    const auto tmp42 = tmp41 + tmp36;
    const auto tmp43 = -1 * tmp42;
    const auto tmp44 = tmp43 + tmp35;
    const auto tmp45 = tmp44 / tmp6;
    const auto tmp46 = tmp0[ 0 ] * tmp45;
    const auto tmp47 = tmp46 + tmp33;
    const auto tmp48 = tmp0[ 2 ] + tmp0[ 2 ];
    const auto tmp49 = tmp48 / tmp8;
    const auto tmp50 = tmp49 * tmp7;
    const auto tmp51 = -1 * tmp50;
    const auto tmp52 = tmp51 / tmp6;
    const auto tmp53 = tmp49 * tmp13;
    const auto tmp54 = -1 * tmp53;
    const auto tmp55 = tmp10 * tmp52;
    const auto tmp56 = 2 * tmp49;
    const auto tmp57 = tmp56 * tmp10;
    const auto tmp58 = -1 * tmp57;
    const auto tmp59 = tmp58 / tmp8;
    const auto tmp60 = tmp59 * tmp7;
    const auto tmp61 = tmp60 + tmp55;
    const auto tmp62 = -1 * tmp61;
    const auto tmp63 = tmp62 + tmp54;
    const auto tmp64 = tmp63 / tmp6;
    const auto tmp65 = tmp0[ 0 ] * tmp64;
    const auto tmp66 = tmp65 + tmp52;
    const auto tmp67 = -1 * tmp36;
    const auto tmp68 = tmp16 * tmp30;
    const auto tmp69 = -1 * tmp68;
    const auto tmp70 = tmp69 / tmp8;
    const auto tmp71 = tmp70 * tmp7;
    const auto tmp72 = tmp71 + tmp34;
    const auto tmp73 = -1 * tmp72;
    const auto tmp74 = tmp73 + tmp67;
    const auto tmp75 = tmp74 / tmp6;
    const auto tmp76 = tmp0[ 0 ] * tmp75;
    const auto tmp77 = tmp76 + tmp33;
    const auto tmp78 = tmp30 * tmp33;
    const auto tmp79 = -1 * tmp78;
    const auto tmp80 = tmp37 * tmp30;
    const auto tmp81 = -1 * tmp80;
    const auto tmp82 = 2 + tmp81;
    const auto tmp83 = tmp82 / tmp8;
    const auto tmp84 = tmp83 * tmp7;
    const auto tmp85 = tmp84 + tmp78;
    const auto tmp86 = -1 * tmp85;
    const auto tmp87 = tmp86 + tmp79;
    const auto tmp88 = tmp87 / tmp6;
    const auto tmp89 = tmp0[ 0 ] * tmp88;
    const auto tmp90 = tmp49 * tmp33;
    const auto tmp91 = -1 * tmp90;
    const auto tmp92 = tmp30 * tmp52;
    const auto tmp93 = tmp56 * tmp30;
    const auto tmp94 = -1 * tmp93;
    const auto tmp95 = tmp94 / tmp8;
    const auto tmp96 = tmp95 * tmp7;
    const auto tmp97 = tmp96 + tmp92;
    const auto tmp98 = -1 * tmp97;
    const auto tmp99 = tmp98 + tmp91;
    const auto tmp100 = tmp99 / tmp6;
    const auto tmp101 = tmp0[ 0 ] * tmp100;
    const auto tmp102 = -1 * tmp55;
    const auto tmp103 = tmp16 * tmp49;
    const auto tmp104 = -1 * tmp103;
    const auto tmp105 = tmp104 / tmp8;
    const auto tmp106 = tmp105 * tmp7;
    const auto tmp107 = tmp106 + tmp53;
    const auto tmp108 = -1 * tmp107;
    const auto tmp109 = tmp108 + tmp102;
    const auto tmp110 = tmp109 / tmp6;
    const auto tmp111 = tmp0[ 0 ] * tmp110;
    const auto tmp112 = tmp111 + tmp52;
    const auto tmp113 = -1 * tmp92;
    const auto tmp114 = tmp37 * tmp49;
    const auto tmp115 = -1 * tmp114;
    const auto tmp116 = tmp115 / tmp8;
    const auto tmp117 = tmp116 * tmp7;
    const auto tmp118 = tmp117 + tmp90;
    const auto tmp119 = -1 * tmp118;
    const auto tmp120 = tmp119 + tmp113;
    const auto tmp121 = tmp120 / tmp6;
    const auto tmp122 = tmp0[ 0 ] * tmp121;
    const auto tmp123 = tmp49 * tmp52;
    const auto tmp124 = -1 * tmp123;
    const auto tmp125 = tmp56 * tmp49;
    const auto tmp126 = -1 * tmp125;
    const auto tmp127 = 2 + tmp126;
    const auto tmp128 = tmp127 / tmp8;
    const auto tmp129 = tmp128 * tmp7;
    const auto tmp130 = tmp129 + tmp123;
    const auto tmp131 = -1 * tmp130;
    const auto tmp132 = tmp131 + tmp124;
    const auto tmp133 = tmp132 / tmp6;
    const auto tmp134 = tmp0[ 0 ] * tmp133;
    const auto tmp135 = tmp0[ 1 ] * tmp25;
    const auto tmp136 = tmp0[ 1 ] * tmp45;
    const auto tmp137 = tmp136 + tmp13;
    const auto tmp138 = tmp0[ 1 ] * tmp64;
    const auto tmp139 = tmp0[ 1 ] * tmp75;
    const auto tmp140 = tmp139 + tmp13;
    const auto tmp141 = tmp0[ 1 ] * tmp88;
    const auto tmp142 = tmp141 + tmp33;
    const auto tmp143 = tmp142 + tmp33;
    const auto tmp144 = tmp0[ 1 ] * tmp100;
    const auto tmp145 = tmp144 + tmp52;
    const auto tmp146 = tmp0[ 1 ] * tmp110;
    const auto tmp147 = tmp0[ 1 ] * tmp121;
    const auto tmp148 = tmp147 + tmp52;
    const auto tmp149 = tmp0[ 1 ] * tmp133;
    const auto tmp150 = tmp0[ 2 ] * tmp25;
    const auto tmp151 = tmp0[ 2 ] * tmp45;
    const auto tmp152 = tmp0[ 2 ] * tmp64;
    const auto tmp153 = tmp152 + tmp13;
    const auto tmp154 = tmp0[ 2 ] * tmp75;
    const auto tmp155 = tmp0[ 2 ] * tmp88;
    const auto tmp156 = tmp0[ 2 ] * tmp100;
    const auto tmp157 = tmp156 + tmp33;
    const auto tmp158 = tmp0[ 2 ] * tmp110;
    const auto tmp159 = tmp158 + tmp13;
    const auto tmp160 = tmp0[ 2 ] * tmp121;
    const auto tmp161 = tmp160 + tmp33;
    const auto tmp162 = tmp0[ 2 ] * tmp133;
    const auto tmp163 = tmp162 + tmp52;
    const auto tmp164 = tmp163 + tmp52;
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp28;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp47;
    ((result[ 0 ])[ 0 ])[ 2 ] = tmp66;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp77;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp89;
    ((result[ 0 ])[ 1 ])[ 2 ] = tmp101;
    ((result[ 0 ])[ 2 ])[ 0 ] = tmp112;
    ((result[ 0 ])[ 2 ])[ 1 ] = tmp122;
    ((result[ 0 ])[ 2 ])[ 2 ] = tmp134;
    ((result[ 1 ])[ 0 ])[ 0 ] = tmp135;
    ((result[ 1 ])[ 0 ])[ 1 ] = tmp137;
    ((result[ 1 ])[ 0 ])[ 2 ] = tmp138;
    ((result[ 1 ])[ 1 ])[ 0 ] = tmp140;
    ((result[ 1 ])[ 1 ])[ 1 ] = tmp143;
    ((result[ 1 ])[ 1 ])[ 2 ] = tmp145;
    ((result[ 1 ])[ 2 ])[ 0 ] = tmp146;
    ((result[ 1 ])[ 2 ])[ 1 ] = tmp148;
    ((result[ 1 ])[ 2 ])[ 2 ] = tmp149;
    ((result[ 2 ])[ 0 ])[ 0 ] = tmp150;
    ((result[ 2 ])[ 0 ])[ 1 ] = tmp151;
    ((result[ 2 ])[ 0 ])[ 2 ] = tmp153;
    ((result[ 2 ])[ 1 ])[ 0 ] = tmp154;
    ((result[ 2 ])[ 1 ])[ 1 ] = tmp155;
    ((result[ 2 ])[ 1 ])[ 2 ] = tmp157;
    ((result[ 2 ])[ 2 ])[ 0 ] = tmp159;
    ((result[ 2 ])[ 2 ])[ 1 ] = tmp161;
    ((result[ 2 ])[ 2 ])[ 2 ] = tmp164;
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
  ConstantTupleType constants_;
  std::tuple<  > coefficients_;
};

} // namespace UFLLocalFunctions_c359ec539158ceaebf6442fb226f656f

PYBIND11_MODULE( localfunction_c359ec539158ceaebf6442fb226f656f_d3813aca2f63526b9a4d82f7fa6d93d1, module )
{
  typedef UFLLocalFunctions_c359ec539158ceaebf6442fb226f656f::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_c359ec539158ceaebf6442fb226f656f::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_c359ec539158ceaebf6442fb226f656f_d3813aca2f63526b9a4d82f7fa6d93d1.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
