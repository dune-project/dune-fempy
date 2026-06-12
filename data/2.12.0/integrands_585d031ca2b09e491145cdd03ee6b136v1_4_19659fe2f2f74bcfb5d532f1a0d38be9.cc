#ifndef GuardIntegrands_585d031ca2b09e491145cdd03ee6b136v1_4
#define GuardIntegrands_585d031ca2b09e491145cdd03ee6b136v1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_585d031ca2b09e491145cdd03ee6b136v1_4
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffcopy_of_u, class Coeffcopy_of_v >
  struct Integrands
  {
    typedef GridPart GridPartType;
    typedef typename GridPartType::GridViewType GridView;
    typedef typename GridView::ctype ctype;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename GridPartType::IntersectionType IntersectionType;
    typedef typename EntityType::Geometry Geometry;
    typedef typename Geometry::GlobalCoordinate GlobalCoordinateType;
    typedef Dune::Fem::IntersectionSide Side;
    typedef std::tuple<> ConstantTupleType;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffcopy_of_u, Coeffcopy_of_v > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffcopy_of_u>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffcopy_of_v>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffcopy_of_u &coeffcopy_of_u, const Coeffcopy_of_v &coeffcopy_of_v, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< Coeffcopy_of_u >( coeffcopy_of_u ), Dune::Fem::ConstLocalFunction< Coeffcopy_of_v >( coeffcopy_of_v ) )
    {}

    void decouple ()
    {}

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      std::get< 0 >( coefficients_ ).bind( this->entity() );
      std::get< 1 >( coefficients_ ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_ ).unbind();
      std::get< 1 >( coefficients_ ).unbind();
    }

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      return (intersection.boundary() && init( intersection.inside() ));
    }
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
        assert(tmp0 == tmp0);
      const auto tmp1 = 12.5 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::get< 0 >( u );
        assert(tmp2 == tmp2);
      const auto tmp3 = -1 * tmp2[ 0 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = 1 + tmp3;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 * tmp1;
        assert(tmp5 == tmp5);
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp6 = evaluateCoefficient< 1 >( x );
        assert(tmp6 == tmp6);
      const auto tmp7 = 0.02 + tmp6[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 / 0.75;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp0[ 0 ] + tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * tmp5;
        assert(tmp11 == tmp11);
      const auto tmp12 = 12.5 * tmp2[ 0 ];
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * tmp0[ 0 ];
        assert(tmp13 == tmp13);
      const auto tmp14 = 1 + tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp12;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp10 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = -1 * (tmp0[ 0 ] < tmp8 ? tmp16 : tmp11);
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp2[ 0 ] + tmp13;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp18 + tmp17;
        assert(tmp19 == tmp19);
      const auto tmp20 = std::get< 1 >( u );
        assert(tmp20 == tmp20);
      const auto tmp21 = 0.0025 * (tmp20[ 0 ])[ 0 ];
        assert(tmp21 == tmp21);
      const auto tmp22 = 0.0025 * (tmp20[ 0 ])[ 1 ];
        assert(tmp22 == tmp22);
      return RangeValueType{ { tmp19 }, { { tmp21, tmp22 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
        assert(tmp0 == tmp0);
      const auto tmp1 = 12.5 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = -1 * tmp1;
        assert(tmp2 == tmp2);
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp3 = evaluateCoefficient< 1 >( x );
        assert(tmp3 == tmp3);
      const auto tmp4 = 0.02 + tmp3[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 / 0.75;
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp0[ 0 ] + tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 * tmp2;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp0[ 0 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = 1 + tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = 12.5 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp7 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * (tmp0[ 0 ] < tmp5 ? tmp12 : tmp8);
        assert(tmp13 == tmp13);
      const auto tmp14 = 1 + tmp13;
        assert(tmp14 == tmp14);
      return [ tmp14 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp14 * (std::get< 0 >( phi ))[ 0 ] }, { { 0.0025 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], 0.0025 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }
    typedef Dune::FieldVector< double, 1 >  RRangeType;
    typedef Dune::FieldMatrix< double, 1, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,1> DirichletComponentType;

    bool hasDirichletBoundary () const
    {
      return false;
    }

    bool isDirichletIntersection ( const IntersectionType &intersection, DirichletComponentType &dirichletComponent ) const
    {
      return false;
    }

    template< class Point >
    void dirichlet ( int bndId, const Point &x, RRangeType &result ) const
    {}

    template< class Point >
    void dDirichlet ( int bndId, const Point &x, RJacobianRangeType &result ) const
    {}

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

    const EntityType &entity () const
    {
      return entity_;
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
    EntityType entity_;
    IntersectionType intersection_;
    ConstantTupleType constants_;
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffcopy_of_u >, Dune::Fem::ConstLocalFunction< Coeffcopy_of_v > > coefficients_;
  };

} // namespace Integrands_585d031ca2b09e491145cdd03ee6b136v1_4

PYBIND11_MODULE( integrands_585d031ca2b09e491145cdd03ee6b136v1_4_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_585d031ca2b09e491145cdd03ee6b136v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_585d031ca2b09e491145cdd03ee6b136v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_585d031ca2b09e491145cdd03ee6b136v1_4_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffcopy_of_u, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffcopy_of_v ) { return new Integrands( coeffcopy_of_u, coeffcopy_of_v ); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_585d031ca2b09e491145cdd03ee6b136v1_4
