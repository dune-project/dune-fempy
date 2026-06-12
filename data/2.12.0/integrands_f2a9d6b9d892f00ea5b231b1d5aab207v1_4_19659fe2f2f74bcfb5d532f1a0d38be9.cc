#ifndef GuardIntegrands_f2a9d6b9d892f00ea5b231b1d5aab207v1_4
#define GuardIntegrands_f2a9d6b9d892f00ea5b231b1d5aab207v1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <cmath>
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

namespace Integrands_f2a9d6b9d892f00ea5b231b1d5aab207v1_4
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffu_h >
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
    typedef double Coneps;
    typedef std::tuple< std::shared_ptr< Coneps > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffu_h > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffu_h>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffu_h &coeffu_h, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< Coeffu_h >( coeffu_h ) )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Coneps >( (Coneps(0)) );
    }

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< Coneps >( (constant<0>()) );
    }

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      geometry_.emplace( this->entity().geometry() );
      std::get< 0 >( coefficients_ ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_ ).unbind();
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
      using std::pow;
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
        assert(tmp0 == tmp0);
      const auto tmp1 = 1 + tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::pow( tmp1, 2 );
        assert(tmp2 == tmp2);
      const auto tmp3 = std::get< 0 >( u );
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3[ 0 ] * tmp2;
        assert(tmp4 == tmp4);
      GlobalCoordinateType tmp5 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5[ 1 ] * tmp5[ 1 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp5[ 0 ] * tmp5[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 + tmp6;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 + tmp4;
        assert(tmp10 == tmp10);
      double tmp11 = constant< 0 >();
        assert(tmp11 == tmp11);
      const auto tmp12 = std::get< 1 >( u );
        assert(tmp12 == tmp12);
      const auto tmp13 = (tmp12[ 0 ])[ 0 ] * tmp11;
        assert(tmp13 == tmp13);
      const auto tmp14 = (tmp12[ 0 ])[ 1 ] * tmp11;
        assert(tmp14 == tmp14);
      return RangeValueType{ { tmp10 }, { { tmp13, tmp14 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
        assert(tmp0 == tmp0);
      const auto tmp1 = 1 + tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::pow( tmp1, 2 );
        assert(tmp2 == tmp2);
      double tmp3 = constant< 0 >();
        assert(tmp3 == tmp3);
      return [ tmp2, tmp3 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp2 * (std::get< 0 >( phi ))[ 0 ] }, { { tmp3 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], tmp3 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
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

    const Coneps &coneps () const
    {
      return *std::get< 0 >( constants_ );
    }

    Coneps &coneps ()
    {
      return *std::get< 0 >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_;
    }

    const Geometry &geometry () const
    {
      return *geometry_;
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
    std::optional< Geometry > geometry_;
    ConstantTupleType constants_;
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffu_h > > coefficients_;
  };

} // namespace Integrands_f2a9d6b9d892f00ea5b231b1d5aab207v1_4

PYBIND11_MODULE( integrands_f2a9d6b9d892f00ea5b231b1d5aab207v1_4_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_f2a9d6b9d892f00ea5b231b1d5aab207v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_f2a9d6b9d892f00ea5b231b1d5aab207v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_f2a9d6b9d892f00ea5b231b1d5aab207v1_4_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffu_h ) { return new Integrands( coeffu_h ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property( "eps", [] ( Integrands &self ) -> Integrands::Coneps { return self.coneps(); }, [] ( Integrands &self, const Integrands::Coneps &v ) { self.coneps() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_f2a9d6b9d892f00ea5b231b1d5aab207v1_4
