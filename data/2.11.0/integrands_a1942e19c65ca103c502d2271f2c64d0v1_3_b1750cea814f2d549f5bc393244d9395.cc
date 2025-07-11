#ifndef GuardIntegrands_a1942e19c65ca103c502d2271f2c64d0v1_3
#define GuardIntegrands_a1942e19c65ca103c502d2271f2c64d0v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_a1942e19c65ca103c502d2271f2c64d0v1_3
{




  // Integrands
  // ----------

  template< class GridPart >
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
    typedef double Conmu;
    typedef double Connu;
    typedef std::tuple< std::shared_ptr< Conmu >, std::shared_ptr< Connu > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple<> CoefficientTupleType;
    static constexpr bool gridPartValid = true;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Conmu >( (Conmu(0)) );
      std::get< 1 >( constants_ ) = std::make_shared< Connu >( (Connu(0)) );
    }

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      geometry_.emplace( this->entity().geometry() );
      return true;
    }

    void unbind ()
    {}

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      return (intersection.boundary() && init( intersection.inside() ));
    }
    typedef std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = -1 * tmp0[ 1 ];
      const auto tmp2 = 1.0 + tmp1;
      const auto tmp3 = tmp0[ 1 ] * tmp2;
      double tmp4 = constant< 1 >();
      const auto tmp5 = tmp4 * tmp3;
      const auto tmp6 = -1 * tmp5;
      const auto tmp7 = std::get< 0 >( u );
      const auto tmp8 = tmp7[ 0 ] * tmp4;
      const auto tmp9 = tmp8 + tmp6;
      const auto tmp10 = tmp7[ 1 ] * tmp4;
      const auto tmp11 = std::get< 1 >( u );
      const auto tmp12 = (tmp11[ 0 ])[ 0 ] + (tmp11[ 0 ])[ 0 ];
      double tmp13 = constant< 0 >();
      const auto tmp14 = tmp13 * tmp12;
      const auto tmp15 = (tmp11[ 0 ])[ 1 ] + (tmp11[ 1 ])[ 0 ];
      const auto tmp16 = tmp13 * tmp15;
      const auto tmp17 = (tmp11[ 1 ])[ 1 ] + (tmp11[ 1 ])[ 1 ];
      const auto tmp18 = tmp13 * tmp17;
      return RangeValueType{ { tmp9, tmp10 }, { { tmp14, tmp16 }, { tmp16, tmp18 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 1 >();
      double tmp1 = constant< 0 >();
      const auto tmp2 = 2 * tmp1;
      return [ tmp0, tmp1, tmp2 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp0 * (std::get< 0 >( phi ))[ 0 ], tmp0 * (std::get< 0 >( phi ))[ 1 ] }, { { tmp2 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], tmp1 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] + tmp1 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ] }, { tmp1 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] + tmp1 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ], tmp2 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ] } } };
        };
    }
    typedef Dune::FieldVector< double, 2 >  RRangeType;
    typedef Dune::FieldMatrix< double, 2, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,2> DirichletComponentType;

    bool hasDirichletBoundary () const
    {
      return true;
    }

    bool isDirichletIntersection ( const IntersectionType &intersection, DirichletComponentType &dirichletComponent ) const
    {
      const int bndId = BoundaryIdProviderType::boundaryId( intersection );
      {
        std::fill( dirichletComponent.begin(), dirichletComponent.end(), bndId );
      }
      switch( bndId )
      {
      default:
        {
          {
            std::fill( dirichletComponent.begin(), dirichletComponent.end(), 1 );
            return true;
          }
          return false;
        }
      }
    }

    template< class Point >
    void dirichlet ( int bndId, const Point &x, RRangeType &result ) const
    {
      switch( bndId )
      {
      default:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = -1 * tmp0[ 1 ];
          const auto tmp2 = 1.0 + tmp1;
          const auto tmp3 = tmp0[ 1 ] * tmp2;
          result[ 0 ] = tmp3;
          result[ 1 ] = 0.0;
        }
      }
    }

    template< class Point >
    void dDirichlet ( int bndId, const Point &x, RJacobianRangeType &result ) const
    {
      switch( bndId )
      {
      default:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = -1 * tmp0[ 1 ];
          const auto tmp2 = 1.0 + tmp1;
          const auto tmp3 = tmp2 + tmp1;
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = tmp3;
          (result[ 1 ])[ 0 ] = 0.0;
          (result[ 1 ])[ 1 ] = 0.0;
        }
      }
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

    const Conmu &conmu () const
    {
      return *std::get< 0 >( constants_ );
    }

    Conmu &conmu ()
    {
      return *std::get< 0 >( constants_ );
    }

    const Connu &connu () const
    {
      return *std::get< 1 >( constants_ );
    }

    Connu &connu ()
    {
      return *std::get< 1 >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_;
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }
    EntityType entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    ConstantTupleType constants_;
    std::tuple<  > coefficients_;
  };

} // namespace Integrands_a1942e19c65ca103c502d2271f2c64d0v1_3

PYBIND11_MODULE( integrands_a1942e19c65ca103c502d2271f2c64d0v1_3_b1750cea814f2d549f5bc393244d9395, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::cube >::LeafGridView > GridPart;
  typedef Integrands_a1942e19c65ca103c502d2271f2c64d0v1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_a1942e19c65ca103c502d2271f2c64d0v1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_a1942e19c65ca103c502d2271f2c64d0v1_3_b1750cea814f2d549f5bc393244d9395.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property( "mu", [] ( Integrands &self ) -> Integrands::Conmu { return self.conmu(); }, [] ( Integrands &self, const Integrands::Conmu &v ) { self.conmu() = v; } );
      cls.def_property( "nu", [] ( Integrands &self ) -> Integrands::Connu { return self.connu(); }, [] ( Integrands &self, const Integrands::Connu &v ) { self.connu() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_a1942e19c65ca103c502d2271f2c64d0v1_3
