#ifndef GuardIntegrands_d7dbcbcedb2d33d16780a3b9c57b029av1_3
#define GuardIntegrands_d7dbcbcedb2d33d16780a3b9c57b029av1_3
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
#include <dune/fempy/py/integrands.hh>

namespace Integrands_d7dbcbcedb2d33d16780a3b9c57b029av1_3
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
    typedef Dune::FieldVector< double, 2 > Conforcing;
    typedef std::tuple< std::shared_ptr< Conforcing > > ConstantTupleType;
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
      std::get< 0 >( constants_ ) = std::make_shared< Conforcing >( (Conforcing(0)) );
    }

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      return true;
    }

    void unbind ()
    {}

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      return (intersection.boundary() && init( intersection.inside() ));
    }
    typedef std::tuple< Dune::FieldMatrix< double, 2, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      Dune::FieldVector< double, 2 > tmp0 = constant< 0 >();
      const auto tmp1 = -1 * tmp0[ 0 ];
      const auto tmp2 = -1 * tmp0[ 1 ];
      const auto tmp3 = std::get< 0 >( u );
      const auto tmp4 = (tmp3[ 0 ])[ 0 ] + (tmp3[ 0 ])[ 0 ];
      const auto tmp5 = 0.5 * tmp4;
      const auto tmp6 = 2 * tmp5;
      const auto tmp7 = (tmp3[ 0 ])[ 0 ] + (tmp3[ 1 ])[ 1 ];
      const auto tmp8 = 0.1 * tmp7;
      const auto tmp9 = tmp8 + tmp6;
      const auto tmp10 = (tmp3[ 0 ])[ 1 ] + (tmp3[ 1 ])[ 0 ];
      const auto tmp11 = 0.5 * tmp10;
      const auto tmp12 = 2 * tmp11;
      const auto tmp13 = 0.5 * tmp12;
      const auto tmp14 = tmp13 + tmp13;
      const auto tmp15 = (tmp3[ 1 ])[ 1 ] + (tmp3[ 1 ])[ 1 ];
      const auto tmp16 = 0.5 * tmp15;
      const auto tmp17 = 2 * tmp16;
      const auto tmp18 = tmp8 + tmp17;
      return RangeValueType{ { tmp1, tmp2 }, { { tmp9, tmp14 }, { tmp14, tmp18 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0, 0 }, { { 2.1 * ((std::get< 0 >( phi ))[ 0 ])[ 0 ] + 0.1 * ((std::get< 0 >( phi ))[ 1 ])[ 1 ], 1.0 * ((std::get< 0 >( phi ))[ 0 ])[ 1 ] + 1.0 * ((std::get< 0 >( phi ))[ 1 ])[ 0 ] }, { 1.0 * ((std::get< 0 >( phi ))[ 0 ])[ 1 ] + 1.0 * ((std::get< 0 >( phi ))[ 1 ])[ 0 ], 0.1 * ((std::get< 0 >( phi ))[ 0 ])[ 0 ] + 2.1 * ((std::get< 0 >( phi ))[ 1 ])[ 1 ] } } };
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
          int domainId;
          {
              auto tmp0 = intersection.geometry().center( );
              domainId = tmp0[ 0 ] < 1e-10;
              if (domainId)
              {
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 2 );
                return true;
              }
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
      case 2:
        {
          result[ 0 ] = 0.0;
          result[ 1 ] = 0.0;
        }
        break;
      default:
        {
          result = RRangeType( 0 );
        }
      }
    }

    template< class Point >
    void dDirichlet ( int bndId, const Point &x, RJacobianRangeType &result ) const
    {
      switch( bndId )
      {
      case 2:
        {
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = 0.0;
          (result[ 1 ])[ 0 ] = 0.0;
          (result[ 1 ])[ 1 ] = 0.0;
        }
        break;
      default:
        {
          result = RJacobianRangeType( 0 );
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

    const Conforcing &conforcing () const
    {
      return *std::get< 0 >( constants_ );
    }

    Conforcing &conforcing ()
    {
      return *std::get< 0 >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_;
    }
    EntityType entity_;
    IntersectionType intersection_;
    ConstantTupleType constants_;
    std::tuple<  > coefficients_;
  };

} // namespace Integrands_d7dbcbcedb2d33d16780a3b9c57b029av1_3

PYBIND11_MODULE( integrands_d7dbcbcedb2d33d16780a3b9c57b029av1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_d7dbcbcedb2d33d16780a3b9c57b029av1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_d7dbcbcedb2d33d16780a3b9c57b029av1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_d7dbcbcedb2d33d16780a3b9c57b029av1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property( "forcing", [] ( Integrands &self ) -> Integrands::Conforcing { return self.conforcing(); }, [] ( Integrands &self, const Integrands::Conforcing &v ) { self.conforcing() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_d7dbcbcedb2d33d16780a3b9c57b029av1_3
