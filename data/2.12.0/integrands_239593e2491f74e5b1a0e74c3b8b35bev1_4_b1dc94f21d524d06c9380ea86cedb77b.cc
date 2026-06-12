#ifndef GuardIntegrands_239593e2491f74e5b1a0e74c3b8b35bev1_4
#define GuardIntegrands_239593e2491f74e5b1a0e74c3b8b35bev1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/alugrid/hierarchical.hh>
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

namespace Integrands_239593e2491f74e5b1a0e74c3b8b35bev1_4
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffu_h_n, class Coeffp_h_n >
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
    typedef double Contau;
    typedef std::tuple< std::shared_ptr< Contau > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffu_h_n, Coeffp_h_n > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffu_h_n>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffp_h_n>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffu_h_n &coeffu_h_n, const Coeffp_h_n &coeffp_h_n, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< Coeffu_h_n >( coeffu_h_n ), Dune::Fem::ConstLocalFunction< Coeffp_h_n >( coeffp_h_n ) )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Contau >( (Contau(0)) );
    }

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< Contau >( (constant<0>()) );
    }

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
    typedef std::tuple< Dune::FieldMatrix< double, 1, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = 1 / tmp0;
        assert(tmp1 == tmp1);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp2 = jacobianCoefficient< 0 >( x );
        assert(tmp2 == tmp2);
      const auto tmp3 = (tmp2[ 0 ])[ 0 ] + (tmp2[ 1 ])[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 * tmp1;
        assert(tmp4 == tmp4);
      typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp5 = jacobianCoefficient< 1 >( x );
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * (tmp5[ 0 ])[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = std::get< 0 >( u );
        assert(tmp7 == tmp7);
      const auto tmp8 = (tmp7[ 0 ])[ 0 ] + tmp6;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * (tmp5[ 0 ])[ 1 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = (tmp7[ 0 ])[ 1 ] + tmp9;
        assert(tmp10 == tmp10);
      return RangeValueType{ { tmp4 }, { { tmp8, tmp10 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { ((std::get< 0 >( phi ))[ 0 ])[ 0 ], ((std::get< 0 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }
    typedef Dune::FieldVector< double, 1 >  RRangeType;
    typedef Dune::FieldMatrix< double, 1, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,1> DirichletComponentType;

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
                assert(tmp0 == tmp0);
              const auto tmp1 = -1 * tmp0[ 0 ];
                assert(tmp1 == tmp1);
              const auto tmp2 = 2.2 + tmp1;
                assert(tmp2 == tmp2);
              domainId = tmp2 <= 1e-08;
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

    const Contau &contau () const
    {
      return *std::get< 0 >( constants_ );
    }

    Contau &contau ()
    {
      return *std::get< 0 >( constants_ );
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
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffu_h_n >, Dune::Fem::ConstLocalFunction< Coeffp_h_n > > coefficients_;
  };

} // namespace Integrands_239593e2491f74e5b1a0e74c3b8b35bev1_4

PYBIND11_MODULE( integrands_239593e2491f74e5b1a0e74c3b8b35bev1_4_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > GridPart;
  typedef Integrands_239593e2491f74e5b1a0e74c3b8b35bev1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_239593e2491f74e5b1a0e74c3b8b35bev1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_239593e2491f74e5b1a0e74c3b8b35bev1_4_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > > &coeffu_h_n, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffp_h_n ) { return new Integrands( coeffu_h_n, coeffp_h_n ); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property( "tau", [] ( Integrands &self ) -> Integrands::Contau { return self.contau(); }, [] ( Integrands &self, const Integrands::Contau &v ) { self.contau() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_239593e2491f74e5b1a0e74c3b8b35bev1_4
