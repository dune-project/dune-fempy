#ifndef GuardIntegrands_04f4f97faccd7a0aad67c4dc02cee589v1_3
#define GuardIntegrands_04f4f97faccd7a0aad67c4dc02cee589v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/geometrygridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_04f4f97faccd7a0aad67c4dc02cee589v1_3
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffcopy_of_uh >
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
    typedef double ContimeStep;
    typedef std::tuple< std::shared_ptr< ContimeStep > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 3 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffcopy_of_uh > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffcopy_of_uh>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffcopy_of_uh &coeffcopy_of_uh, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< Coeffcopy_of_uh >( coeffcopy_of_uh ) )
    {
      std::get< 0 >( constants_ ) = std::make_shared< ContimeStep >( (ContimeStep(0)) );
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
    typedef std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 3 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 3 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
      const auto tmp1 = -1 * tmp0[ 0 ];
      const auto tmp2 = std::get< 0 >( u );
      const auto tmp3 = tmp2[ 0 ] + tmp1;
      const auto tmp4 = -1 * tmp0[ 1 ];
      const auto tmp5 = tmp2[ 1 ] + tmp4;
      const auto tmp6 = -1 * tmp0[ 2 ];
      const auto tmp7 = tmp2[ 2 ] + tmp6;
      double tmp8 = constant< 0 >();
      const auto tmp9 = std::get< 1 >( u );
      const auto tmp10 = (tmp9[ 0 ])[ 0 ] * tmp8;
      const auto tmp11 = (tmp9[ 0 ])[ 1 ] * tmp8;
      const auto tmp12 = (tmp9[ 0 ])[ 2 ] * tmp8;
      const auto tmp13 = (tmp9[ 1 ])[ 0 ] * tmp8;
      const auto tmp14 = (tmp9[ 1 ])[ 1 ] * tmp8;
      const auto tmp15 = (tmp9[ 1 ])[ 2 ] * tmp8;
      const auto tmp16 = (tmp9[ 2 ])[ 0 ] * tmp8;
      const auto tmp17 = (tmp9[ 2 ])[ 1 ] * tmp8;
      const auto tmp18 = (tmp9[ 2 ])[ 2 ] * tmp8;
      return RangeValueType{ { tmp3, tmp5, tmp7 }, { { tmp10, tmp11, tmp12 }, { tmp13, tmp14, tmp15 }, { tmp16, tmp17, tmp18 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 0 >();
      return [ tmp0 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (std::get< 0 >( phi ))[ 0 ], (std::get< 0 >( phi ))[ 1 ], (std::get< 0 >( phi ))[ 2 ] }, { { tmp0 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], tmp0 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], tmp0 * ((std::get< 1 >( phi ))[ 0 ])[ 2 ] }, { tmp0 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ], tmp0 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ], tmp0 * ((std::get< 1 >( phi ))[ 1 ])[ 2 ] }, { tmp0 * ((std::get< 1 >( phi ))[ 2 ])[ 0 ], tmp0 * ((std::get< 1 >( phi ))[ 2 ])[ 1 ], tmp0 * ((std::get< 1 >( phi ))[ 2 ])[ 2 ] } } };
        };
    }
    typedef Dune::FieldVector< double, 3 >  RRangeType;
    typedef Dune::FieldMatrix< double, 3, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,3> DirichletComponentType;

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
          result[ 0 ] = tmp0[ 0 ];
          result[ 1 ] = tmp0[ 1 ];
          result[ 2 ] = tmp0[ 2 ];
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
          (result[ 0 ])[ 0 ] = 1;
          (result[ 0 ])[ 1 ] = 0.0;
          (result[ 0 ])[ 2 ] = 0.0;
          (result[ 1 ])[ 0 ] = 0.0;
          (result[ 1 ])[ 1 ] = 1;
          (result[ 1 ])[ 2 ] = 0.0;
          (result[ 2 ])[ 0 ] = 0.0;
          (result[ 2 ])[ 1 ] = 0.0;
          (result[ 2 ])[ 2 ] = 1;
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

    const ContimeStep &contimeStep () const
    {
      return *std::get< 0 >( constants_ );
    }

    ContimeStep &contimeStep ()
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
    EntityType entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    ConstantTupleType constants_;
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffcopy_of_uh > > coefficients_;
  };

} // namespace Integrands_04f4f97faccd7a0aad67c4dc02cee589v1_3

PYBIND11_MODULE( integrands_04f4f97faccd7a0aad67c4dc02cee589v1_3_981d293fc017fb4ad04dc4ea0cca632b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::GeometryGridPart< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > > > GridPart;
  typedef Integrands_04f4f97faccd7a0aad67c4dc02cee589v1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::GeometryGridPart< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > >>, Dune::FieldVector< double, 3 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_04f4f97faccd7a0aad67c4dc02cee589v1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::GeometryGridPart< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > >>, Dune::FieldVector< double, 3 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_04f4f97faccd7a0aad67c4dc02cee589v1_3_981d293fc017fb4ad04dc4ea0cca632b.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::GeometryGridPart< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > >>, Dune::FieldVector< double, 3 > > &coeffcopy_of_uh ) { return new Integrands( coeffcopy_of_uh ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property( "timeStep", [] ( Integrands &self ) -> Integrands::ContimeStep { return self.contimeStep(); }, [] ( Integrands &self, const Integrands::ContimeStep &v ) { self.contimeStep() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_04f4f97faccd7a0aad67c4dc02cee589v1_3
