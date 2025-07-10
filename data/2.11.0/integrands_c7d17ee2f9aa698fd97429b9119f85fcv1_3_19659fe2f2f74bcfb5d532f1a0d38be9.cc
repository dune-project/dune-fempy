#ifndef GuardIntegrands_c7d17ee2f9aa698fd97429b9119f85fcv1_3
#define GuardIntegrands_c7d17ee2f9aa698fd97429b9119f85fcv1_3
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
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_c7d17ee2f9aa698fd97429b9119f85fcv1_3
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffbnd >
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
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffbnd > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffbnd>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffbnd &coeffbnd, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< Coeffbnd >( coeffbnd ) )
    {}

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
    typedef std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 > > RangeValueType;
    static constexpr bool _nonlinear = true;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      using std::sin;
      const auto tmp0 = std::get< 0 >( u );
      const auto tmp1 = 2 * tmp0[ 1 ];
      const auto tmp2 = -1 * tmp1;
      const auto tmp3 = 1 + tmp2;
      const auto tmp4 = 10 * tmp3;
      const auto tmp5 = tmp0[ 0 ] + tmp4;
      GlobalCoordinateType tmp6 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp7 = 6.283185307179586 * tmp6[ 1 ];
      const auto tmp8 = std::sin( tmp7 );
      const auto tmp9 = 10 * tmp8;
      const auto tmp10 = 2 * tmp0[ 0 ];
      const auto tmp11 = -1 * tmp10;
      const auto tmp12 = 1 + tmp11;
      const auto tmp13 = std::pow( tmp12, 2 );
      const auto tmp14 = tmp0[ 1 ] * tmp13;
      const auto tmp15 = tmp14 + tmp9;
      const auto tmp16 = std::get< 1 >( u );
      return RangeValueType{ { tmp5, tmp15 }, { { (tmp16[ 0 ])[ 0 ], (tmp16[ 0 ])[ 1 ] }, { (tmp16[ 1 ])[ 0 ], (tmp16[ 1 ])[ 1 ] } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      const auto tmp0 = std::get< 0 >( u );
      const auto tmp1 = 2 * tmp0[ 0 ];
      const auto tmp2 = -1 * tmp1;
      const auto tmp3 = 1 + tmp2;
      const auto tmp4 = -4 * tmp3;
      const auto tmp5 = tmp0[ 1 ] * tmp4;
      const auto tmp6 = std::pow( tmp3, 2 );
      return [ tmp5, tmp6 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (std::get< 0 >( phi ))[ 0 ] + -20 * (std::get< 0 >( phi ))[ 1 ], tmp5 * (std::get< 0 >( phi ))[ 0 ] + tmp6 * (std::get< 0 >( phi ))[ 1 ] }, { { ((std::get< 1 >( phi ))[ 0 ])[ 0 ], ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { ((std::get< 1 >( phi ))[ 1 ])[ 0 ], ((std::get< 1 >( phi ))[ 1 ])[ 1 ] } } };
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
              domainId = tmp0[ 1 ] < 1e-10;
              if (domainId)
              {
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 2 );
                return true;
              }
          }
          {
            std::fill( dirichletComponent.begin(), dirichletComponent.end(), 1 );
            dirichletComponent[1] = 0;
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
      case 2:
        {
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 0 ] + tmp0[ 1 ];
          const auto tmp2 = 6.283185307179586 * tmp1;
          const auto tmp3 = std::sin( tmp2 );
          result[ 0 ] = tmp3;
          result[ 1 ] = 1;
        }
        break;
      default:
        {
          typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
          result[ 0 ] = tmp0[ 0 ];
          result[ 1 ] = 0.0;
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
          using std::cos;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = tmp0[ 0 ] + tmp0[ 1 ];
          const auto tmp2 = 6.283185307179586 * tmp1;
          const auto tmp3 = std::cos( tmp2 );
          const auto tmp4 = 6.283185307179586 * tmp3;
          (result[ 0 ])[ 0 ] = tmp4;
          (result[ 0 ])[ 1 ] = tmp4;
          (result[ 1 ])[ 0 ] = 0.0;
          (result[ 1 ])[ 1 ] = 0.0;
        }
        break;
      default:
        {
          typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp0 = jacobianCoefficient< 0 >( x );
          (result[ 0 ])[ 0 ] = (tmp0[ 0 ])[ 0 ];
          (result[ 0 ])[ 1 ] = (tmp0[ 0 ])[ 1 ];
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
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffbnd > > coefficients_;
  };

} // namespace Integrands_c7d17ee2f9aa698fd97429b9119f85fcv1_3

PYBIND11_MODULE( integrands_c7d17ee2f9aa698fd97429b9119f85fcv1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_c7d17ee2f9aa698fd97429b9119f85fcv1_3::Integrands< GridPart, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_c7d17ee2f9aa698fd97429b9119f85fcv1_3::Integrands< GridPart, Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_c7d17ee2f9aa698fd97429b9119f85fcv1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::Python::SimpleGridFunction< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, Dune::Python::detail::PyGridFunctionEvaluator< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView, 0, std::function<double(const Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView::template Codim<0>::Entity&,const Dune::FieldVector< double, 2 >&)> > > &coeffbnd ) { return new Integrands( coeffbnd ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_c7d17ee2f9aa698fd97429b9119f85fcv1_3
