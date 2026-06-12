#ifndef GuardIntegrands_1f73213e700a160a8b1134dbaf710657v1_4
#define GuardIntegrands_1f73213e700a160a8b1134dbaf710657v1_4
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

namespace Integrands_1f73213e700a160a8b1134dbaf710657v1_4
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
    typedef double Conmu;
    typedef double Conrho;
    typedef Dune::FieldVector< double, 2 > Conf;
    typedef std::tuple< std::shared_ptr< Contau >, std::shared_ptr< Conmu >, std::shared_ptr< Conrho >, std::shared_ptr< Conf > > ConstantTupleType;
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
      std::get< 1 >( constants_ ) = std::make_shared< Conmu >( (Conmu(0)) );
      std::get< 2 >( constants_ ) = std::make_shared< Conrho >( (Conrho(0)) );
      std::get< 3 >( constants_ ) = std::make_shared< Conf >( (Conf(0)) );
    }

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< Contau >( (constant<0>()) );
      std::get< 1 >( constants_ ) = std::make_shared< Conmu >( (constant<1>()) );
      std::get< 2 >( constants_ ) = std::make_shared< Conrho >( (constant<2>()) );
      std::get< 3 >( constants_ ) = std::make_shared< Conf >( (constant<3>()) );
    }

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      geometry_.emplace( this->entity().geometry() );
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
      double tmp0 = constant< 2 >();
        assert(tmp0 == tmp0);
      Dune::FieldVector< double, 2 > tmp1 = constant< 3 >();
        assert(tmp1 == tmp1);
      const auto tmp2 = tmp1[ 0 ] * tmp0;
        assert(tmp2 == tmp2);
      const auto tmp3 = -1 * tmp2;
        assert(tmp3 == tmp3);
      double tmp4 = constant< 0 >();
        assert(tmp4 == tmp4);
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp5 = evaluateCoefficient< 0 >( x );
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = std::get< 0 >( u );
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7[ 0 ] + tmp6;
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp8 / tmp4;
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 * tmp0;
        assert(tmp10 == tmp10);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp11 = jacobianCoefficient< 0 >( x );
        assert(tmp11 == tmp11);
      const auto tmp12 = (tmp11[ 0 ])[ 1 ] * tmp5[ 1 ];
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp5[ 0 ] * (tmp11[ 0 ])[ 0 ];
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp13 + tmp12;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp0;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 + tmp10;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 + tmp3;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp1[ 1 ] * tmp0;
        assert(tmp18 == tmp18);
      const auto tmp19 = -1 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = -1 * tmp5[ 1 ];
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp7[ 1 ] + tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 / tmp4;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 * tmp0;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp5[ 1 ] * (tmp11[ 1 ])[ 1 ];
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp5[ 0 ] * (tmp11[ 1 ])[ 0 ];
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp25 + tmp24;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 * tmp0;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 + tmp23;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 + tmp19;
        assert(tmp29 == tmp29);
      const auto tmp30 = std::get< 1 >( u );
        assert(tmp30 == tmp30);
      const auto tmp31 = (tmp30[ 0 ])[ 0 ] + (tmp11[ 0 ])[ 0 ];
        assert(tmp31 == tmp31);
      const auto tmp32 = 0.5 * tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 + tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp33 / 2;
        assert(tmp34 == tmp34);
      double tmp35 = constant< 1 >();
        assert(tmp35 == tmp35);
      const auto tmp36 = 2 * tmp35;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 * tmp34;
        assert(tmp37 == tmp37);
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp38 = evaluateCoefficient< 1 >( x );
        assert(tmp38 == tmp38);
      const auto tmp39 = -1 * tmp38[ 0 ];
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 + tmp37;
        assert(tmp40 == tmp40);
      const auto tmp41 = (tmp30[ 1 ])[ 0 ] + (tmp11[ 1 ])[ 0 ];
        assert(tmp41 == tmp41);
      const auto tmp42 = 0.5 * tmp41;
        assert(tmp42 == tmp42);
      const auto tmp43 = (tmp30[ 0 ])[ 1 ] + (tmp11[ 0 ])[ 1 ];
        assert(tmp43 == tmp43);
      const auto tmp44 = 0.5 * tmp43;
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp44 + tmp42;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 / 2;
        assert(tmp46 == tmp46);
      const auto tmp47 = tmp36 * tmp46;
        assert(tmp47 == tmp47);
      const auto tmp48 = 0.5 * tmp47;
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp48 + tmp48;
        assert(tmp49 == tmp49);
      const auto tmp50 = (tmp30[ 1 ])[ 1 ] + (tmp11[ 1 ])[ 1 ];
        assert(tmp50 == tmp50);
      const auto tmp51 = 0.5 * tmp50;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 + tmp51;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp52 / 2;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp36 * tmp53;
        assert(tmp54 == tmp54);
      const auto tmp55 = tmp39 + tmp54;
        assert(tmp55 == tmp55);
      return RangeValueType{ { tmp17, tmp29 }, { { tmp40, tmp49 }, { tmp49, tmp55 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 2 >();
        assert(tmp0 == tmp0);
      double tmp1 = constant< 0 >();
        assert(tmp1 == tmp1);
      const auto tmp2 = 1 / tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 * tmp0;
        assert(tmp3 == tmp3);
      double tmp4 = constant< 1 >();
        assert(tmp4 == tmp4);
      const auto tmp5 = 2 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = 0.5 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 0.25 * tmp5;
        assert(tmp7 == tmp7);
      const auto tmp8 = 0.5 * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp8 + tmp8;
        assert(tmp9 == tmp9);
      return [ tmp3, tmp6, tmp9 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp3 * (std::get< 0 >( phi ))[ 0 ], tmp3 * (std::get< 0 >( phi ))[ 1 ] }, { { tmp6 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], tmp9 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] + tmp9 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ] }, { tmp9 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] + tmp9 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ], tmp6 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ] } } };
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
                assert(tmp0 == tmp0);
              domainId = tmp0[ 0 ] <= 1e-08;
              if (domainId)
              {
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 2 );
                return true;
              }
          }
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
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 3 );
                dirichletComponent[0] = 0;
                dirichletComponent[1] = 0;
                return true;
              }
          }
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
      case 2:
        {
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = 6 * tmp0[ 1 ];
            assert(tmp1 == tmp1);
          const auto tmp2 = -1 * tmp0[ 1 ];
            assert(tmp2 == tmp2);
          const auto tmp3 = 0.41 + tmp2;
            assert(tmp3 == tmp3);
          const auto tmp4 = tmp3 * tmp1;
            assert(tmp4 == tmp4);
          const auto tmp5 = tmp4 / 0.16809999999999997;
            assert(tmp5 == tmp5);
          result[ 0 ] = tmp5;
          result[ 1 ] = 0.0;
        }
        break;
      case 3:
        {
          result[ 0 ] = 0.0;
          result[ 1 ] = 0.0;
        }
        break;
      default:
        {
          result[ 0 ] = 0.0;
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
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = 6 * tmp0[ 1 ];
            assert(tmp1 == tmp1);
          const auto tmp2 = -1 * tmp1;
            assert(tmp2 == tmp2);
          const auto tmp3 = -1 * tmp0[ 1 ];
            assert(tmp3 == tmp3);
          const auto tmp4 = 0.41 + tmp3;
            assert(tmp4 == tmp4);
          const auto tmp5 = 6 * tmp4;
            assert(tmp5 == tmp5);
          const auto tmp6 = tmp5 + tmp2;
            assert(tmp6 == tmp6);
          const auto tmp7 = tmp6 / 0.16809999999999997;
            assert(tmp7 == tmp7);
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = tmp7;
          (result[ 1 ])[ 0 ] = 0.0;
          (result[ 1 ])[ 1 ] = 0.0;
        }
        break;
      case 3:
        {
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = 0.0;
          (result[ 1 ])[ 0 ] = 0.0;
          (result[ 1 ])[ 1 ] = 0.0;
        }
        break;
      default:
        {
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = 0.0;
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

    const Contau &contau () const
    {
      return *std::get< 0 >( constants_ );
    }

    Contau &contau ()
    {
      return *std::get< 0 >( constants_ );
    }

    const Conmu &conmu () const
    {
      return *std::get< 1 >( constants_ );
    }

    Conmu &conmu ()
    {
      return *std::get< 1 >( constants_ );
    }

    const Conrho &conrho () const
    {
      return *std::get< 2 >( constants_ );
    }

    Conrho &conrho ()
    {
      return *std::get< 2 >( constants_ );
    }

    const Conf &conf () const
    {
      return *std::get< 3 >( constants_ );
    }

    Conf &conf ()
    {
      return *std::get< 3 >( constants_ );
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
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffu_h_n >, Dune::Fem::ConstLocalFunction< Coeffp_h_n > > coefficients_;
  };

} // namespace Integrands_1f73213e700a160a8b1134dbaf710657v1_4

PYBIND11_MODULE( integrands_1f73213e700a160a8b1134dbaf710657v1_4_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > GridPart;
  typedef Integrands_1f73213e700a160a8b1134dbaf710657v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_1f73213e700a160a8b1134dbaf710657v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_1f73213e700a160a8b1134dbaf710657v1_4_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 2 > > &coeffu_h_n, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffp_h_n ) { return new Integrands( coeffu_h_n, coeffp_h_n ); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property( "tau", [] ( Integrands &self ) -> Integrands::Contau { return self.contau(); }, [] ( Integrands &self, const Integrands::Contau &v ) { self.contau() = v; } );
      cls.def_property( "mu", [] ( Integrands &self ) -> Integrands::Conmu { return self.conmu(); }, [] ( Integrands &self, const Integrands::Conmu &v ) { self.conmu() = v; } );
      cls.def_property( "rho", [] ( Integrands &self ) -> Integrands::Conrho { return self.conrho(); }, [] ( Integrands &self, const Integrands::Conrho &v ) { self.conrho() = v; } );
      cls.def_property( "f", [] ( Integrands &self ) -> Integrands::Conf { return self.conf(); }, [] ( Integrands &self, const Integrands::Conf &v ) { self.conf() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_1f73213e700a160a8b1134dbaf710657v1_4
