
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCPluginFactory.hh"
#include "NCPhysicsModel.hh"

namespace NCPluginNamespace {

  class PluginScatter final : public NC::ProcImpl::ScatterIsotropicMat {
  public:

    // Adapt PhysicsModel to NCrystal's isotropic scatter interface.

    const char * name() const noexcept override
    {
      return NCPLUGIN_NAME_CSTR "Model";
    }

    PluginScatter( PhysicsModel && pm ) : m_pm(std::move(pm)) {}

    NC::CrossSect
    crossSectionIsotropic( NC::CachePtr& cache,
                           NC::NeutronEnergy ekin ) const override
    {
      return NC::CrossSect{ m_pm.calcCrossSection(cache, ekin.dbl()) };
    }

    NC::ScatterOutcomeIsotropic
    sampleScatterIsotropic( NC::CachePtr& cache,
                            NC::RNG& rng,
                            NC::NeutronEnergy ekin ) const override
    {
      auto outcome = m_pm.sampleScatteringEvent( cache, rng, ekin.dbl() );
      return { NC::NeutronEnergy{outcome.ekin_final},
               NC::CosineScatAngle{outcome.mu} };
    }

  private:
    PhysicsModel m_pm;
  };

}

const char * NCP::PluginFactory::name() const noexcept
{
  return NCPLUGIN_NAME_CSTR "Factory";
}

NC::Priority
NCP::PluginFactory::query( const NC::FactImpl::ScatterRequest& cfg ) const
{
  if ( cfg.get_inelas()=="none"
       || cfg.get_inelas()=="0"
       || cfg.get_inelas()=="false"
       || cfg.get_inelas()=="sterile" )
    return NC::Priority::Unable;

  if ( ! PhysicsModel::isApplicable( cfg.info() ) )
    return NC::Priority::Unable;

  // Override NCrystal's standard inelastic factory.
  return NC::Priority{999};
}

NC::ProcImpl::ProcPtr
NCP::PluginFactory::produce( const NC::FactImpl::ScatterRequest& cfg ) const
{
  auto physics_model = PhysicsModel::createFromInfo( cfg.info() );
  const bool keep_standard_elastic = physics_model.usesSolid();
  auto sc_ourmodel
    = NC::makeSO<PluginScatter>( std::move(physics_model) );

  // DirectGA supplies the inelastic process. Retain NCrystal's standard
  // elastic processes for solids, whose DirectGA kernel is inelastic only.
  auto sc_std = globalCreateScatter( cfg.modified(
      keep_standard_elastic
        ? "inelas=0"
        : "inelas=0;incoh_elas=0;coh_elas=0" ) );

  return combineProcs( sc_std, sc_ourmodel );
}
