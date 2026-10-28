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

#ifndef NCPlugin_PhysicsModel_hh
#define NCPlugin_PhysicsModel_hh

#include "NCrystal/NCPluginBoilerplate.hh"
#include "NCYKSkold.hh"
#include <unordered_map>

namespace NCPluginNamespace {

  // Parsed settings for one @CUSTOM_DIRECTGA element block.
  struct DirectGAInfo {
    bool solid = false;
    double max_energy_transfer_meV = 2000.0;
    double numerical_time_end_ps = 2.0; // Liquid default; solids use 5 ps.
    double time_step_ps = 0.0;
    bool time_step_auto = true;
    NC::VectD energy_meV;
    NC::VectD density_meV_inv; // f(E) in meV^-1 per atom
    YKSkoldInfo corrections;
  };

  class PhysicsModel final : public NC::MoveOnly {
  public:
    static bool isApplicable( const NC::Info& );
    static PhysicsModel createFromInfo( const NC::Info& );

    PhysicsModel( std::unordered_map<std::string, DirectGAInfo>,
                  const NC::Info& );

    bool usesSolid() const noexcept;
    double calcCrossSection( NC::CachePtr&, double neutron_ekin ) const;

    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::CachePtr&, NC::RNG& rng,
                                     double neutron_ekin ) const;

  private:
    std::unordered_map<std::string, DirectGAInfo> m_direct_ga_data;
    NC::ProcImpl::OptionalProcPtr m_proc;
  };

}
#endif
