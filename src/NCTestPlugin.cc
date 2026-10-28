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

#include "NCTestPlugin.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

#include <cmath>
#include <string>

// Check cross sections and sampled events for one liquid and one solid.
void NCP::customPluginTest()
{
  const double energies[] = {
    1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 0.5
  };

  struct TestCase {
    const char * material;
    const char * cfg_suffix;
  };

  const TestCase test_cases[] = {
    { "tsl-para-H_20K_directga.ncmat", ";temp=20K" },
    { "Al_sg225_directga.ncmat", "" }
  };

  for (const auto& test_case : test_cases) {
    auto scatter = NC::createScatter(
        std::string("plugins::DirectGA/") + test_case.material
        + test_case.cfg_suffix);
    for (double energy : energies) {
      const double xs = scatter.crossSectionIsotropic(
          NC::NeutronEnergy{energy}).dbl();
      if (!std::isfinite(xs) || !(xs > 0.0)) {
        NCPLUGIN_MSG("Cross-section check failed for "
                     << test_case.material << " at E=" << energy
                     << " eV: xs=" << xs);
        nc_assert_always(false);
      }

      for (unsigned sample = 0; sample < 10; ++sample) {
        const auto event = scatter.sampleScatterIsotropic(
            NC::NeutronEnergy{energy});
        const double output_energy = event.ekin.dbl();
        const double mu = event.mu.dbl();
        if (!std::isfinite(output_energy) || output_energy < 0.0
            || !std::isfinite(mu) || mu < -1.0 || mu > 1.0) {
          NCPLUGIN_MSG("Sampling check failed for " << test_case.material
                       << " at E=" << energy << " eV: eout="
                       << output_energy << ", mu=" << mu);
          nc_assert_always(false);
        }
      }
    }
  }
}
