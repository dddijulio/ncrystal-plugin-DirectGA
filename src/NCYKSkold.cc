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

#include "NCYKSkold.hh"

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/utils/NCMath.hh"

#include <algorithm>
#include <cmath>

namespace {

  // Linearly interpolate S(Q), holding its endpoint values outside the table.
  double interpolateStructureFactor(const NC::VectD& grid,
                                    const NC::VectD& values,
                                    double x)
  {
    if (x <= grid.front())
      return values.front();
    if (x >= grid.back())
      return values.back();
    const auto upper = std::upper_bound(grid.begin(), grid.end(), x);
    const size_t hi = static_cast<size_t>(upper - grid.begin());
    const size_t lo = hi - 1;
    const double fraction = (x - grid[lo]) / (grid[hi] - grid[lo]);
    return (1.0 - fraction) * values[lo] + fraction * values[hi];
  }

  // Evaluate spherical Bessel j_n(x) using a series or stable recurrence.
  double sphericalBessel(int n, double x)
  {
    if (n < 0 || x < 0.0)
      return 0.0;
    if (x == 0.0)
      return n == 0 ? 1.0 : 0.0;

    if (x < 0.1) {
      double leading = 1.0;
      for (int k = 1; k <= n; ++k)
        leading *= x / double(2 * k + 1);

      double term = 1.0;
      double series = 1.0;
      const double y = -0.25 * x * x;
      for (int k = 1; k < 32; ++k) {
        term *= y / (double(k) * (double(n + k) + 0.5));
        series += term;
        if (std::abs(term) <= 1.0e-16 * std::max(1.0, std::abs(series)))
          break;
      }
      return leading * series;
    }

    const double j0 = std::sin(x) / x;
    if (n == 0)
      return j0;
    const double j1 = std::sin(x) / (x * x) - std::cos(x) / x;
    if (n == 1)
      return j1;

    if (x > n) {
      double jm1 = j0;
      double j = j1;
      for (int ell = 1; ell < n; ++ell) {
        const double jp1 = (2 * ell + 1) * j / x - jm1;
        jm1 = j;
        j = jp1;
      }
      return j;
    }

    const int m = std::max(n + 32, int(x) + 32);
    double jp1 = 0.0;
    double j = 1.0;
    double target = 0.0;
    for (int ell = m; ell > 0; --ell) {
      const double jm1 = (2 * ell + 1) * j / x - jp1;
      if (ell - 1 == n)
        target = jm1;
      jp1 = j;
      j = jm1;
      if (std::abs(j) > 1.0e100) {
        jp1 *= 1.0e-100;
        j *= 1.0e-100;
        target *= 1.0e-100;
      }
    }
    return j0 * target / j;
  }

  // Return |<j1 0, j2 0 | l 0>|^2.
  double cg00Squared(int j1, int j2, int l)
  {
    if (l < std::abs(j1-j2) || l > j1+j2)
      return 0.0;
    if ((j1+j2+l) % 2 != 0)
      return 0.0;
    const int s = (j1+j2+l) / 2;
    const double log_c2 =
        std::log(2*l+1)
      + std::lgamma(j1+j2-l+1)
      + std::lgamma(j1-j2+l+1)
      + std::lgamma(-j1+j2+l+1)
      + 2.0*std::lgamma(s+1)
      - std::lgamma(j1+j2+l+2)
      - 2.0*std::lgamma(s-j1+1)
      - 2.0*std::lgamma(s-j2+1)
      - 2.0*std::lgamma(s-l+1);
    return std::exp(log_c2);
  }

  std::vector<double> coherentCorrectionDelta(
      const NC::SABData& self_sab,
      const NCP::YKSkoldInfo& input,
      int alpha_index,
      double alpha,
      double q)
  {
    // Vineyard approximation: G. H. Vineyard, "Scattering of Slow Neutrons
    // by a Liquid," Physical Review 110 (1958) 999-1010, Eq. (14).
    // Skold modification: K. Skold, "Small Energy Transfer Scattering of
    // Cold Neutrons from Liquid Argon," Physical Review Letters 19 (1967)
    // 1023-1025, Eq. (2).
    const auto& alpha_grid = self_sab.alphaGrid();
    const auto& beta_grid = self_sab.betaGrid();
    const auto& sab = self_sab.sab();
    const int n_alpha = static_cast<int>(alpha_grid.size());
    const int n_beta = static_cast<int>(beta_grid.size());

    const double sq = interpolateStructureFactor(
        input.q_values, input.structure_factor, q);
    const double sq_effective = std::max(sq, 1.0e-8);
    const double shifted_alpha = input.coherent_model == "SKOLD"
        ? alpha / sq_effective
        : alpha;

    auto interpolate_self = [&](double alpha_value, int beta_index) {
      if (alpha_value <= alpha_grid.front())
        return sab[beta_index * n_alpha];
      if (alpha_value >= alpha_grid.back())
        return sab[beta_index * n_alpha + n_alpha - 1];
      const auto upper = std::upper_bound(
          alpha_grid.begin(), alpha_grid.end(), alpha_value);
      const int hi = static_cast<int>(upper - alpha_grid.begin());
      const int lo = hi - 1;
      const double fraction =
          (alpha_value - alpha_grid[lo]) / (alpha_grid[hi] - alpha_grid[lo]);
      return (1.0 - fraction) * sab[beta_index * n_alpha + lo]
             + fraction * sab[beta_index * n_alpha + hi];
    };

    std::vector<double> result(n_beta);
    for (int ibeta = 0; ibeta < n_beta; ++ibeta) {
      const double self_alpha = sab[ibeta * n_alpha + alpha_index];
      const double self_shifted = interpolate_self(shifted_alpha, ibeta);
      result[ibeta] = sq * self_shifted - self_alpha;
    }
    return result;
  }

  NC::SABData applyCoherentCorrection(const NC::SABData& self_sab,
                                      const NC::AtomData& atom_data,
                                      const NCP::YKSkoldInfo& input)
  {
    const auto& alpha_grid = self_sab.alphaGrid();
    const auto& beta_grid = self_sab.betaGrid();
    const int n_alpha = static_cast<int>(alpha_grid.size());
    const int n_beta = static_cast<int>(beta_grid.size());
    NC::VectD corrected(self_sab.sab().begin(), self_sab.sab().end());

    const double hbar = NC::constant_planck / NC::k2Pi * 1.0e12;
    const double neutron_mass =
        (NC::const_neutron_mass_amu * NC::constant_dalton2eVc2)
        / ((NC::constant_c * 1.0e-12) * (NC::constant_c * 1.0e-12));
    const double coherent_xs = atom_data.scatteringXS().dbl()
                               - atom_data.incoherentXS().dbl();
    const double coherent_fraction =
        coherent_xs / atom_data.scatteringXS().dbl();

    for (int ialpha = 0; ialpha < n_alpha; ++ialpha) {
      const double alpha = alpha_grid[ialpha];
      const double q = std::sqrt(
          2.0 * neutron_mass * self_sab.temperature().kT() * alpha) / hbar;
      const auto bracket = coherentCorrectionDelta(
          self_sab, input, ialpha, alpha, q);
      for (int ibeta = 0; ibeta < n_beta; ++ibeta) {
        const size_t index = ibeta * n_alpha + ialpha;
        corrected[index] += coherent_fraction * bracket[ibeta];
        corrected[index] = std::max(0.0, corrected[index]);
      }
    }

    return NC::SABData(NC::VectD(alpha_grid), NC::VectD(beta_grid),
                       std::move(corrected), self_sab.temperature(),
                       self_sab.boundXS(), self_sab.elementMassAMU(),
                       self_sab.suggestedEmax());
  }

  NC::SABData applyYoungKoppel(const NC::SABData& self_sab,
                               const NC::AtomData& atom_data,
                               const NCP::YKSkoldInfo& input)
  {
    // J. A. Young and J. U. Koppel, "Slow Neutron Scattering by Molecular
    // Hydrogen and Deuterium," Physical Review 135 (1964) A603-A611.
    const double atomic_mass_amu = atom_data.averageMassAMU().dbl();
    if (atomic_mass_amu > 2.5)
      NCRYSTAL_THROW2(BadInput,
                      "yk_model is only supported for H2 or D2, got atom mass "
                      << atomic_mass_amu << " amu");
    const bool is_h2 = atomic_mass_amu < 1.5;
    const bool is_ortho = input.yk_model == "ORTHO";

    const double kt = NC::constant_boltzmann * self_sab.temperature().dbl();
    const double rotational_constant = is_h2 ? 7.35e-3 : 3.71e-3;
    const double atom_to_cm = 0.3707;
    constexpr int j_max = 10;
    const double hbar = NC::constant_planck / NC::k2Pi * 1.0e12;
    const double neutron_mass =
        (NC::const_neutron_mass_amu * NC::constant_dalton2eVc2)
        / ((NC::constant_c * 1.0e-12) * (NC::constant_c * 1.0e-12));

    // Boltzmann populations within the selected ortho or para species.
    // Ortho-H2 and para-D2 have odd J; para-H2 and ortho-D2 have even J.
    const int j_start = ((is_h2 && is_ortho) || (!is_h2 && !is_ortho))
        ? 1 : 0;
    std::vector<double> populations(j_max + 1, 0.0);
    double partition = 0.0;
    for (int j = j_start; j <= j_max; j += 2) {
      populations[j] = (2*j + 1)
                       * std::exp(-j*(j+1)*rotational_constant / kt);
      partition += populations[j];
    }
    for (int j = j_start; j <= j_max; j += 2)
      populations[j] /= partition;

    const double coherent_xs = atom_data.scatteringXS().dbl()
                               - atom_data.incoherentXS().dbl();
    const double incoherent_xs = atom_data.incoherentXS().dbl();
    const double bc2 = coherent_xs / (2.0 * NC::k2Pi);
    const double bi2 = incoherent_xs / (2.0 * NC::k2Pi);

    // Spin-correlation weights from Young and Koppel (1964), Sec. II.1
    // for H2 and Appendix A for D2.
    double even_coefficient;
    double odd_coefficient;
    if (is_h2 && !is_ortho) {
      even_coefficient = bc2;
      odd_coefficient = bi2;
    } else if (is_h2 && is_ortho) {
      even_coefficient = bi2 / 3.0;
      odd_coefficient = bc2 + (2.0 / 3.0) * bi2;
    } else if (!is_h2 && !is_ortho) {
      even_coefficient = (3.0 / 4.0) * bi2;
      odd_coefficient = bc2 + (1.0 / 4.0) * bi2;
    } else {
      even_coefficient = bc2 + (5.0 / 8.0) * bi2;
      odd_coefficient = (3.0 / 8.0) * bi2;
    }

    // Precompute the angular-momentum coupling weights.
    constexpr int l_max = 2 * j_max;
    std::vector<std::vector<std::vector<double>>> cg_squared(
        j_max + 1,
        std::vector<std::vector<double>>(
            j_max + 1, std::vector<double>(l_max + 1, 0.0)));
    for (int j = 0; j <= j_max; ++j)
      for (int jp = 0; jp <= j_max; ++jp)
        for (int l = std::abs(j-jp); l <= j+jp; l += 2)
          cg_squared[j][jp][l] = cg00Squared(j, jp, l);

    const auto& alpha_grid = self_sab.alphaGrid();
    const int n_alpha = static_cast<int>(alpha_grid.size());
    // Calculate the momentum-dependent strength of each rotational transition.
    std::vector<std::vector<std::vector<double>>> momentum_factor(
        j_max + 1,
        std::vector<std::vector<double>>(
            j_max + 1, std::vector<double>(n_alpha, 0.0)));
    for (int ialpha = 0; ialpha < n_alpha; ++ialpha) {
      const double y = atom_to_cm * std::sqrt(
          2.0 * neutron_mass * kt * alpha_grid[ialpha]) / hbar;
      for (int j = 0; j <= j_max; ++j)
        for (int jp = 0; jp <= j_max; ++jp)
          for (int l = std::abs(j-jp); l <= j+jp; l += 2) {
            const double jl = sphericalBessel(l, y);
            momentum_factor[j][jp][ialpha] +=
                4.0 * jl * jl * cg_squared[j][jp][l];
          }
    }

    // Calculate the rotational energy shifts in units of kT.
    std::vector<std::vector<double>> beta_shift(
        j_max + 1, std::vector<double>(j_max + 1, 0.0));
    for (int j = 0; j <= j_max; ++j)
      for (int jp = 0; jp <= j_max; ++jp)
        beta_shift[j][jp] =
            (jp*(jp+1) - j*(j+1)) * rotational_constant / kt;

    const auto& beta_grid = self_sab.betaGrid();
    const auto& self_values = self_sab.sab();
    const int n_beta = static_cast<int>(beta_grid.size());
    // Interpolate the self kernel at the shifted beta.
    auto interpolate_self = [&](int ialpha, double beta_target) {
      if (beta_target <= beta_grid.front()
          || beta_target >= beta_grid.back())
        return 0.0;
      const auto upper = std::upper_bound(
          beta_grid.begin(), beta_grid.end(), beta_target);
      const int hi = static_cast<int>(upper - beta_grid.begin());
      const int lo = hi - 1;
      const double fraction =
          (beta_target - beta_grid[lo]) / (beta_grid[hi] - beta_grid[lo]);
      return (1.0 - fraction) * self_values[lo*n_alpha + ialpha]
             + fraction * self_values[hi*n_alpha + ialpha];
    };

    // Sum the rotational transitions, converting molecular scattering to per atom.
    const double scattering_length_sum = bc2 + bi2;
    NC::VectD yk_sab(n_alpha * n_beta, 0.0);
    for (int ialpha = 0; ialpha < n_alpha; ++ialpha) {
      for (int ibeta = 0; ibeta < n_beta; ++ibeta) {
        double value = 0.0;
        for (int j = j_start; j <= j_max; j += 2) {
          for (int jp = 0; jp <= j_max; ++jp) {
            const double coefficient = jp % 2 == 0
                ? even_coefficient : odd_coefficient;
            const double transition_weight = coefficient * (2*jp + 1);
            const double shifted_self = interpolate_self(
                ialpha, beta_grid[ibeta] + beta_shift[j][jp]);
            value += 0.5 * populations[j] * transition_weight
                     * shifted_self * momentum_factor[j][jp][ialpha];
          }
        }
        yk_sab[ibeta * n_alpha + ialpha] =
            value / scattering_length_sum;
      }
    }

    // Add the intermolecular coherent contribution after the intramolecular
    // Young-Koppel sum (IAEA International Evaluation Co-operation Vol. 42,
    // Eq. 7, p. 20).
    if (!input.q_values.empty()) {
      const double coherent_prefactor = 2.0 * bc2 / scattering_length_sum;
      for (int ialpha = 0; ialpha < n_alpha; ++ialpha) {
        const double alpha = alpha_grid[ialpha];
        const double q = std::sqrt(
            2.0 * neutron_mass * kt * alpha) / hbar;
        const double qb = q * atom_to_cm;
        const double j0_qb = qb > 1.0e-10 ? std::sin(qb) / qb : 1.0;
        const auto bracket = coherentCorrectionDelta(
            self_sab, input, ialpha, alpha, q);
        for (int ibeta = 0; ibeta < n_beta; ++ibeta) {
          const size_t index = ibeta * n_alpha + ialpha;
          yk_sab[index] += coherent_prefactor * j0_qb * j0_qb
                           * bracket[ibeta];
          yk_sab[index] = std::max(0.0, yk_sab[index]);
        }
      }
    }

    return NC::SABData(NC::VectD(alpha_grid), NC::VectD(beta_grid),
                       std::move(yk_sab), self_sab.temperature(),
                       self_sab.boundXS(), self_sab.elementMassAMU(),
                       self_sab.suggestedEmax());
  }

}

bool NCP::usesYKSkoldCorrection(const YKSkoldInfo& input) noexcept
{
  return !input.yk_model.empty()
         || input.coherent_model == "SKOLD"
         || input.coherent_model == "VINEYARD";
}

NC::SABData NCP::applyYKSkold(const NC::SABData& self_sab,
                              const NC::AtomData& atom_data,
                              const YKSkoldInfo& input)
{
  if (!input.yk_model.empty())
    return applyYoungKoppel(self_sab, atom_data, input);
  if (input.coherent_model == "SKOLD"
      || input.coherent_model == "VINEYARD")
    return applyCoherentCorrection(self_sab, atom_data, input);
  NCRYSTAL_THROW(BadInput,
                 "applyYKSkold called without a Young-Koppel, Skold, or "
                 "Vineyard correction");
}
