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

#include "NCPhysicsModel.hh"

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/freegas/NCFreeGas.hh"
#include "NCrystal/internal/sabscatter/NCSABScatter.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"

#include <algorithm>
#include <cmath>
#include <complex>

namespace {

  struct DirectGASpectrum {
    std::vector<double> energy_meV;
    std::vector<double> density_meV_inv; // f(E) in meV^-1 per atom
  };

  struct DirectGAGrids {
    std::vector<double> alpha;
    std::vector<double> beta;
  };

  // Trapezoidal integration weights for a nonuniform grid.
  static std::vector<double> trapezoidWeights(const std::vector<double>& x)
  {
    std::vector<double> w(x.size(), 0.0);
    if (x.size() < 2)
      return w;
    w.front() = 0.5 * (x[1] - x[0]);
    w.back() = 0.5 * (x.back() - x[x.size()-2]);
    for (size_t i = 1; i + 1 < x.size(); ++i)
      w[i] = 0.5 * (x[i+1] - x[i-1]);
    return w;
  }

  static double trapz(const std::vector<double>& y,
                      const std::vector<double>& x)
  {
    if (x.size() != y.size() || x.size() < 2)
      return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i + 1 < x.size(); ++i)
      sum += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]);
    return sum;
  }

  // Construct logarithmic alpha and symmetric beta grids for S(alpha,beta).
  static DirectGAGrids createDirectGAGrids(double beta_max,
                                            size_t grid_points = 400)
  {
    if (grid_points < 2)
      NCRYSTAL_THROW(BadInput,
                     "direct_ga requires at least two alpha and beta grid points");
    if (!std::isfinite(beta_max) || !(beta_max > 1.0e-5))
      NCRYSTAL_THROW(BadInput,
                     "direct_ga requires beta_max greater than 1e-5");

    DirectGAGrids grids;
    const double alpha_max = 4.0 * beta_max;
    if (!(alpha_max > 1.0e-4))
      NCRYSTAL_THROW(BadInput,
                     "direct_ga requires alpha_max greater than 1e-4");
    std::vector<double> positive_beta(grid_points);

    double log_start = std::log10(1.0e-5);
    double log_end = std::log10(beta_max);
    double log_step = (log_end - log_start) / double(grid_points - 1);
    for (size_t i = 0; i < grid_points; ++i)
      positive_beta[i] = std::pow(10.0, log_start + double(i) * log_step);

    grids.alpha.resize(grid_points);
    log_start = std::log10(1.0e-4);
    log_end = std::log10(alpha_max);
    log_step = (log_end - log_start) / double(grid_points - 1);
    for (size_t i = 0; i < grid_points; ++i)
      grids.alpha[i] = std::pow(10.0, log_start + double(i) * log_step);

    grids.beta.reserve(2 * grid_points + 1);
    for (size_t i = 0; i < grid_points; ++i)
      grids.beta.push_back(-positive_beta[grid_points - 1 - i]);
    grids.beta.push_back(0.0);
    grids.beta.insert(grids.beta.end(), positive_beta.begin(), positive_beta.end());
    return grids;
  }

  // Use the absolute liquid spectrum as supplied.
  static DirectGASpectrum liquidSpectrum(const NCP::DirectGAInfo& input)
  {
    DirectGASpectrum spectrum;
    spectrum.energy_meV.assign(input.energy_meV.begin(), input.energy_meV.end());
    spectrum.density_meV_inv.assign(input.density_meV_inv.begin(),
                                    input.density_meV_inv.end());
    const double area = trapz(spectrum.density_meV_inv,
                              spectrum.energy_meV);
    if (!(area > 0.0) || !std::isfinite(area))
      NCRYSTAL_THROW(BadInput,
                     "direct_ga requires a positive absolute frequency-spectrum area");
    return spectrum;
  }

  // Convert NCrystal's normalized solid VDOS to the f(E) grid used by DirectGA.
  static DirectGASpectrum solidSpectrum(const NC::DI_VDOS& di_vdos)
  {
    NC::VDOSEval vdos_eval(di_vdos.vdosData());
    const auto grid = vdos_eval.getGridInfo();
    if (grid.npts_extended < 2 || !(grid.emax > 0.0))
      NCRYSTAL_THROW(BadInput,
                     "solid direct_ga requires a non-empty positive VDOS grid");

    DirectGASpectrum spectrum;
    spectrum.energy_meV.resize(grid.npts_extended);
    spectrum.density_meV_inv.resize(grid.npts_extended);
    const double step_eV = grid.emax / double(grid.npts_extended - 1);
    for (size_t i = 0; i < grid.npts_extended; ++i) {
      const double energy_eV = double(i) * step_eV;
      spectrum.energy_meV[i] = energy_eV * 1000.0;
      spectrum.density_meV_inv[i] =
          i == 0 ? 0.0 : vdos_eval.eval(energy_eV) / 1000.0;
    }

    const double area = trapz(spectrum.density_meV_inv, spectrum.energy_meV);
    if (!(area > 0.0) || !std::isfinite(area))
      NCRYSTAL_THROW(BadInput,
                     "solid direct_ga could not normalize the @DYNINFO VDOS");
    for (double& density : spectrum.density_meV_inv)
      density /= area;
    return spectrum;
  }

  // Construct the complex Gaussian width Gamma(t) from f(E).
  // N. Morishima, "Numerical methods for calculation of low-energy neutron
  // scattering cross sections," Journal of Neutron Research 13 (2005)
  // 225-240, Appendix A, Eqs. (43)-(44).
  // Gamma(t) here is one half of w(t) in Morishima Eq. (44), so Morishima
  // Eq. (43) becomes F_s(Q,t) = exp[-Q^2*Gamma(t)].
  static std::vector<std::complex<double>> buildGamma(
      const std::vector<double>& time_ps,
      double atomic_cv0_a2_ps2,
      const DirectGASpectrum& spectrum,
      double temperature_K)
  {
    const double hbar_meV_ps = (NC::constant_planck / NC::k2Pi) * 1.0e15;
    const double kt_meV = NC::constant_boltzmann * temperature_K * 1000.0;
    const double prefactor =
        (atomic_cv0_a2_ps2 / 3.0) * hbar_meV_ps * hbar_meV_ps;
    const auto energy_weights = trapezoidWeights(spectrum.energy_meV);

    const size_t n_energy = spectrum.energy_meV.size();
    std::vector<double> omega(n_energy);
    std::vector<double> x_coth_x(n_energy);
    for (size_t ie = 0; ie < n_energy; ++ie) {
      const double energy = spectrum.energy_meV[ie];
      const double x = energy / (2.0 * kt_meV);
      omega[ie] = energy / hbar_meV_ps;
      x_coth_x[ie] = std::abs(x) < 1.0e-8 ? 1.0 : x / std::tanh(x);
    }

    std::vector<std::complex<double>> gamma(time_ps.size());
    for (size_t it = 0; it < time_ps.size(); ++it) {
      const double time = time_ps[it];
      double real_integral = 0.0;
      double imaginary_integral = 0.0;
      for (size_t ie = 0; ie < n_energy; ++ie) {
        const double energy = spectrum.energy_meV[ie];
        const double phase = time * omega[ie];
        double real_kernel;
        double imaginary_kernel;
        // Use the exact E -> 0 limits of the real and imaginary kernels.
        if (std::abs(energy) < 1.0e-14) {
          const double time_over_hbar = time / hbar_meV_ps;
          real_kernel = 0.5 * time_over_hbar * time_over_hbar;
          imaginary_kernel = time / (2.0 * kt_meV * hbar_meV_ps);
        } else {
          const double sine_half = std::sin(0.5 * phase);
          real_kernel = 2.0 * sine_half * sine_half / (energy * energy);
          imaginary_kernel = std::sin(phase) / (2.0 * kt_meV * energy);
        }
        const double weighted_density =
            energy_weights[ie] * spectrum.density_meV_inv[ie];
        real_integral += weighted_density * x_coth_x[ie] * real_kernel;
        imaginary_integral += weighted_density * imaginary_kernel;
      }
      gamma[it] = prefactor
                  * std::complex<double>(real_integral, -imaginary_integral);
    }
    return gamma;
  }

  // Calculate the solid width Gamma(infinity), used to remove the elastic
  // contribution before the Fourier transform.
  // Y. Abe et al., "Calculation scheme for thermal neutron scattering
  // cross-sections of hydrogenous molecules in liquid and solid phases by
  // molecular dynamics," Annals of Nuclear Energy 194 (2023) 110101,
  // Eqs. (2)-(5).
  static double buildGammaInfinity(double atomic_cv0_a2_ps2,
                                   const DirectGASpectrum& spectrum,
                                   double temperature_K)
  {
    const double hbar_meV_ps = (NC::constant_planck / NC::k2Pi) * 1.0e15;
    const double kt_meV = NC::constant_boltzmann * temperature_K * 1000.0;
    const double prefactor =
        (atomic_cv0_a2_ps2 / 3.0) * hbar_meV_ps * hbar_meV_ps;
    const auto energy_weights = trapezoidWeights(spectrum.energy_meV);
    const double first_energy = spectrum.energy_meV.at(1);
    // Estimate the finite Debye limit f(E)/E^2 from the first positive point.
    const double zero_limit =
        spectrum.density_meV_inv.at(1) / (first_energy * first_energy);

    double integral = 0.0;
    for (size_t ie = 0; ie < spectrum.energy_meV.size(); ++ie) {
      const double energy = spectrum.energy_meV[ie];
      const double x = energy / (2.0 * kt_meV);
      const double x_coth_x =
          std::abs(x) < 1.0e-8 ? 1.0 : x / std::tanh(x);
      const double density_over_e2 =
          std::abs(energy) < 1.0e-14
              ? zero_limit
              : spectrum.density_meV_inv[ie] / (energy * energy);
      integral += energy_weights[ie] * density_over_e2 * x_coth_x;
    }
    return prefactor * integral;
  }

  // Evaluate exp(value)-1 accurately when value is close to zero.
  static std::complex<double> complexExpm1(
      const std::complex<double>& value)
  {
    if (std::abs(value) > 1.0e-5)
      return std::exp(value) - std::complex<double>(1.0, 0.0);
    return value * (1.0 + value * (0.5 + value / 6.0));
  }

  static NC::SABData buildDirectGASAB(const NCP::DirectGAInfo& input,
                                      const NC::DI_VDOS* di_vdos,
                                      const NC::AtomData& atom_data,
                                      NC::Temperature temperature)
  {
    const double temperature_K = temperature.dbl();
    const double hbar_eV_ps = NC::constant_planck / NC::k2Pi * 1.0e12;
    const double hbar_meV_ps = hbar_eV_ps * 1000.0;
    const double kt_meV = NC::constant_boltzmann * temperature_K * 1000.0;
    const double kt_eV = NC::constant_boltzmann * temperature_K;
    const double neutron_mass =
        (NC::const_neutron_mass_amu * NC::constant_dalton2eVc2)
        / ((NC::constant_c * 1.0e-12) * (NC::constant_c * 1.0e-12));
    const double recoil_eV_a2 =
        hbar_eV_ps * hbar_eV_ps / (2.0 * neutron_mass);

    DirectGASpectrum spectrum;
    if (input.solid) {
      if (!di_vdos)
        NCRYSTAL_THROW(BadInput,
                       "solid direct_ga requires a matching @DYNINFO VDOS");
      spectrum = solidSpectrum(*di_vdos);
    } else {
      spectrum = liquidSpectrum(input);
    }

    const double atom_mass =
        atom_data.averageMassAMU().dbl() * NC::constant_dalton2eVc2
        / ((NC::constant_c * 1.0e-12) * (NC::constant_c * 1.0e-12));
    const double atomic_cv0 =
        3.0 * NC::constant_boltzmann * temperature_K / atom_mass;
    double diffusion_a2_ps = 0.0;
    if (!input.solid) {
      // Compute D from f(0): D = (pi/6) * C_v(0) * hbar * f(0).
      diffusion_a2_ps = (NC::kPi / 6.0) * atomic_cv0
                        * hbar_meV_ps * spectrum.density_meV_inv.front();
      if (!(diffusion_a2_ps > 0.0) || !std::isfinite(diffusion_a2_ps))
        NCRYSTAL_THROW(BadInput,
                       "liquid direct_ga requires a finite positive f(0)");
    }

    // Auto mode uses eight points across the shortest required period.
    const double required_energy_meV =
        std::max(input.max_energy_transfer_meV, spectrum.energy_meV.back());
    const double shortest_period_ps =
        2.0 * NC::kPi * hbar_meV_ps / required_energy_meV;
    const double maximum_time_step_ps = shortest_period_ps / 2.0;
    const double time_step_ps = input.time_step_auto
        ? shortest_period_ps / 8.0
        : input.time_step_ps;
    if (!(time_step_ps > 0.0) || !std::isfinite(time_step_ps))
      NCRYSTAL_THROW(BadInput,
                     "direct_ga time_step_ps must be finite and positive");
    if (time_step_ps > maximum_time_step_ps * (1.0 + 1.0e-12))
      NCRYSTAL_THROW2(BadInput,
                      "direct_ga time_step_ps=" << time_step_ps
                      << " ps is incompatible with the required energy range "
                      << required_energy_meV << " meV; it must not exceed "
                      << maximum_time_step_ps << " ps");
    if (input.numerical_time_end_ps < time_step_ps)
      NCRYSTAL_THROW2(BadInput,
                      "direct_ga numerical_time_end_ps="
                      << input.numerical_time_end_ps
                      << " ps is shorter than the selected time_step_ps="
                      << time_step_ps << " ps");

    // Build the linear time grid for the numerical Fourier integral.
    const size_t final_time_index = std::max<size_t>(
        1, static_cast<size_t>(
            std::ceil(input.numerical_time_end_ps / time_step_ps)));
    std::vector<double> time_ps(final_time_index + 1);
    for (size_t it = 0; it < time_ps.size(); ++it)
      time_ps[it] = double(it) * time_step_ps;

    auto gamma = buildGamma(time_ps, atomic_cv0, spectrum, temperature_K);
    const double gamma_infinity = input.solid
        ? buildGammaInfinity(atomic_cv0, spectrum, temperature_K)
        : 0.0;

    const auto grids = createDirectGAGrids(
        input.max_energy_transfer_meV / kt_meV);
    const auto& alpha_grid = grids.alpha;
    const auto& beta_grid = grids.beta;
    const size_t n_alpha = alpha_grid.size();
    const size_t n_beta = beta_grid.size();

    // Convert NCrystal's alpha coordinate to physical momentum transfer Q.
    std::vector<double> q_grid(n_alpha);
    for (size_t ia = 0; ia < n_alpha; ++ia)
      q_grid[ia] = std::sqrt(alpha_grid[ia] * kt_eV / recoil_eV_a2);

    // Use E_loss = -beta*k_B*T, reversing beta to keep energy loss increasing.
    std::vector<double> loss_grid(n_beta);
    for (size_t iloss = 0; iloss < n_beta; ++iloss)
      loss_grid[iloss] = -beta_grid[n_beta - 1 - iloss] * kt_meV;

    // Precompute exp(-i*E_loss*t/hbar), which is independent of Q.
    const size_t n_time = time_ps.size();
    std::vector<double> omega(n_beta);
    for (size_t iloss = 0; iloss < n_beta; ++iloss)
      omega[iloss] = -loss_grid[iloss] / hbar_meV_ps;

    std::vector<double> phase_real(n_time * n_beta);
    std::vector<double> phase_imag(n_time * n_beta);
    for (size_t it = 0; it < n_time; ++it) {
      const size_t offset = it * n_beta;
      for (size_t iloss = 0; iloss < n_beta; ++iloss) {
        const double phase = omega[iloss] * time_ps[it];
        phase_real[offset + iloss] = std::cos(phase);
        phase_imag[offset + iloss] = std::sin(phase);
      }
    }

    NC::VectD sab(n_alpha * n_beta, 0.0);
    const auto time_weights = trapezoidWeights(time_ps);

    // Evaluate the one-sided Fourier transform of F_s(Q,t).
    // N. Morishima, Journal of Neutron Research 13 (2005) 225-240,
    // Eqs. (3)-(4) and Appendix A, Eq. (43).
    constexpr double negative_sab_numerical_tolerance = 1.0e-4;
    constexpr size_t q_block_size = 8;
    for (size_t ia_start = 0; ia_start < n_alpha;
         ia_start += q_block_size) {
      const size_t ia_stop = std::min(ia_start + q_block_size, n_alpha);
      const size_t n_q = ia_stop - ia_start;
      std::vector<double> q2(n_q);
      for (size_t iq = 0; iq < n_q; ++iq) {
        const double q = q_grid[ia_start + iq];
        q2[iq] = q * q;
      }

      // Accumulate the Fourier sums for all Q values in the current block.
      std::vector<double> integral_real(n_beta * n_q, 0.0);
      std::vector<double> weighted_fs_real(n_q);
      std::vector<double> weighted_fs_imag(n_q);

      for (size_t it = 0; it < n_time; ++it) {
        for (size_t iq = 0; iq < n_q; ++iq) {
          std::complex<double> fs;
          if (input.solid) {
            const double elastic = std::exp(-q2[iq] * gamma_infinity);
            fs = elastic
                 * complexExpm1(-q2[iq] * (gamma[it] - gamma_infinity));
          } else {
            fs = std::exp(-q2[iq] * gamma[it]);
          }
          weighted_fs_real[iq] = time_weights[it] * fs.real();
          weighted_fs_imag[iq] = time_weights[it] * fs.imag();
        }

        const size_t phase_offset = it * n_beta;
        for (size_t iloss = 0; iloss < n_beta; ++iloss) {
          const double real_phase = phase_real[phase_offset + iloss];
          const double imaginary_phase = phase_imag[phase_offset + iloss];
          const size_t integral_offset = iloss * n_q;
          for (size_t iq = 0; iq < n_q; ++iq)
            integral_real[integral_offset + iq] +=
                weighted_fs_real[iq] * real_phase
                - weighted_fs_imag[iq] * imaginary_phase;
        }
      }

      const size_t final_phase_offset = (n_time - 1) * n_beta;
      std::vector<double> row(n_beta);
      for (size_t iq = 0; iq < n_q; ++iq) {
        const size_t ia = ia_start + iq;
        std::complex<double> fs_last;
        double tail_decay = 0.0;
        if (!input.solid) {
          // Continue the liquid correlation analytically beyond the time grid
          // as F_s(Q,T+u) = F_s(Q,T) * exp(-Q^2*D*u).
          // For liquids, the long-time real part of Gamma(t) is linear and its
          // imaginary part is constant: Y. Abe et al., Annals of Nuclear
          // Energy 194 (2023) 110101, Eqs. (23)-(24).
          fs_last = std::exp(-q2[iq] * gamma.back());
          tail_decay = std::max(q2[iq] * diffusion_a2_ps, 1.0e-12);
        }
        for (size_t iloss = 0; iloss < n_beta; ++iloss) {
          double result = integral_real[iloss * n_q + iq];
          if (!input.solid) {
            const std::complex<double> phase_last(
                phase_real[final_phase_offset + iloss],
                phase_imag[final_phase_offset + iloss]);
            const std::complex<double> tail =
                fs_last * phase_last
                / std::complex<double>(tail_decay, -omega[iloss]);
            result += tail.real();
          }
          row[iloss] = result / (NC::kPi * hbar_meV_ps);
        }

        for (size_t ibeta = 0; ibeta < n_beta; ++ibeta) {
          const double value = kt_meV * row[n_beta - 1 - ibeta];
          if (!std::isfinite(value))
            NCRYSTAL_THROW2(BadInput,
                            "direct_ga produced a non-finite S(alpha,beta) "
                            "value at alpha=" << alpha_grid[ia]
                            << ", beta=" << beta_grid[ibeta]);
          if (value < -negative_sab_numerical_tolerance)
            NCRYSTAL_THROW2(BadInput,
                            "direct_ga produced negative S(alpha,beta)="
                            << value << " at alpha=" << alpha_grid[ia]
                            << ", beta=" << beta_grid[ibeta]
                            << "; tolerance is "
                            << negative_sab_numerical_tolerance
                            << ". Try decreasing time_step_ps and/or increasing "
                               "max_energy_transfer_meV.");
          sab[ibeta * n_alpha + ia] = std::max(0.0, value);
        }
      }
    }

    return NC::SABData(NC::VectD(alpha_grid), NC::VectD(beta_grid),
                       std::move(sab), temperature,
                       atom_data.scatteringXS(),
                       atom_data.averageMassAMU(),
                       input.max_energy_transfer_meV * 1.0e-3);
  }

}

bool NCP::PhysicsModel::isApplicable(const NC::Info& info)
{
  return info.countCustomSections("DIRECTGA") > 0;
}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo(const NC::Info& info)
{
  const auto section_count = info.countCustomSections("DIRECTGA");
  if (section_count != 1)
    NCRYSTAL_THROW(BadInput,
                   "Exactly one @CUSTOM_DIRECTGA section is required");

  // Parse one DirectGA settings-and-data block for each element.
  struct SeenFields {
    bool phase = false;
    bool max_energy_transfer = false;
    bool numerical_time_end = false;
    bool time_step = false;
    bool model = false;
    bool yk_model = false;
    bool sq = false;
    bool spectrum = false;
    enum class Table { None, Spectrum, StructureFactor } table = Table::None;
  };

  std::unordered_map<std::string, DirectGAInfo> input_data;
  std::unordered_map<std::string, SeenFields> seen_fields;
  std::string current_label;

  for (const auto& line : info.getCustomSection("DIRECTGA")) {
    if (line.empty())
      continue;
    if (line.at(0) == "element") {
      if (line.size() != 2)
        NCRYSTAL_THROW(BadInput,
                       "@CUSTOM_DIRECTGA element lines must be 'element <label>'");
      current_label = line.at(1);
      if (input_data.find(current_label) != input_data.end())
        NCRYSTAL_THROW2(BadInput,
                        "Duplicate @CUSTOM_DIRECTGA block for element '"
                        << current_label << "'");
      input_data.emplace(current_label, DirectGAInfo{});
      seen_fields.emplace(current_label, SeenFields{});
      continue;
    }

    if (current_label.empty())
      NCRYSTAL_THROW(BadInput,
                     "@CUSTOM_DIRECTGA data must follow an 'element <label>' line");
    auto& input = input_data.at(current_label);
    auto& seen = seen_fields.at(current_label);

    if (line.at(0) == "phase") {
      if (seen.table != SeenFields::Table::None)
        NCRYSTAL_THROW(BadInput,
                       "direct_ga settings must precede table rows");
      if (line.size() != 2 || seen.phase)
        NCRYSTAL_THROW(BadInput,
                       "phase must be specified once as 'liquid' or 'solid'");
      const std::string value = NC::upperCase(line.at(1));
      if (value != "LIQUID" && value != "SOLID")
        NCRYSTAL_THROW(BadInput,"phase must be 'liquid' or 'solid'");
      input.solid = value == "SOLID";
      input.numerical_time_end_ps = input.solid ? 5.0 : 2.0;
      seen.phase = true;
      continue;
    }
    if (!seen.phase)
      NCRYSTAL_THROW2(BadInput,
                      "phase must precede other direct_ga settings for element '"
                      << current_label << "'");

    if (line.at(0) == "model") {
      if (seen.table != SeenFields::Table::None)
        NCRYSTAL_THROW(BadInput,
                       "direct_ga settings must precede table rows");
      if (line.size() != 2 || seen.model)
        NCRYSTAL_THROW(BadInput,
                       "model must be specified at most once as incoherent, "
                       "skold, or vineyard");
      const std::string value = NC::upperCase(line.at(1));
      if (value != "INCOHERENT" && value != "SKOLD"
          && value != "VINEYARD")
        NCRYSTAL_THROW(BadInput,
                       "model must be incoherent, skold, or vineyard");
      input.corrections.coherent_model = value;
      seen.model = true;
      continue;
    }

    if (line.at(0) == "yk_model") {
      if (seen.table != SeenFields::Table::None)
        NCRYSTAL_THROW(BadInput,
                       "direct_ga settings must precede table rows");
      if (line.size() != 2 || seen.yk_model)
        NCRYSTAL_THROW(BadInput,
                       "yk_model must be specified at most once as ortho or para");
      const std::string value = NC::upperCase(line.at(1));
      if (value != "ORTHO" && value != "PARA")
        NCRYSTAL_THROW(BadInput,"yk_model must be ortho or para");
      input.corrections.yk_model = value;
      seen.yk_model = true;
      continue;
    }

    if (line.at(0) == "max_energy_transfer_meV") {
      if (seen.table != SeenFields::Table::None)
        NCRYSTAL_THROW(BadInput,
                       "direct_ga settings must precede table rows");
      double value;
      if (line.size() != 2 || seen.max_energy_transfer
          || !NC::safe_str2dbl(line.at(1), value)
          || !std::isfinite(value) || !(value > 0.0))
        NCRYSTAL_THROW(BadInput,
                       "max_energy_transfer_meV must be a finite positive "
                       "value, specified at most once");
      input.max_energy_transfer_meV = value;
      seen.max_energy_transfer = true;
      continue;
    }

    if (line.at(0) == "numerical_time_end_ps") {
      if (seen.table != SeenFields::Table::None)
        NCRYSTAL_THROW(BadInput,
                       "direct_ga settings must precede table rows");
      double value;
      if (line.size() != 2 || seen.numerical_time_end
          || !NC::safe_str2dbl(line.at(1), value)
          || !std::isfinite(value) || !(value > 0.0))
        NCRYSTAL_THROW(BadInput,
                       "numerical_time_end_ps must be a finite positive "
                       "value, specified at most once");
      input.numerical_time_end_ps = value;
      seen.numerical_time_end = true;
      continue;
    }

    if (line.at(0) == "time_step_ps") {
      if (seen.table != SeenFields::Table::None)
        NCRYSTAL_THROW(BadInput,
                       "direct_ga settings must precede table rows");
      if (line.size() != 2 || seen.time_step)
        NCRYSTAL_THROW(BadInput,
                       "time_step_ps must be 'auto' or a finite positive "
                       "value, specified at most once");
      const std::string value_text = NC::upperCase(line.at(1));
      if (value_text == "AUTO") {
        input.time_step_auto = true;
      } else {
        double value;
        if (!NC::safe_str2dbl(line.at(1), value)
            || !std::isfinite(value) || !(value > 0.0))
          NCRYSTAL_THROW(BadInput,
                         "time_step_ps must be 'auto' or a finite positive value");
        input.time_step_auto = false;
        input.time_step_ps = value;
      }
      seen.time_step = true;
      continue;
    }

    // A liquid block may contain an S(Q) table followed by its f(E) table.
    if (line.at(0) == "sq") {
      if (line.size() != 1 || input.solid || seen.sq || seen.spectrum)
        NCRYSTAL_THROW(BadInput,
                       "sq must occur once before spectrum in a liquid element block");
      seen.table = SeenFields::Table::StructureFactor;
      seen.sq = true;
      continue;
    }

    if (line.at(0) == "spectrum") {
      if (line.size() != 1 || input.solid || seen.spectrum)
        NCRYSTAL_THROW(BadInput,
                       "spectrum must occur once in a liquid element block");
      seen.table = SeenFields::Table::Spectrum;
      seen.spectrum = true;
      continue;
    }

    if (input.solid)
      NCRYSTAL_THROW2(BadInput,
                      "solid direct_ga uses its matching @DYNINFO VDOS; custom "
                      "spectrum rows are not allowed for element '"
                      << current_label << "'");
    if (seen.table == SeenFields::Table::None)
      NCRYSTAL_THROW2(BadInput,
                      "Expected sq or spectrum before table rows for element '"
                      << current_label << "'");

    double first_value;
    double second_value;
    if (line.size() != 2
        || !NC::safe_str2dbl(line.at(0), first_value)
        || !NC::safe_str2dbl(line.at(1), second_value)
        || !std::isfinite(first_value) || !std::isfinite(second_value)
        || first_value < 0.0 || second_value < 0.0)
      NCRYSTAL_THROW2(BadInput,
                      "Invalid table row for element '" << current_label
                      << "'; expected two non-negative finite values");
    if (seen.table == SeenFields::Table::StructureFactor) {
      if (!input.corrections.q_values.empty()
          && !(first_value > input.corrections.q_values.back()))
        NCRYSTAL_THROW2(BadInput,
                        "Structure-factor Q values for element '"
                        << current_label << "' must be strictly increasing");
      input.corrections.q_values.push_back(first_value);
      input.corrections.structure_factor.push_back(second_value);
    } else {
      if (!input.energy_meV.empty() && !(first_value > input.energy_meV.back()))
        NCRYSTAL_THROW2(BadInput,
                        "Spectrum energies for element '" << current_label
                        << "' must be strictly increasing");
      input.energy_meV.push_back(first_value);
      input.density_meV_inv.push_back(second_value);
    }
  }

  if (input_data.empty())
    NCRYSTAL_THROW(BadInput,
                   "@CUSTOM_DIRECTGA must contain at least one element block");

  bool have_solid = false;
  bool have_liquid = false;
  // Validate each element block against its phase and dynamic information.
  for (const auto& entry : input_data) {
    const auto& label = entry.first;
    const auto& input = entry.second;
    const auto& seen = seen_fields.at(label);
    if (!seen.phase)
      NCRYSTAL_THROW2(BadInput,
                      "Incomplete @CUSTOM_DIRECTGA block for element '"
                      << label << "'");
    have_solid = have_solid || input.solid;
    have_liquid = have_liquid || !input.solid;

    bool matching_dynamic_info = false;
    bool matching_vdos = false;
    for (const auto& di : info.getDynamicInfoList()) {
      if (info.displayLabel(di->atom().index) != label)
        continue;
      matching_dynamic_info = true;
      matching_vdos = matching_vdos
                      || dynamic_cast<const NC::DI_VDOS*>(di.get());
    }
    if (!matching_dynamic_info)
      NCRYSTAL_THROW2(BadInput,
                      "@CUSTOM_DIRECTGA has no matching @DYNINFO entry for "
                      "element '" << label << "'");

    if (input.solid) {
      if (usesYKSkoldCorrection(input.corrections)
          || !input.corrections.q_values.empty())
        NCRYSTAL_THROW2(BadInput,
                        "Young-Koppel, Skold, and Vineyard corrections are "
                        "only supported for liquid direct_ga elements; element '"
                        << label << "' is solid");
      if (!matching_vdos)
        NCRYSTAL_THROW2(BadInput,
                        "solid direct_ga requires a matching @DYNINFO VDOS for "
                        "element '" << label << "'");
    } else {
      const bool coherent_model =
          input.corrections.coherent_model == "SKOLD"
          || input.corrections.coherent_model == "VINEYARD";
      if (coherent_model && input.corrections.q_values.size() < 2)
        NCRYSTAL_THROW2(BadInput,
                        "The " << input.corrections.coherent_model
                        << " model requires at least two S(Q) rows for element '"
                        << label << "'");
      if (!coherent_model && !input.corrections.q_values.empty())
        NCRYSTAL_THROW2(BadInput,
                        "S(Q) rows require model skold or vineyard for element '"
                        << label << "'");
      if (!seen.spectrum)
        NCRYSTAL_THROW2(BadInput,
                        "liquid direct_ga requires a spectrum marker for element '"
                        << label << "'");
      if (input.energy_meV.size() < 2)
        NCRYSTAL_THROW2(BadInput,
                        "liquid direct_ga requires at least two spectrum rows "
                        "for element '" << label << "'");
      if (std::abs(input.energy_meV.front()) > 1.0e-14)
        NCRYSTAL_THROW2(BadInput,
                        "The direct_ga spectrum for element '" << label
                        << "' must start with an explicit E=0 point");
      if (!(input.density_meV_inv.front() > 0.0))
        NCRYSTAL_THROW2(BadInput,
                        "The direct_ga spectrum for element '" << label
                        << "' requires a positive f(0)");
    }
  }

  // Elastic scattering is handled differently for liquid and solid phases.
  if (have_solid && have_liquid)
    NCRYSTAL_THROW(BadInput,
                   "Solid and liquid direct_ga elements cannot be mixed in "
                   "one material");

  if (have_solid) {
    // Every solid VDOS element must be handled by DirectGA.
    for (const auto& di : info.getDynamicInfoList()) {
      if (!dynamic_cast<const NC::DI_VDOS*>(di.get()))
        continue;
      const auto label = info.displayLabel(di->atom().index);
      const auto it = input_data.find(label);
      if (it == input_data.end() || !it->second.solid)
        NCRYSTAL_THROW2(BadInput,
                        "solid direct_ga requires a @CUSTOM_DIRECTGA block for "
                        "every @DYNINFO VDOS element; missing element '"
                        << label << "'");
    }
  }

  return PhysicsModel(std::move(input_data), info);
}

NCP::PhysicsModel::PhysicsModel(
    std::unordered_map<std::string, DirectGAInfo> input_data,
    const NC::Info& info)
  : m_direct_ga_data(std::move(input_data))
{
  NC::ProcImpl::ProcComposition::ComponentList components;
  // Build one scattering component for each @DYNINFO entry.
  for (const auto& di : info.getDynamicInfoList()) {
    const auto label = info.displayLabel(di->atom().index);
    const auto input_it = m_direct_ga_data.find(label);
    if (input_it == m_direct_ga_data.end()) {
      // Use free-gas scattering for elements without DirectGA data.
      components.push_back({di->fraction(),
                            NC::makeSO<NC::FreeGas>(di->temperature(),
                                                    di->atomData())});
      continue;
    }

    const auto* di_vdos = dynamic_cast<const NC::DI_VDOS*>(di.get());
    auto sab = buildDirectGASAB(input_it->second, di_vdos,
                                di->atomData(), di->temperature());
    if (usesYKSkoldCorrection(input_it->second.corrections))
      sab = applyYKSkold(sab, di->atomData(),
                         input_it->second.corrections);
    components.push_back({di->fraction(),
                          NC::makeSO<NC::SABScatter>(std::move(sab))});
  }

  auto process = NC::ProcImpl::ProcComposition::consumeAndCombine(
      std::move(components), NC::ProcessType::Scatter);
  m_proc = process.optional();
}

bool NCP::PhysicsModel::usesSolid() const noexcept
{
  for (const auto& entry : m_direct_ga_data)
    if (entry.second.solid)
      return true;
  return false;
}

double NCP::PhysicsModel::calcCrossSection(NC::CachePtr& cache,
                                           double neutron_ekin) const
{
  return m_proc->crossSectionIsotropic(cache,
      NC::NeutronEnergy{neutron_ekin}).dbl();
}

NCP::PhysicsModel::ScatEvent
NCP::PhysicsModel::sampleScatteringEvent(NC::CachePtr& cache,
                                         NC::RNG& rng,
                                         double neutron_ekin) const
{
  const auto result = m_proc->sampleScatterIsotropic(
      cache, rng, NC::NeutronEnergy{neutron_ekin});
  return {result.ekin.dbl(), result.mu.dbl()};
}
