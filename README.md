# NCrystal DirectGA plugin

This plugin constructs the self-scattering kernel directly from a frequency
spectrum using the Gaussian approximation in the time domain. It supports
liquid spectra supplied on an absolute per-atom scale and normalized solid
VDOS data supplied through `@DYNINFO`.

## Input

Each participating element has a block in `@CUSTOM_DIRECTGA`:

### Liquid

A liquid H2 example (tables abbreviated):

```text
@DYNINFO
  element H
  fraction 1
  type sterile
@CUSTOM_DIRECTGA
element H
phase liquid
model skold
yk_model para
sq
# Q [angstrom^-1], S(Q)
0.0 0.1
...
spectrum
# energy [meV], f(E) [meV^-1 per atom]
0.0 0.05
...
```

`model` accepts `incoherent`, `skold`, or `vineyard` and defaults to
`incoherent`. The Skold and Vineyard models require an `sq` table with
strictly increasing Q values in angstrom^-1 and dimensionless S(Q) values.
For liquid H2 and D2, `yk_model` optionally selects the `ortho` or `para`
Young-Koppel rotational model. The `spectrum` line begins the f(E) table.

The complete liquid test input is in
[`data/tsl-para-H_20K_directga.ncmat`](data/tsl-para-H_20K_directga.ncmat).

### Aluminium

Append this block to an aluminium NCMAT file with its existing
`@DYNINFO type vdos` entry:

```text
@CUSTOM_DIRECTGA
element Al
phase solid
```

The complete example is in
[`data/Al_sg225_directga.ncmat`](data/Al_sg225_directga.ncmat).

For a solid, set `phase solid` and omit the spectrum rows. The element must
instead have a matching `@DYNINFO` entry of `type vdos`. Every solid VDOS
element must have its own `@CUSTOM_DIRECTGA` block. Young-Koppel, Skold, and
Vineyard corrections are currently liquid-only.

The direct_ga spectrum uses energy in meV and `f(E)` in meV^-1 per atom.

## Numerical Settings

These parameters are optional and apply separately to each element block.
The table lists their default values, used when a parameter is omitted:

| Setting | Default | Purpose |
|---|---|---|
| `max_energy_transfer_meV` | `2000` | Maximum absolute energy transfer represented by the kernel, in meV. |
| `numerical_time_end_ps` | `2` for liquids; `5` for solids | End of the numerical Fourier integral, in ps. |
| `time_step_ps` | `auto` | Automatically selected time step; a positive value specifies it explicitly in ps. |

For liquids, the analytical diffusion tail continues beyond
`numerical_time_end_ps`; later times are not discarded. For solids, the
numerical integral stops at this time. A longer interval can help with slowly
relaxing liquids or narrow solid-state spectral features, at increased
calculation time and memory use.

In auto mode, the timestep is selected from the larger of the requested
maximum energy transfer and the highest spectrum energy. Increasing
`max_energy_transfer_meV` therefore usually reduces the automatic timestep
and increases calculation time. The transfer range is not the same as the
incident neutron energy or the input spectrum cutoff.

To override the defaults, include any of these lines after `phase` and before
`sq` or `spectrum`. For example, a liquid block can include:

```text
max_energy_transfer_meV 5000
numerical_time_end_ps 4
time_step_ps 0.0001
```

These controls can be used to check convergence and numerical effects:
increase the transfer range, extend the time interval, or reduce the timestep.
Change one setting at a time to distinguish its effect. Normally leave
`time_step_ps` as `auto` unless checking timestep sensitivity.

The defaults are starting values, not an automatic accuracy guarantee for
every material. Explicit numbers must be finite and positive. An incompatible
timestep, or a time interval shorter than the selected timestep, raises an
error.

## Test

With the plugin installed:

```bash
ncrystal-pluginmanager --test DirectGA
```
