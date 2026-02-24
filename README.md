# Codes for Entropy-Driven Physical Amplification in Multivalent Biosensing

This repository contains C++ codes for the theoretical and simulation studies of entropy-driven physical amplification in multivalent biosensing systems. The codes compute surface coverage (θ) and the amplification factor (α) for multivalent guest particles binding to receptor-decorated substrates, using both analytical self-consistent field theory and Monte Carlo (MC) simulations.

## Repository Structure

```
new_code/
├── code/                      # Main binding model (no misbinding)
│   ├── theory/
│   │   ├── system.h           # System class with self-consistent equations and derivative calculations
│   │   └── theta_mu.cpp       # Theoretical θ and α vs. μ_linker (Poisson-averaged over nr)
│   └── simulation/
│       ├── guest.h            # Guest particle class (position, linker state, MC moves)
│       ├── substrate.h        # Substrate/receptor class
│       ├── array_computation.h# Utility functions (norms, periodic boundary, etc.)
│       ├── system.h           # MC simulation system class (cell lists, MC moves)
│       ├── theta_rho.cpp      # MC simulation: θ vs. μ_linker
│       ├── alpha_rho.cpp      # MC simulation: α vs. μ_linker (finite difference over ρ_r)
│       ├── alpha_xi.cpp       # MC simulation: α vs. μ_linker (finite difference over ξ)
│       ├── reverse_ad.cpp     # MC simulation: find adsorption threshold μ_linker
│       └── reverse_dp.cpp     # MC simulation: find desorption threshold μ_linker
└── misbinding_code/           # Extended model with non-specific (misbinding) interactions
    ├── theory/
    │   ├── binder_system.h    # System class extended for misbinding
    │   └── binder_mu.cpp      # Theoretical θ and α vs. μ_linker with misbinding
    └── simulation/
        ├── guest.h            # Guest particle class (with misbinding)
        ├── substrate.h        # Substrate class (with misbinding)
        ├── array_computation.h# Utility functions
        ├── system.h           # MC simulation system class (with misbinding)
        └── theta_rho.cpp      # MC simulation: θ vs. μ_linker with misbinding
```

## Physical Model

The system models multivalent guest particles (e.g., nanoparticles carrying multiple ligands) binding to a substrate decorated with receptors. Each guest carries `nl` ligands, each capable of forming up to `kl` bonds. Each receptor on the substrate can form up to `kr` bonds with ligands.

Key quantities computed:

- **θ** (theta): Surface coverage fraction (fraction of substrate patches occupied by a guest)
- **α** (alpha): Amplification factor — the logarithmic sensitivity of θ to changes in `μ_linker`
- **F_eff**: Effective free energy of binding for a guest particle near the substrate
- **q**: Mayer-f function defined by `q = exp(-F_eff) - 1`, related to θ by `θ = z_g·q / (1 + z_g·q)`

## Parameters

| Parameter    | Description |
|-------------|-------------|
| `xi_l`      | Boltzmann weight for ligand bond formation, `exp(-f_l)` where `f_l` is the ligand binding free energy |
| `kl`        | Maximum number of bonds a single ligand can form |
| `nl`        | Number of ligands per guest particle |
| `xi_r`      | Boltzmann weight for receptor bond formation, `exp(-f_r)` |
| `kr`        | Maximum number of bonds a single receptor can form |
| `nr`        | Number of receptors within reach of a guest particle (Poisson-distributed with mean `nr` in theory) |
| `mu_linker` | Chemical potential of the linker (log concentration; varied to sweep θ from 0 to 1) |
| `mu_g`      | Chemical potential of the guest particle in bulk |
| `rho_r`     | Surface density of receptors |
| `fcnf`      | Confinement free energy cost for bridging (guest-substrate bond formation) |
| `area`      | Simulation box area |
| `xi_lm`, `xi_rm`, `mu_m` | *(misbinding model only)* Misbinding Boltzmann weights and chemical potential |

## Code Descriptions

### Theory (`new_code/code/theory/`)

**`theta_mu.cpp`** — Computes θ and α as a function of `μ_linker` using analytical self-consistent equations. For each value of `μ_linker`, results are averaged over a Poisson distribution of receptor numbers `nr`. Outputs `theta.txt` and `alpha.txt`.

Command-line arguments:
```
./theta_mu <-ln(xi_l)> <kl> <nl> <-ln(xi_r)> <kr> <nr>
```

### Simulations (`new_code/code/simulation/`)

All simulation codes use a grand-canonical Monte Carlo scheme: guest particles are inserted/removed from the simulation box based on the Metropolis criterion, with the effective free energy `F_eff` computed from the self-consistent equations at each step.

**`theta_rho.cpp`** — Computes θ as a function of `μ_linker`. Outputs `theta.txt`.

**`alpha_rho.cpp`** — Computes α as a function of `μ_linker` using a four-point finite difference with respect to `μ_linker`. Outputs `alpha.txt`.

**`alpha_xi.cpp`** — Computes α as a function of `μ_linker` using a two-point finite difference with respect to `ln(ξ)`. Outputs `alpha.txt`.

**`reverse_ad.cpp`** — Finds the adsorption threshold value of `μ_linker` (where θ crosses 0.1 from below) using bisection. Outputs `adsorption.txt`.

**`reverse_dp.cpp`** — Finds the desorption threshold value of `μ_linker` (where θ crosses 0.1 from above) using bisection. Outputs `desorption.txt`.

Command-line arguments for all simulation codes:
```
./program <mu_g> <rho_r> <kl> <kr> <nl> <area> <xi>
```

### Misbinding Model (`new_code/misbinding_code/`)

These codes extend the main model to include non-specific (misbinding) interactions, where ligands and receptors can form weaker off-target bonds characterized by `xi_lm`, `xi_rm`, and `mu_m`.

- **`theory/binder_mu.cpp`**: Analytical θ and α vs. `μ_linker` with misbinding. Command-line arguments: `<-ln(xi_l)> <kl> <nl> <-ln(xi_r)> <kr> <nr> <mu_m> <-ln(xi_m)>`
- **`simulation/theta_rho.cpp`**: MC simulation of θ vs. `μ_linker` with misbinding. Command-line arguments: `<mu_g> <rho_r> <kl> <kr> <nl> <area> <xi> <mu_m> <xi_m>`

## Building

The codes are written in standard C++11. Compile each program individually, for example:

```bash
g++ -O2 -o theta_mu new_code/code/theory/theta_mu.cpp -lm
g++ -O2 -o theta_rho new_code/code/simulation/theta_rho.cpp -lm
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
