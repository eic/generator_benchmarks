# Afterburner benchmark

Compares the beam effects of the afterburner (`abconv`) with those of Pythia8's
`eicBeamShape` (eicSimuBeamEffects), the reference it was validated against in
2021, for the beam energies and beam regimes that both define.

For every setup one Pythia8 steering file and one seed give two samples:

| sample | how | role |
| --- | --- | --- |
| `py8_full` | Pythia8 with `eicBeamShape`, crossing angle −25 mrad | reference |
| `py8_headon` | the same, momentum and vertex spread off | abconv's input |
| `afterburner` | `abconv -p <preset> py8_headon` | under test |

Because both beam-effects samples come from the same generator events, any
difference is a difference of the beam-effects models. The checks:

* **Beam geometry**, event by event (`geometry.json`, `geometry.png`): mean beam
  directions, divergences, size of the luminous region in x, y, z and t, the
  crossing term of the interaction point (IP x against the lepton's position in
  its bunch, z + ct) and the crab kicks (beam angle against the bunch position).
  Both codes set ct = (z_e − z_h)/2 and z = (z_e + z_h)/2 from the two bunch
  positions, so the slopes show which crossing geometry each part was built in.
* **Validation histograms** (`distributions.json`, `distributions.png`): abconv's
  own histograms of the three samples, overlaid with the ratio, as in
  `afterburner/python/comparison.ipynb`. Information only.
* **Displaced vertices** (`vertices.json`): decays more than 1 mm from their
  production vertex must point along the parent's momentum after the transform.
* **Seed** (`seed.json`): abconv runs on head-on files with different Pythia8
  seeds must not repeat the same beam configurations.
* **Collision time from the vertex** (in `geometry.json`): how well a linear fit
  on the vertex coordinates gives ct. This is information only; the x term carries
  the sign of the crossing geometry.

A check passes when the afterburner value agrees with the Pythia8 value within
5 combined statistical errors, or within a relative tolerance: 2 % for
horizontal beam directions, 5 % for vertical ones, 3 % for widths, 5 % for
slopes. This floor keeps large samples from flagging tiny differences.
`scripts/verdict.py` holds the rule, and `scripts/geometry.py` (`CHECKS`) the
tolerances.

## Layout

| Path | Description |
| --- | --- |
| `setups.yml` | setups: steering file, eicBeamShape configuration, abconv preset, CI flag |
| `Snakefile` | generation, abconv, analysis and summary rules |
| `config.yml` | GitLab CI jobs |
| `benchmark.json` | benchmark definition (name, target pass fraction) |
| `scripts/hepmc_text.py` | HepMC3 text readers (numpy only) |
| `scripts/geometry.py` | beam geometry, crossing slopes, crab kicks, collision time |
| `scripts/distributions.py` | abconv validation histograms, overlays and ratios |
| `scripts/vertices.py` | displaced decay vertices |
| `scripts/seed.py` | seed dependence of abconv and of Pythia8 |
| `scripts/report.py` | per-setup `report.md` and `tests.json` (common_bench test format) |
| `scripts/summary.py` | `summary.md` over all setups and `results/afterburner.json` |

## Setups

| setup | Pythia8 steering file | eicBeamShape | abconv preset | CI |
| --- | --- | --- | --- | --- |
| 18x275_hidiv | dis_eicBeam_hiDiv_18x275_10to100 | 1 | ip6_hidiv_275x18 | yes |
| 10x275_hidiv | dis_eicBeam_hiDiv_10x275_10to100 | 1 | ip6_hidiv_275x10 | |
| 10x100_hidiv | dis_eicBeam_hiDiv_10x100_10to100 | 1 | ip6_hidiv_100x10 | yes |
| 5x100_hidiv | dis_eicBeam_hiDiv_5x100_10to100 | 1 | ip6_hidiv_100x5 | |
| 5x41_hidiv | dis_eicBeam_hiDiv_5x41_10to100 | 1 | ip6_hidiv_41x5 | yes |
| 18x275_hiacc | dis_eicBeam_hiAcc_18x275 | 2 | ip6_hiacc_275x18 | |
| 10x275_hiacc | dis_eicBeam_hiAcc_10x275 | 2 | ip6_hiacc_275x10 | yes |
| 10x100_hiacc | dis_eicBeam_hiAcc_10x100 | 2 | ip6_hiacc_100x10 | |
| 5x100_hiacc | dis_eicBeam_hiAcc_5x100 | 2 | ip6_hiacc_100x5 | |
| 5x41_hiacc | dis_eicBeam_hiAcc_5x41 | 2 | ip6_hiacc_41x5 | |
| 18x110_eA | dis_eicBeam_eA_18x110 | 3 | ip6_eAu_110x18 | yes |
| 10x110_eA | dis_eicBeam_eA_10x110 | 3 | ip6_eAu_110x10 | |
| 5x110_eA | dis_eicBeam_eA_5x110 | 3 | ip6_eAu_110x5 | |
| 5x41_eA | dis_eicBeam_eA_5x41 | 3 | ip6_eAu_41x5 | |
| 10x130_hiacc | dis_eicBeam_hiAcc_10x100 at 10x130 | 2 | ip6_ep_130x10 | yes |
| 9x130_hiacc | dis_eicBeam_hiAcc_10x100 at 9x130 | 2 | ip6_ep_130x9 | |
| 10x250_hiacc | dis_eicBeam_hiAcc_10x275 at 10x250 | 2 | ip6_ep_250x10 | |
| 9x275_hidiv | dis_eicBeam_hiDiv_9x275_10to100 | 1 | ip6_hidiv_275x9 | |
| 9x275_hiacc | dis_eicBeam_hiAcc_10x275 at 9x275 | 2 | ip6_hiacc_275x9 | |

Notes on the setups:

* The beam energies always come from the setup name; the benchmark appends them
  to the steering file. eicBeamShape parametrizes the early-science energies as
  copies of 10x100 and 10x275 (its comments say "Copy ... for early science"),
  and the setups use the steering file of the energy it copies.
* abconv's `ip6_ep_130x10`, `ip6_ep_130x9` and `ip6_ep_250x10` carry the
  high-acceptance parameters. eicBeamShape's high-divergence 10x130 and 9x130,
  for which steering files exist, have no abconv counterpart, so they are not
  setups. abconv's `ip6_hidiv_100x9` and `ip6_hiacc_100x9` have no eicBeamShape
  counterpart.
* For eA, Pythia8 runs e+p at the per-nucleon energy with the e+Au beam
  parameters. Only the beam effects are compared.

## Running

In `eic-shell` (or the `eic_xl` container), from the repository root:

```console
snakemake --cores 4 afterburner
```

This runs the CI selection (`ci: true` in `setups.yml`) with 2000 events per
sample and writes `results/afterburner/summary.md`. The full matrix with
analysis-size samples:

```console
snakemake --cores 16 afterburner --config AFTERBURNER_SETUPS=all AFTERBURNER_N_EVENTS=50000
```

One setup, or a list of them:

```console
snakemake --cores 4 results/afterburner/10x100_hidiv/report.md
snakemake --cores 8 afterburner --config AFTERBURNER_SETUPS=10x100_hidiv,10x100_hiacc
```

Outside the container, the rules carry the repository's `singularity:` image
(`EIC_SINGULARITY_CONTAINER`, default `eic_xl:nightly` from cvmfs), so
`snakemake --cores 4 afterburner --use-apptainer` should work where cvmfs and
apptainer are available. This has not been tested; the tests ran inside
eic_xl.

The steering files are downloaded from GitHub. Without network access, put
them in `sim_output/afterburner/steerFiles/` by hand (from eicSimuBeamEffects
`Pythia8/steerFiles/` at `AFTERBURNER_STEERING_REF`); the download rule is
then skipped.

Options (`--config KEY=VALUE`):

| key | default | meaning |
| --- | --- | --- |
| `AFTERBURNER_SETUPS` | `ci` | `ci`, `all`, or a comma-separated list of setup names |
| `AFTERBURNER_N_EVENTS` | 2000 | events per Pythia8 sample |
| `AFTERBURNER_N_SEED_CHECK` | min(1000, N) | events of the second-seed samples |
| `AFTERBURNER_SEED` | `BENCHMARK_RNG_SEED` (1) | Pythia8 seed; the seed check uses seed + 1 |
| `AFTERBURNER_STEERING_REF` | `v1.1` | eicSimuBeamEffects ref of the steering files; keep it at the version the container ships, since eicBeamShape stops on steering values it does not expect |

Changing an option re-runs the affected rules (Snakemake tracks the parameters).

For orientation (eic_xl 26.11.0-stable on a 24-core workstation):

* the CI selection (6 setups, 2000 events) takes 1–2 minutes on 6 cores;
* all 19 setups at 20000 events take 5 minutes on 10 cores.

The HepMC3 text files take about 25 MB per 1000 events per setup in
`sim_output/` (9 GB for the full matrix at 20000 events).

## Outputs

| file | content |
| --- | --- |
| `results/afterburner/summary.md` | one row per setup: checks passed and verdict per group, IP crossing slope |
| `results/afterburner.json` | benchmark result: every test, pass fraction, pass/fail against `target` |
| `results/afterburner/<setup>/report.md` | the setup's tables and figures |
| `results/afterburner/<setup>/tests.json` | the setup's checks, common_bench test format |
| `results/afterburner/<setup>/geometry.png` | profiles: IP x vs z + ct, hadron and electron crab kicks |
| `results/afterburner/<setup>/distributions.png` | abconv validation histograms and ratios |
| `results/afterburner/versions.txt` | abconv, Pythia8 and container versions |
| `sim_output/afterburner/<setup>/` | steering files, HepMC3 samples, abconv histogram files |
| `logs/afterburner/<setup>/` | Pythia8 and abconv logs |

The analysis scripts run on any HepMC3 text files; see each script's `--help`.

## Status

Written against eic_xl 26.11.0-stable (abconv 0.2.1, Pythia8 8.316,
eicSimuBeamEffects 1.1). With all 19 setups at 20000 events, 260 of 361
checks pass.

What agrees in every setup:

* mean beam directions, horizontal and vertical (hadron at −25.0 mrad);
* all four divergences;
* horizontal, longitudinal and time widths of the luminous region;
* both crab kicks, sign and size.

What fails:

* **Interaction point crossing slope: opposite sign, every setup.** abconv
  gives IP x = +(θ/2)(z + ct), Pythia8 −(θ/2)(z + ct): +12.6…+12.9 against
  −12.2…−12.7 mrad. abconv builds the IP in the configured +25 mrad frame,
  while its beams leave at −25 mrad like Pythia8's. The IP code is
  eicBeamShape's, run with the opposite angle.
* **Beam energy spread: absent, every setup.** abconv applies none; the
  relative rms is 1–6 × 10⁻⁶, from the angles alone. Pythia8 applies
  (6–11) × 10⁻⁴ as set in the steering files. This stays a failing check, not
  information: the spread is part of the beam effects abconv is meant to
  reproduce (38 of the 101 failures).
* **Displaced vertices: misaligned, every setup.** Decay vertices stay in the
  head-on frame, 25 mrad off their momenta. afterburner PR #38 fixes this.
* **Seed: repeated, every setup.** The smearer has a fixed seed, so abconv runs
  on different inputs repeat the same interaction points and beam four-vectors.
* **Vertical luminous width: 4–8 % apart in 6 setups** (10x100 both regimes,
  18x275 hidiv, the eA setups except 18x110). abconv computes it from β*·ε of
  the two beams; eicBeamShape takes the steering file's `sigmaVertexY` (the CDR
  tables). For example, 10x100_hidiv gives 8.70 µm in abconv and 8.04 µm in
  Pythia8.

Information worth knowing:

* abconv shifts the beam energies, by the same amount in every setup: hadron
  +2.34 × 10⁻⁴, electron −0.78 × 10⁻⁴. Together they make θ²/4. Pythia8 does
  not shift them.
* A vertex-based collision time is 23–30 ps from z alone (the electron bunch
  length) and 14–26 ps with x added. A (z, x) estimator fitted on Pythia8 gives
  30–50 ps on abconv samples, worse than z alone, because of the IP sign.

The collect job is therefore `allow_failure: true`: `afterburner_check` exits
non-zero while the benchmark result is fail, without blocking the pipeline.
