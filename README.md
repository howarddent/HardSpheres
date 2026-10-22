# HardSpheres

An ideal gas, and at higher density a non-ideal one, simulated as hard spheres in a box with a piston. The simulation runs in OpenCL on the GPU, and OpenGL draws the spheres in real time from the same buffer.

![100,000 hard spheres at packing fraction 0.05, coloured by speed, with the speed distribution against Maxwell-Boltzmann](docs/hardspheres.png)

*100,000 spheres at packing fraction 0.05, coloured by speed. The orange top face is the piston. On the right are the controls, the readouts (pressure against Carnahan–Starling), and the speed histogram against Maxwell–Boltzmann.*

The sticky version in [`Sticky/`](Sticky/) adds attraction, binding and gravity:

![StickySpheres: a Lennard-Jones liquid wetting the floor under its vapour](docs/sticky-pool.png)

## Build and run

```
lazbuild --lazarusdir=$HOME/Lazarus --ws=gtk2 HardSpheres.lpi
./run-nvidia.sh            # OpenGL on the NVIDIA GPU: CL-GL buffer sharing
./HardSpheres              # default GL (Intel iGPU): positions copied each frame
fpc -O2 gascheck.lpr && ./gascheck [N]     # headless physics checks
```

`gas.cl` must sit beside the executable; it is compiled at start-up.

On a hybrid-graphics laptop, OpenGL runs on the integrated GPU unless NVIDIA PRIME offload is requested. OpenCL–OpenGL sharing needs both on the same device. `run-nvidia.sh` sets `__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia`. Without it the program still runs, copying positions through the host, and says so in the status line.

## Physics

Units: sphere diameter 1, mass 1, k_B = 1.

- **Time stepping:** each substep is short enough that a typical sphere moves 0.02 diameters. A frame is 0.05 time units, so a hotter gas takes more substeps per frame and visibly moves faster.
- **Collisions** (`gas.cl`):
  - Spheres are binned into a uniform grid.
  - Each sphere picks the overlapping, approaching neighbour it touched earliest.
  - Mutual picks collide exactly as they would have at contact: the impulse acts along the line of centres at the moment of touching, and each sphere gets the displacement since then.
  - One collision per sphere per substep keeps energy exact. Summing simultaneous impulses does not: an early version gained 78 % in temperature at packing fraction 0.3.
- **Walls** reflect. The piston is the top face (z = Lz), so only the box's height changes with its volume; the floor and side walls stay put. It reflects in its own moving frame, so the gas heats when it advances and cools when it withdraws.
- **Pressure** is the momentum given to the walls per unit time and area. The box dimensions used are those open to the centres, L − 1, so an ideal gas gives exactly P = NT/V.

`gascheck` results (20k–160k spheres):

| check | result |
|---|---|
| energy, no thermostat | T drift 3e-6 (packing 0.05), 3e-4 (0.3) over 20 time units |
| equal speeds → Maxwell–Boltzmann | ⟨v⁴⟩/⟨v²⟩² 1.000 → 1.664 (5/3) |
| P/(n T) vs Carnahan–Starling | −0.7 % at packing 0.05, −0.5 % at 0.3 (160k, bulk density) |
| collision rate vs Enskog | −1.3 %, −2.4 % |
| slow adiabatic compression, V/2 | T × 1.636–1.645; ideal 1.598, hard-sphere correction about ×1.027 → 1.641 |

Near the walls a dense hard-sphere fluid piles up, so the bulk density is below N/V. At packing fraction 0.3 in a 34-diameter box this puts P V/(N T) 9 % below Carnahan–Starling. The check therefore uses the density four diameters in from the walls.

A compression is quasi-static only if the piston is slow compared with sound (√(5T/3) ≈ 1.3 at T = 1), and if the move takes longer than a sound wave needs to cross the box. A faster piston heats the gas beyond the adiabat. The `TV^(2/3)` readout shows this.

## Files

- `gas.cl`: kernels (`move`, `clear`, `bin`, `pairs`, `apply`, `scale`).
- `ugassim.pas`: `TGasSim`, the simulation: buffers, frames, measurements, the piston and the thermostat.
- `clgl.pas`: runtime OpenCL bindings, including `cl_khr_gl_sharing`.
- `umain.pas`: the window: shader-drawn sphere sprites, orbit camera, sliders, readouts, and the speed histogram against Maxwell–Boltzmann.
- `gascheck.lpr`: headless checks.

Command line: `--n=`, `--phi=`, `--maxwell`, `--hold`, `--speed=`, `--compress=<f>`, `--snapshot=<png>` with `--frames=`.

## Next

- Non-ideal gas with attraction: done as a square well in `Sticky/` (StickySpheres), with binding, inelastic capture, breaking, clusters and condensation. See its README.
- A counting-sort cell list in place of fixed-capacity cells, so memory follows N rather than the box volume.

## Licence

MIT; see [LICENSE](LICENSE). Covers `Sticky/` too.
