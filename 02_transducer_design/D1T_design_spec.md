# Dual-polarization EMAT — design specification (D1T stack, as used in the manuscript)

All dimensions in mm; z along the pipe axis, z = 0 at the probe centre; r radial.
Specimen: 7075-T6 tube, outer radius 6.5, wall thickness 1.0 (Table 1).

## Magnet stack (Section 3.1, Fig. 4(a), (b))

| Item | Value |
|---|---|
| Magnets | 3 NdFeB rings, N45H (B_r ≈ 1.35 T), **radially magnetized** (not diametrically), NiCuNi plated |
| Ring size | inner Ø15, outer Ø25, axial width 8 (tolerance: ID +0.1/0, OD 0/−0.1, width ±0.1) |
| Polarity | alternating: outward / inward / outward (+ − +) |
| Pitch | 10.4 = λ_w (centre to centre); non-magnetic spacer 2.4 thick between rings |
| Ring centres | z = −10.4, 0, +10.4 |
| Attraction between adjacent rings | ≈ 65 N (stack is self-locking; assemble with a non-magnetic push rod) |

## Coils (Section 3.1, Fig. 4(c), (d), Supplementary S3.1)

Both coil sets are wound on one thin-walled bobbin (inner Ø13.2, outer Ø14.0, wall 0.4, length 52; PTFE,
polyimide, FR4 or machined POM) with 0.15 mm polyurethane enamelled copper wire; wire centre radius ≈ 7.0.

| Set | Position | Elements | Element width | Turns per element | Winding |
|---|---|---|---|---|---|
| Wall-mode set | under each magnet (radial field → axial force) | 3, at z = −10.4, 0, +10.4 | 5.25 (= λ_w/2) | 35, single layer | adjacent elements opposite (+ − +), in series |
| Liquid-mode set | under each gap (axial field → radial force) | 2, at z = ±5.2 | 2.85 | 19 per layer, 2 layers (38) | the two elements opposite, in series |

- Gap between neighbouring wall and liquid elements: 1.15.
- Force spectra in `data/E4_paperfig_data.csv` and Fig. 4(c), (d) are per ampere for one layer of the liquid
  set (19 turns); the second layer doubles the liquid-set force.
- Design values at the operating wavelengths (one layer): F_z,wall(10.4) = 1.051 N/A; F_r,liq(5.2) = 0.444 N/A;
  F_r,liq(6.3 / 4.5) = 0.284 / 0.173 N/A.
- Expected impedance on the pipe: wall set ≈ 31–43 µH at 520 kHz; liquid set (2 × 38 turns) ≈ 30–42 µH at
  280 kHz; measure with an LCR meter on the pipe before matching.

## Electronics and acquisition (Section 4.1, Supplementary S3.3, S4.1)

| Item | Setting |
|---|---|
| Excitation | AFG31022 (Tektronix), 5-cycle sine bursts, 520 kHz (wall-mode set) and 280 kHz (liquid-mode set), fired in turn |
| Power amplifier | GA2500 (Ritec), gated; series-resonant matching TEM-128 |
| Reception | wall-mode coil of the receiving probe for both families; matching REMP-128; preamplifier PXPA6, 46 dB |
| Oscilloscope | MSO24 (Tektronix), 125 MS/s, high-resolution mode, 512 averages, 200 µs record with 10 µs pre-trigger; CH1 transmit monitor, CH2 wall-mode receive coil, CH3 liquid-mode receive coil |
| Gauge length | 100 (probe centre to centre), fixed by an acrylic connecting tube |

## Holder (Fig. 4(a), `cad/`)

Two printed half-shells clamp the bobbin and the magnet stack around the pipe and are closed by
4 × M3 screws with nuts (screw axes at z = ±12.0 from the probe centre). Ribs in the magnet pockets fix the
10.4 pitch; the shells span 49.6 axially and reach r = 24.2. The lower half-shell carries the coaxial socket.
Check the dimensions against the parts before machining.
