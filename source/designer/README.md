# DNA Moiré Designer prototype

The second workflow page now implements a staged, auditable structure build:

1. when Scaffold routing is generated, silently create and validate the fixed
   two-layer capture SST on virtual helices 48–63 at scaffold positions
   48–175 and 208–335; there is no separate SST review stage;
2. generate two balanced left/right Path-view scaffold bands (7548 nt each).
   Native Square 10/11-bp endpoint phases redistribute capacity between the
   two sides: the left endpoints span only bases 35--43 and the right endpoints
   only 348--356. All 48 Seed helices remain covered, including the complete
   capture regions, before cadnano review or expert JSON replacement;
3. only after explicit acceptance, generate staple crossovers, run Autobreak,
   and install the validated two-layer capture topology and pair colors from
   `S-seed-pore_3L.json`;
4. after sequences are saved in cadnano, export isolated capture and intact
   output XLSX/JSON/SVG bundles. The intact output topology is a temporary
   same-coordinate snapshot, so no extra SST copy appears in the live design.

Returning to an earlier stage preserves downstream results until the edited
file is accepted again.  If its topology fingerprint is unchanged, the prior
downstream design resumes automatically.  If routing, crossovers, indels,
coordinates, or effective ranges changed, the UI warns before invalidating the
dependent stages.  Files remain on disk as versioned audit records.

Capacity balancing may extend or trim either edge, including non-capture
helices, but only at native Square 10/11-bp phases. Two edge modules are
recognized by two complete legal 32-bp crossover domains, not by the raw
inclusive endpoint difference (the apparent 63/64 ambiguity). A single
extended helix pair stays one edge domain. Only a run of at least two adjacent
helix pairs that each contains two complete 32-bp domains is split, and its
second module is joined to the nearest main scaffold by a legal seam. Orphan
32-bp components fail validation.

This left/right capacity balancing and edge-module seam policy belongs only
to the Moiré Designer structure worker. It never changes, wraps, or supplies
defaults to cadnano's normal `AutoCS_scaffolds` command.

The independent data-analysis module provides TEM/FFT image analysis and does
not require parameter, structure, or sequence acceptance. TEM mode detects the scale
bar with OCR, estimates the mean real-space lattice constant and moire period,
and calculates the equal-square-lattice twist as
`2*asin(lattice_constant/(2*moire_period))`.  FFT mode separates the dominant
reciprocal peaks into two square-lattice orientation families.  Combined mode
uses the TEM real-space period for the reported twist and retains the FFT fit
as an independent cross-check.  Annotated PNGs and a JSON analysis record are
saved next to the uploaded image; every automatic value remains editable.

The first prototype implements one calibrated workflow:

- Seed S(F), 48-helix square frame
- Square–Square bilayer
- live twist-angle / moiré-period / lattice-constant coupling
- phase-coupled Z1–Z2–Z3 selection (8-bp domains, 32-bp repeat)
- S8-R4x4C calibrated twist/indel conversion
- fully cooperative capture validation
- interactive lightweight 3D preview
- `.moire.json` project persistence
- reference Seed-S JSON and capture-map export
- manual simulation/experiment measurement comparison

Run from a Python environment containing the dependencies listed by the
platform build instructions:

```bash
PYTHONPATH=. python run_moire_designer.py
```

The 0.2 solver uses the direct S8-R4x4C calibration. Native twist is
`0.1030236095 degrees/base`, so 32 bp corresponds to `3.2967555 degrees`.
Z lengths remain independent inputs: changing Z2 preserves the target twist and
updates the required mean indel. Sequence assignment and chain-level
regeneration remain pending; changing Z metadata does not yet rebuild strand
indices in the reference cadnano JSON.
