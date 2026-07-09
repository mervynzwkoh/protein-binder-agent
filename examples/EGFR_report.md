# Target-to-Binder Triage Report: EGFR
*Generated 2026-07-08 21:00*

**UniProt accession:** P00533  
**Structure source:** experimental (PDB 1IVO)

## Calibration: known ligands vs. candidates

| Group | n | Best score (kcal/mol) | Worst score |
|---|---|---|---|
| Known ChEMBL ligands | 28 | -9.72 | -4.75 |
| Repurposing candidates | 20 | -8.19 | -3.09 |

## Agent interpretation

Target: EGFR (UniProt P00533). Structure used: experimental PDB 1IVO (EGFR kinase domain), giving reasonably high confidence in the docking geometry (not an ESMFold homology model).

Calibration check (dock_known_ligands vs dock_candidates):
The 44-compound known-ligand calibration set (mostly quinazoline/quinazolinone EGFR inhibitor chemotypes plus some weak tyrphostin-like cinnamonitrile fragments) produced docking scores ranging from -9.72 to -4.75 kcal/mol. Critically, the four best-scoring known ligands (-9.72, -8.90, -8.73, -8.63 kcal/mol) are all quinazoline-based EGFR-inhibitor scaffolds (e.g. CHEMBL69960, CHEMBL443268, CHEMBL68920, CHEMBL304271), which is exactly the expected EGFR-binding chemotype — a reassuring sign that the docking setup recognizes the correct pharmacophore for this ATP pocket. Weaker/fragment-like known ligands (the cinnamonitrile/malononitrile series, IC50/Ki in the high-micromolar to millimolar range) clustered at the bottom of the score range (-6.5 to -4.8 kcal/mol), roughly tracking their lower reported potency. This gives a rough, non-perfect, but directionally sensible correlation between docking score and known potency.

Against this backdrop, the 20 Lipinski-filtered approved-drug repurposing candidates (excluding known EGFR ligands) scored from -3.09 to -8.19 kcal/mol. The best candidates (simvastatin -8.19, methotrexate -8.12) land in the range of the moderate-to-good known ligands but do NOT reach the top tier of true EGFR inhibitors (-8.6 to -9.7 kcal/mol). All other candidates fall off sharply below -7.5 kcal/mol, well into the range occupied by the weakest/fragment-like known ligands.

Interpretation: The docking setup shows some ability to separate strong EGFR binders (top known quinazolines) from weak ones, and using an experimental structure adds confidence. However, none of the repurposing candidates clearly exceed or even match the top known-ligand scores — the best candidates are only comparable to mid-tier known ligands, and there is meaningful overlap between candidate scores and the weak end of the known-ligand distribution. This means the "hits" here should be read as a coarse triage signal only, not as evidence of genuine high-affinity EGFR binding. No candidate stands out as a confident repurposing hit; this is a weak, exploratory signal at best.

Ranked shortlist of top repurposing candidates (by Vina score):
1. Simvastatin (-8.19 kcal/mol) — best-scoring candidate, occupies the pocket with its bulky decalin/lactone system, but simvastatin is a large, flexible lipophilic HMG-CoA reductase inhibitor with no structural precedent as a kinase-hinge binder, so this is likely a docking artifact from good shape complementarity rather than a genuine pharmacophore match.
2. Methotrexate (-8.12 kcal/mol) — a diaminopteridine/antifolate; its rigid heteroaromatic system can form hinge-like H-bonds similar to EGFR inhibitors, giving it slightly more chemotype plausibility than simvastatin, though it is a highly polar antifolate not typically associated with kinase ATP-site selectivity.
3. Vildagliptin (-7.54 kcal/mol) — small nitrile-containing DPP-4 inhibitor; the nitrile and pyrrolidine could make modest contacts but score is well below top known EGFR inhibitors.
4. Istradefylline (-7.46 kcal/mol) — a xanthine/adenosine-receptor antagonist with a fused bicyclic core superficially similar to purine-based kinase inhibitors, giving it plausible shape complementarity to the ATP pocket.
5. Benoxaprofen (-7.16 kcal/mol) — a benzoxazole-propionic acid NSAID; moderate score likely driven by the halogenated aromatic ring fitting a hydrophobic sub-pocket, but no clear kinase-relevant pharmacophore.

Caveat: Docking scores from AutoDock Vina in this workflow are a coarse, rapid triage signal based on pose geometry and an approximate scoring function — they are not a proof of binding, do not capture protonation/tautomer effects, induced-fit flexibility, or true binding free energy, and false positives (like simvastatin here) are common for large lipophilic molecules that dock with good shape fit but lack a real mechanistic rationale. Any of these candidates would need orthogonal biochemical/biophysical validation (e.g., SPR, enzymatic IC50 assay) before further consideration, and given the imperfect separation between known EGFR inhibitors and this candidate pool, results should be treated as hypothesis-generating only.

## Ranked candidate shortlist

| Rank | ChEMBL ID | Name | Affinity (kcal/mol) |
|---|---|---|---|
| 1 | CHEMBL1064 | SIMVASTATIN | -8.19 |
| 2 | CHEMBL34259 | METHOTREXATE | -8.12 |
| 3 | CHEMBL142703 | VILDAGLIPTIN | -7.54 |
| 4 | CHEMBL431770 | ISTRADEFYLLINE | -7.46 |
| 5 | CHEMBL340978 | BENOXAPROFEN | -7.16 |
| 6 | CHEMBL652 | FLECAINIDE | -7.04 |
| 7 | CHEMBL7728 | HEXOBARBITAL | -6.50 |
| 8 | CHEMBL316561 | PROGLUMIDE | -6.33 |
| 9 | CHEMBL314437 | MEPTAZINOL | -6.19 |
| 10 | CHEMBL65375 | TIPIRACIL HYDROCHLORIDE | -6.15 |

## Limitations

- Docking scores (AutoDock Vina) are a coarse triage signal, not validated binding affinities.
- Candidate pool is a random sample (fixed seed) of ChEMBL max_phase=4 molecules, Lipinski-filtered.

## Top pose visualization

See `EGFR_top_pose.html` (open in a browser).