# v12 database-versus-resource claim firewall

Generated: 2026-07-20T06:56:34.928504+00:00

## Decision

The v12 Data Descriptor is restricted to data-resource claims. The previous
GPCR-centred pocket-reuse atlas has been removed. The previous family-organized
gateway distribution has been removed from the main figures. Figure 3 retains
orthosteric recovery only as a technical positive control; Figure 4 contains
only repeat-level and reduced-record technical validation.

## Ownership matrix

| Topic | Database Data Descriptor | Resource-atlas paper | AI paper |
|---|---|---|---|
| Cohort construction, identifiers and exclusions | owns | cites | cites |
| Simulation protocol and provenance | owns | summarizes as input | summarizes as input |
| Three computational repeats and QC | owns | uses as input | uses as input |
| Reduced PDB/XTC records, checksums and access | owns | links | links |
| Pocket method and 49/58 positive control | owns limited validation | may cite | may cite |
| 2,149 pocket rows; 205 detected + 2 zero | owns completeness | uses as input | uses as input |
| Recurrent pocket atlas and 65 clusters | excluded | owns | excluded |
| Pocket rankings and biological recurrence | excluded | owns | excluded |
| C7–Gq enrichment | excluded | owns | excluded |
| Partner-associated pocket changes | excluded | owns | excluded |
| Gateway method and processed S6 records | owns methods/intermediate | uses with disclosure | excluded |
| Gateway family distribution and interpretation | excluded | owns | excluded |
| α5 and family-selectivity analyses | excluded | owns | model input only if disclosed |
| Mechanism or druggability conclusions | excluded | owns only if supported | excluded |
| Learned models and architecture | excluded | excluded | owns |
| Precog3D benchmarking and predictions | excluded | excluded | owns |
| Prospective panels, BRET, rescue or epistasis | excluded | excluded | owns |

## Shared-data disclosure

The resource paper may reuse the final-207 cohort table, source-simulation
provenance, processed pocket records and processed gateway intermediates.
Reuse must cite the Data Descriptor and disclose that the database paper owns
the cohort, data-record and technical-validation descriptions. The database
paper does not reproduce the resource paper's biological figures, rankings,
cluster headline, family enrichment or mechanistic interpretation.

## Automated prohibited-content gate

The v12 manuscript is checked for model architectures, Precog3D, BRET,
prospective challenge panels, model-selected mutations, rescue, epistasis,
the 65-cluster headline, C7–Gq enrichment, pocket rankings and gateway-family
biological claims. Incidental use of “model” for starting structural models is
not an AI claim and is reviewed in context.
