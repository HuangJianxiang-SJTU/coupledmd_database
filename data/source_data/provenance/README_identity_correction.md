# Gq_8HCQ identity-correction provenance

The companion database-paper inventory
`CoupledMD_Supplementary_Data_v13/Supplementary_Data_S1_system_inventory.csv`
is a frozen release and is not edited in place. Its `Gq_8HCQ` row labels the
receptor EDNRB/P24530.

`Supplementary_Data_S1_system_inventory_v14_identity_corrected.csv` is a
byte-for-byte copy of that table except for the audited correction:

| system_id | old identity | corrected identity |
|---|---|---|
| Gq_8HCQ | Endothelin receptor type B / P24530 | Endothelin-1 receptor / P25101 |

The correction is supported by
`Gq_8HCQ_RECEPTOR_IDENTITY_SEQUENCE_AUDIT_20260728.md`. Atlas rigor-v3 code
applies the same correction through `paper1_metadata_corrections_v3.csv`
when reading the frozen v13 table, so old releases remain traceable while
current analyses use EDNRA/P25101.
