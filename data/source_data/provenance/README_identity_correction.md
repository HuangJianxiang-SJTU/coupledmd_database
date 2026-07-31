# Gq_8HCQ identity-correction provenance

The companion database-paper inventory
`Supplementary_Data_S1_system_inventory.csv` is a frozen release and is
not edited in place. Its `Gq_8HCQ` row labels the receptor EDNRB/P24530.

`Supplementary_Data_S1_system_inventory_identity_corrected.csv` is a
byte-for-byte copy of that table except for the audited correction:

| system_id | old identity | corrected identity |
|---|---|---|
| Gq_8HCQ | Endothelin receptor type B / P24530 | Endothelin-1 receptor / P25101 |

The correction is supported by
`Gq_8HCQ_RECEPTOR_IDENTITY_SEQUENCE_AUDIT_20260728.md`. The atlas
analysis applies the same correction through the metadata-corrections
table when reading the frozen inventory, so old releases remain
traceable while current analyses use EDNRA/P25101.
