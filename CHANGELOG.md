### Metalink 0.3.12

Release notes for Metalink 0.3.12.

- **Updated default dRep version**: `metalink/envs/metalink-annotation.yaml` now defaults to `drep =4.0.3` (previously 3.5.0), speeding up mMAG dereplication; the corresponding conda environment name is updated from `metalink-annotation-0.3` to `metalink-annotation-0.3.1`.
- **Fixed sample- and project-level merge R scripts**: table merging could fail in some cases.
