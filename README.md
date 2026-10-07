# MLST-Typer (extracted and generalized from Pneumo-Typer)

This package isolates the MLST sequence-typing workflow from the supplied Pneumo-Typer code and generalizes its database updater using the supplied `data_mlst` scheme table.

## Key changes

1. **Select an MLST database by species name** from `data/mlst_schemes_all.tab`.
2. If a species has multiple MLST schemes, choose the scheme explicitly (or interactively in a terminal).
3. **Only the selected scheme is updated**, and each scheme has its own database directory.
4. **PubMLST updating now follows the original Pneumo-Typer OAuth workflow**: consumer key/secret -> request token (only if no access token) -> browser authorization -> persistent access token -> renewable session token -> signed downloads of scheme metadata, `alleles_fasta`, and `profiles_csv`.
5. **FASTA + GBK/GenBank + mixed directories** are supported. Each input is detected from file content and converted internally to normalized FASTA. Original files are not renamed or modified.
6. Locus number is read dynamically from the selected scheme. It is no longer hard-coded to seven loci.
7. Output preserves the original simple `ST_out.txt` idea and adds detailed tables.

## Requirements

- Perl 5
- NCBI BLAST+ (`makeblastdb`, `blastn`)
- Network access when updating a database
- The OAuth Perl dependencies are bundled from the supplied Pneumo-Typer package under `lib/`

No BioPerl is required for the MLST-only workflow.


## PubMLST authorization: one global token set

The first extracted version incorrectly replaced the authenticated update with anonymous HTTP downloads. This revision restores the original Pneumo-Typer mechanism.

The package uses one shared directory:

```text
auth/
├── access_token
├── session_token      # generated/renewed automatically
└── request_token      # temporary; only used during first authorization
```

`access_token` and `session_token` are **not stored under each species or scheme**. The same global token set is reused when switching among PubMLST databases. The supplied existing `access_token` from Pneumo-Typer is retained in `auth/access_token`; therefore, when it remains valid, the updater normally only has to obtain/renew `session_token`. If a protected request returns HTTP 401, only the global `session_token` is discarded and renewed from the persistent `access_token`.

The original PubMLST OAuth client credentials are also retained inside `RESTful_API_download_pubmlst.pl`; they are not printed in normal output.

Database files themselves remain separated by scheme, for example `database/abaumannii/` and `database/abaumannii_2/`. This separation concerns MLST data only, not authorization.

### BIGSdb-Pasteur entries

The supplied scheme table also contains a small number of schemes hosted at `bigsdb.pasteur.fr`. Those are a different BIGSdb site and require their own authenticated API credentials; the PubMLST OAuth token is not silently reused for them. The updater stops with an explicit message for such entries rather than downloading an incomplete anonymous database.

## 1. View supported species

```bash
perl update_mlst_database.pl --list-species
```

View schemes for a species:

```bash
perl update_mlst_database.pl \
  --species "Acinetobacter baumannii" \
  --list-schemes
```

## 2. Update only one species/scheme

Species with one scheme:

```bash
perl update_mlst_database.pl \
  --species "Streptococcus pneumoniae" \
  --threads 4
```

Species with multiple schemes:

```bash
# Oxford
perl update_mlst_database.pl \
  --species "Acinetobacter baumannii" \
  --scheme abaumannii

# Pasteur
perl update_mlst_database.pl \
  --species "Acinetobacter baumannii" \
  --scheme abaumannii_2
```

The update is staged first. The previous database is replaced only after all loci and the profile table download successfully.

## 3. Batch MLST typing

The input directory can contain files such as:

```text
genomes/
├── strain01.fasta
├── strain02.fna
├── strain03.gbk
├── strain04.gbff
└── sample_without_standard_extension
```

Run:

```bash
perl mlst_typing.pl \
  -d genomes \
  -s "Streptococcus pneumoniae" \
  -o MLST_result \
  -t 8 \
  --update T
```

For *A. baumannii* Pasteur MLST:

```bash
perl mlst_typing.pl \
  -d genomes \
  -s "Acinetobacter baumannii" \
  --scheme abaumannii_2 \
  -o MLST_result \
  -t 8 \
  --update T
```

`--update T` updates **only the scheme selected for this run**. Use `--update F` (default) to reuse the local database.

## 4. Calling rules

The allele search uses each downloaded MLST allele as the BLAST query against the complete genome sequence. A known allele is called only when the allele has a **100% identity, full-length exact match**.

- all loci exact + profile exists -> known ST
- all loci exact + profile combination absent -> `Novel`
- one or more loci have no exact database allele -> `Unknown`
- conflicting exact alleles for one locus -> `Unknown` / ambiguous allele

This retains the conservative behavior of the original Pneumo-Typer MLST code while avoiding dependence on gene prediction for FASTA genomes.

## 5. Output

### `ST_out.txt`

Keeps the original compact style:

```text
Strain_name    ST       Profile
sample1.gbk    1        aroE_1-gdh_1-gki_1-recP_1-spi_1-xpt_1-ddl_1
sample2.fasta  Unknown  aroE_1-gdh_#-gki_1-recP_1-spi_1-xpt_1-ddl_1
```

### `MLST_detail.tsv`

Includes input format, selected species/scheme, ST status, profile string, and one allele column per locus.

### `MLST_alleles.tsv`

Long-format table with one row per sample/locus and `exact`, `missing`, or `ambiguous` call status.

### `run_info.tsv`

Records species, scheme, source URI, locus count, input directory and run time.

## Notes on species names

The species labels are retained from the supplied `data_mlst.tar.gz`. Most are Latin names, but some Pasteur entries in that source table are database-level labels rather than full binomials. In those cases, select the scheme directly with `--scheme` if needed.
