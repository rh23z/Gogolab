# Novel-System Discovery Pipeline

A modular four-step pipeline (HMM search, sequence filtering & clustering, neighborhood extraction, detailed annotation, and reporting) orchestrated by `main.sh`.  

---

## Project Structure

```
config/                     # All YAML configuration files
  ├── config.yaml            # Global config, points to step configs
  ├── step1_config_g2.yaml   # Step 1: HMM search + filters
  ├── step2_config.yaml      # Step 2: Sequence extraction + MMseqs clustering
  ├── step3_config.yaml      # Step 3: Flanking gene extraction + annotation
  └── step4_config.yaml      # Step 4: Reporting
scripts/                     # All step scripts (bash + python)
main.sh                      # Main entrypoint
result/                      # Final outputs
log/                         # Log files
data/                        # Input HMM & FASTA files
```

---

## Dependencies

- Core tools
  - [HMMER](http://hmmer.org/) (`hmmsearch`)
  - [MMseqs2](https://github.com/soedinglab/MMseqs2)
  - [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) (emapper)
  - `GNU parallel`, `yq`
- Python (>=3.8)
  - `pandas`, `biopython`, `tqdm`, `numpy`

---

## Configuration

### Global (`config/config.yaml`)
Points to:
- `root_dir`, `log_dir`, `data_dir`, `scripts_dir`, `result_dir`
- Paths to each step config (`step1_config`, `step2_config`, etc.)

### Step 1 (`step1_config_g2.yaml`)
- `threads`, `parallel_jobs`
- `hmmsearch`: list of {`source`: source tag, `hmm_dir`: path to hmm files (note: the filename must be consistent with the internal HMM name), `fasta_dir`: path to database}.  
  To use different datasets, specify these three parameters in the same format. The `source` name must be unique and clearly distinguishable from other sources.
- `filter.score_cutoff`: score cutoff used in initial filtering
- `function_filter`:  
  - `cov_threshold`: proportion of the hit length relative to the domain length  
  - `min_target_len`: minimum target protein length  
  - `and_filters`: a list of HMM or domain names; when multiple HMMs are provided, the protein must contain all of them simultaneously to pass  
  - `any_filters`: a list of HMM or domain names; if the protein contains any one of them, it passes  
  - `threads`: threads for filtering

### Step 2 (`step2_config.yaml`)
- `threads`, `output_root`
- `clustering`:  
  - `identity`: a list of sequence identity thresholds; clustering is performed iteratively, and each protein will be assigned to clusters at different thresholds  
  - `coverage`: coverage threshold for clustering

### Step 3 (`step3_config.yaml`)
- `threads`, `output_root`
- `s1_extract`: `summary_file`, `upstream`, `downstream`
- `s2_split.parts`
- `s3_detailed_hmm`: `hmm_dir`, `Score_threshold`
- `s5_emapper.emapper_py`

### Step 4 (`step4_config.yaml`)
- `threads`, `output_root`
- `s1_merge_all_info`:  
  - `output_prefix`  
  - `subset_threshold`: a list specifying which sequence identity thresholds should be used for outputting clustering reports; each value must come from the `identity` list in Step 2  
- `s2_generate_report.parts`

---

## Quick Start

1. Prepare configs  
   - Update `config/config.yaml` with your paths (root, log, result, etc.)
   - Edit each step’s YAML to point to your HMM/FASTA directories and set thresholds.

2. Run pipeline
   ```bash
   bash main.sh
   ```

3. Monitor logs (all under `./log/`)  
   Example: `step1_hmmsearch_<source>.log`, `step2_mmseqs_cluster.log`, `step3_s1_extract.log`

---

## Workflow Summary

1. Step 1: HMM Search & Filtering  
   - Run `hmmsearch` in parallel on multiple FASTA sources  
   - Integrate domtblout results  
   - Apply score cutoff and functional filters  

2. Step 2: Sequence Extraction & Clustering  
   - Extract sequences from FASTA using filtered IDs  
   - Perform iterative MMseqs2 clustering at multiple identity thresholds  

3. Step 3: Neighborhood Analysis & Annotation  
   - Extract upstream/downstream gene neighborhoods  
   - Merge/split FASTA for batch processing  
   - Run detailed HMM search and eggNOG-mapper  
   - Integrate PFAM/EGG annotations  

4. Step 4: Reporting  
   - Merge clustering and annotation results  
   - Generate split reports (flanking and target genes)  

---

## Outputs

- Step 1 → `step1_filtered_df.tsv`, `step1_function_filtered_df.tsv`
- Step 2 → `step2_function_filtered_seqs.fasta`, MMseqs clustering results
- Step 3 → `step3_df.tsv`, `step3_merged_pfam.tsv`, `step3_merged_egg.tsv`, `step3_df_position.tsv`
- Step 4 → Final merged results in `result/`, with split reports for target/flanking genes.  
  You can also specify which clustering thresholds to output reports for. Since final reports are usually large, they will be split into multiple files in the corresponding directories. Final outputs include annotations of target genes and their flanking genes.

---

## Tips

- Increase or reduce `threads` and `parallel_jobs` according to your machine resources  
- Adjust `filter.score_cutoff` and `function_filter` criteria to change filtering stringency  
- Use `s2_split.parts` (Step 3) to handle very large FASTA datasets  
- Logs are the first place to check for debugging
