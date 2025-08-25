#!/bin/bash

# ========== Read global config path ==========
GLOBAL_CONFIG="config/config.yaml"

ROOT_DIR=$(cat "$GLOBAL_CONFIG" | yq -r '.root_dir')
LOG_DIR=$(cat "$GLOBAL_CONFIG" | yq -r '.log_dir')
DATA_DIR=$(cat "$GLOBAL_CONFIG" | yq -r '.data_dir')
RESULT_DIR=$(cat "$GLOBAL_CONFIG" | yq -r '.result_dir')
SCRIPTS_DIR=$(cat "$GLOBAL_CONFIG" | yq -r '.scripts_dir')
STEP1_CONFIG=$(cat "$GLOBAL_CONFIG" | yq -r '.step1_config')
STEP2_CONFIG=$(cat "$GLOBAL_CONFIG" | yq -r '.step2_config')
STEP3_CONFIG=$(cat "$GLOBAL_CONFIG" | yq -r '.step3_config')
STEP4_CONFIG=$(cat "$GLOBAL_CONFIG" | yq -r '.step4_config')

# ========== Extract parameters from step1_config ==========
STEP1_THREADS=$(cat "$STEP1_CONFIG" | yq -r '.threads')
STEP1_JOBS=$(cat "$STEP1_CONFIG" | yq -r '.parallel_jobs')
STEP1_OUTPUT=$(cat "$STEP1_CONFIG" | yq -r '.output_root')
STEP1_HMMSEARCH_COUNT=$(cat "$STEP1_CONFIG" | yq -r '.hmmsearch | length')

mkdir -p "$STEP1_OUTPUT" "$LOG_DIR"

# ========== Step 1.1: HMMSEARCH loop ==========
for ((i = 0; i < STEP1_HMMSEARCH_COUNT; i++)); do
  STEP1_SOURCE=$(cat "$STEP1_CONFIG" | yq -r ".hmmsearch[$i].source")
  STEP1_HMM_DIR=$(cat "$STEP1_CONFIG" | yq -r ".hmmsearch[$i].hmm_dir")
  STEP1_FASTA_DIR=$(cat "$STEP1_CONFIG" | yq -r ".hmmsearch[$i].fasta_dir")
  STEP1_OUTPUT_DIR="$STEP1_OUTPUT/$STEP1_SOURCE"
  mkdir -p "$STEP1_OUTPUT_DIR"

  echo "Source: $STEP1_SOURCE | HMM_DIR: $STEP1_HMM_DIR | FASTA_DIR: $STEP1_FASTA_DIR -> OUTPUT: $STEP1_OUTPUT_DIR"
  bash "$SCRIPTS_DIR/step1_s1_batch_hmmsearch_parallel.sh" \
       "$STEP1_HMM_DIR" "$STEP1_FASTA_DIR" "$STEP1_OUTPUT_DIR" "$STEP1_THREADS" "$STEP1_JOBS" \
       2>&1 | tee "$LOG_DIR/step1_hmmsearch_${STEP1_SOURCE}.log"
done

echo "Step1 finished"


# ========== Step 1.2: Integrate domtblout results ==========
echo "Step1.2: Integrating domtblout results"

python "$SCRIPTS_DIR/step1_s2_integrate_hmmsearch_domtblout.py" \
  --input_dir "$STEP1_OUTPUT" \
  --output_dir "$STEP1_OUTPUT" \
  --threads "$STEP1_THREADS" \
  2>&1 | tee "$LOG_DIR/step1_integrate.log"

echo "Step1.2 merge completed"

# ========== Step 1.3: Filter & simplify ==========
echo "Step1.3: Filter & simplify"

STEP1_SCORE_CUTOFF=$(cat "$STEP1_CONFIG" | yq -r '.filter.score_cutoff')

python "$SCRIPTS_DIR/step1_s3_filter_simplify.py" \
  --input_dir "$STEP1_OUTPUT" \
  --output_prefix "$STEP1_OUTPUT/step1_filtered_df" \
  --score_cutoff "$STEP1_SCORE_CUTOFF" \
  2>&1 | tee "$LOG_DIR/step1_filter.log"

echo "Step1.3 filtering completed"

# ========== Step 1.4: Functional filtering ==========
echo "Step1.4: Functional domain-based filtering"

STEP1_FUNC_COV=$(cat "$STEP1_CONFIG" | yq -r '.function_filter.cov_threshold')
STEP1_FUNC_LEN=$(cat "$STEP1_CONFIG" | yq -r '.function_filter.min_target_len')
STEP1_FUNC_THREADS=$(cat "$STEP1_CONFIG" | yq -r '.function_filter.threads')
STEP1_FUNC_AND=$(cat "$STEP1_CONFIG" | yq -r '.function_filter.and_filters | join(" ")')
STEP1_FUNC_ANY=$(cat "$STEP1_CONFIG" | yq -r '.function_filter.any_filters | join(" ")')

python "$SCRIPTS_DIR/step1_s4_function_filter_parallel.py" \
  --input_file "$STEP1_OUTPUT/step1_filtered_df.tsv" \
  --output_file "$STEP1_OUTPUT/step1_function_filtered_df.tsv" \
  --cov_threshold "$STEP1_FUNC_COV" \
  --min_target_len "$STEP1_FUNC_LEN" \
  --and_filters $STEP1_FUNC_AND \
  --any_filters $STEP1_FUNC_ANY \
  --threads "$STEP1_FUNC_THREADS" \
  2>&1 | tee "$LOG_DIR/step1_function_filter.log"

echo "Step1.4 functional filtering completed"

# ========== Extract parameters from step2_config ==========
STEP2_THREADS=$(cat "$STEP2_CONFIG" | yq -r '.threads')
STEP2_OUTPUT_ROOT=$(cat "$STEP2_CONFIG" | yq -r '.output_root')

mkdir -p "$STEP2_OUTPUT_ROOT" "$LOG_DIR"

# ========== Step 2.1: Extract and filter sequences ==========
echo "Step2.1: Extract target sequences and filter start"

python "$SCRIPTS_DIR/step2_s1_extract_seqs.py" \
  --input "$STEP1_OUTPUT/step1_function_filtered_df.tsv" \
  --output_prefix "$STEP2_OUTPUT_ROOT/step2_function_filtered_seqs" \
  --threads "$STEP2_THREADS" \
  --log_dir "$LOG_DIR" \
  2>&1 | tee "$LOG_DIR/step2_s1_extract_seqs.log"

echo "Step2.1 completed"

# ========== Step 2.2: MMseqs2 clustering ==========
echo "Step2.2: Iterative MMseqs2 clustering"

MMSEQS_THREADS=$(cat "$STEP2_CONFIG" | yq -r '.threads')
MMSEQS_COV=$(cat "$STEP2_CONFIG" | yq -r '.clustering.coverage')
MMSEQS_ID_LIST=$(cat "$STEP2_CONFIG" | yq -r '.clustering.identity | join(" ")')
MMSEQS_TMP_DIR="$STEP2_OUTPUT_ROOT/tmp"
MMSEQS_OUTPUT_DIR="$STEP2_OUTPUT_ROOT/mmseqs_output"

mkdir -p "$MMSEQS_OUTPUT_DIR" "$MMSEQS_TMP_DIR"

python "$SCRIPTS_DIR/step2_s2_mmseqs_iterative.py" \
  "$STEP2_OUTPUT_ROOT/step2_function_filtered_seqs.fasta" \
  "$MMSEQS_OUTPUT_DIR" \
  "$MMSEQS_TMP_DIR" \
  --identity $MMSEQS_ID_LIST \
  --coverage "$MMSEQS_COV" \
  --threads "$MMSEQS_THREADS" \
  2>&1 | tee "$LOG_DIR/step2_mmseqs_cluster.log"

echo "Step2.2 MMseqs clustering completed"


# ========== Extract parameters from step3_config ==========
STEP3_THREADS=$(cat "$STEP3_CONFIG" | yq -r '.threads')
STEP3_OUTPUT_ROOT=$(cat "$STEP3_CONFIG" | yq -r '.output_root')

# S1: flanking extraction
STEP3_S1_SUMMARY_FILE=$(cat "$STEP3_CONFIG" | yq -r '.s1_extract.summary_file')
STEP3_S1_UPSTREAM=$(cat "$STEP3_CONFIG" | yq -r '.s1_extract.upstream')
STEP3_S1_DOWNSTREAM=$(cat "$STEP3_CONFIG" | yq -r '.s1_extract.downstream')

# S2: merge & split FASTA
STEP3_S2_SPLIT_PARTS=$(cat "$STEP3_CONFIG" | yq -r '.s2_split.parts')

# S3: detailed HMM search (support full-seq threshold)
STEP3_S3_HMM_DIR=$(cat "$STEP3_CONFIG" | yq -r '.s3_detailed_hmm.hmm_dir')
STEP3_S3_SCORE_THRESHOLD=$(cat "$STEP3_CONFIG" | yq -r '.s3_detailed_hmm.Score_threshold')

# S5: emapper
STEP3_S5_EMAPPER_PY=$(cat "$STEP3_CONFIG" | yq -r '.s5_emapper.emapper_py')

# Create directories
mkdir -p "$STEP3_OUTPUT_ROOT" \
         "$STEP3_OUTPUT_ROOT/flanking_fasta" \
         "$STEP3_OUTPUT_ROOT/flanking_split" \
         "$STEP3_OUTPUT_ROOT/flanking_hmm" \
         "$STEP3_OUTPUT_ROOT/flanking_emapper" \
         "$LOG_DIR"

# ========== Step 3.1: Flanking extraction ==========
echo "Step3.1: Flanking sequence extraction start"
python "$SCRIPTS_DIR/step3_s1_extract_flanking_genes_parallel.py" \
  -i "$STEP3_S1_SUMMARY_FILE" \
  -o "$STEP3_OUTPUT_ROOT/flanking_fasta" \
  --df-out "$STEP3_OUTPUT_ROOT/step3_df.tsv" \
  --workers "$STEP3_THREADS" \
  --upstream "$STEP3_S1_UPSTREAM" --downstream "$STEP3_S1_DOWNSTREAM" \
  2>&1 | tee "$LOG_DIR/step3_s1_extract.log"
echo "Step3.1 completed"

# ========== Step 3.2: Merge & split FASTA ==========
echo "Step3.2: Merge & split FASTA"
bash "$SCRIPTS_DIR/step3_s2_split.sh" \
  "$STEP3_OUTPUT_ROOT/flanking_fasta" \
  "$STEP3_OUTPUT_ROOT/flanking_split" \
  "$STEP3_S2_SPLIT_PARTS" \
  step3_flanking_prot_merged.fasta \
  2>&1 | tee "$LOG_DIR/step3_s2_split.log"
echo "Step3.2 completed"
rm -r "$STEP3_OUTPUT_ROOT/flanking_fasta"

# ========== Step 3.3: Detailed HMM search ==========
echo "Step3.3: Detailed HMM search"
bash "$SCRIPTS_DIR/step3_s3_detailed_hmm.sh" \
  "$STEP3_S3_HMM_DIR" \
  "$STEP3_OUTPUT_ROOT/flanking_split" \
  "$STEP3_OUTPUT_ROOT/flanking_hmm" \
  "$STEP3_S2_SPLIT_PARTS" \
  "$STEP3_S3_SCORE_THRESHOLD" \
  2>&1 | tee "$LOG_DIR/step3_s3_hmm.log"
echo "Step3.3 completed"

# ========== Step 3.4: Integrate domtblout ==========
echo "Step3.4: Integrate domtblout"
python "$SCRIPTS_DIR/step3_s4_detailed_hmm_inte.py" \
  -i "$STEP3_OUTPUT_ROOT/flanking_hmm" \
  -o "$STEP3_OUTPUT_ROOT/step3_merged_pfam.tsv" \
  --workers "$STEP3_THREADS" \
  --recursive --include-gz --relative-path \
  2>&1 | tee "$LOG_DIR/step3_s4_inte.log"
echo "Step3.4 completed"

# ========== Step 3.5: emapper ==========
echo "Step3.5: emapper batch annotation"
bash "$SCRIPTS_DIR/step3_s5_emmaper.sh" \
  -i "$STEP3_OUTPUT_ROOT/flanking_split" \
  -o "$STEP3_OUTPUT_ROOT/flanking_emapper" \
  -t "$STEP3_THREADS" \
  -j "$STEP3_S2_SPLIT_PARTS" \
  --emapper "$STEP3_S5_EMAPPER_PY" \
  2>&1 | tee "$LOG_DIR/step3_s5_emapper.log"
echo "Step3.5 completed"

# ========== Step 3.6: Merge emapper results ==========
echo "Step3.6: Merge emapper results"
python "$SCRIPTS_DIR/step3_s6_emmapper_merge.py" \
  -i "$STEP3_OUTPUT_ROOT/flanking_emapper" \
  -o "$STEP3_OUTPUT_ROOT/step3_merged_egg.tsv" \
  --recursive --workers "$STEP3_THREADS" \
  2>&1 | tee "$LOG_DIR/step3_s6_merge.log"
echo "Step3.6 completed"

# ========== Step 3.7: (Optional) Map annotations back to flanking table ==========
echo "Step3.7: Integrate PFAM/EGG annotations into flanking table"
python "$SCRIPTS_DIR/step3_s7_integrate_flanking_parallel.py" \
  --pfam "$STEP3_OUTPUT_ROOT/step3_merged_pfam.tsv" \
  --emapper "$STEP3_OUTPUT_ROOT/step3_merged_egg.tsv" \
  --step3-df "$STEP3_OUTPUT_ROOT/step3_df.tsv" \
  -o "$STEP3_OUTPUT_ROOT/step3_df_position.tsv" \
  --workers "$STEP3_THREADS" --chunksize 5000 \
  2>&1 | tee "$LOG_DIR/step3_s7_integrate.log"
echo "Step3.7 completed"

echo "Step3 all done -> $STEP3_OUTPUT_ROOT"


# ========== Extract parameters from step4_config ==========
STEP4_THREADS=$(cat "$STEP4_CONFIG" | yq -r '.threads')
STEP4_OUTPUT_ROOT=$(cat "$STEP4_CONFIG" | yq -r '.output_root')
STEP4_OUTPUT_PREFIX=$(cat "$STEP4_CONFIG" | yq -r '.s1_merge_all_info.output_prefix')
STEP4_SELECT_IDENTITY=$(cat "$STEP4_CONFIG" | yq -r '.s1_merge_all_info.subset_threshold | join(" ")')
STEP4_REPORT_PARTS=$(cat "$STEP4_CONFIG" | yq -r '.s2_generate_report.parts')

# ========== Step4.1: Merge clustering info ==========
echo "Step4.1: Merge clustering info start"
python "$SCRIPTS_DIR/step4_s1_merge_all_info.py" \
  --input_df "$STEP3_OUTPUT_ROOT/step3_df.tsv" \
  --cluster_dir "$MMSEQS_OUTPUT_DIR" \
  --identity_thresholds $MMSEQS_ID_LIST \
  --select_identity $STEP4_SELECT_IDENTITY \
  --output_prefix "$RESULT_DIR/$STEP4_OUTPUT_PREFIX" \
  --input_fasta "$STEP2_OUTPUT_ROOT/step2_function_filtered_seqs.fasta" \
  2>&1 | tee "$LOG_DIR/step4_s1_merge_all_info.log"
echo "Step4.1 completed"

# ========== Step4.2: Generate split reports ==========
echo "Step4.2: Generate flanking gene annotation report"
python "$SCRIPTS_DIR/step4_s2_generate_report.py" \
  --input "$STEP3_OUTPUT_ROOT/step3_df_position.tsv" \
  --num_parts "$STEP4_REPORT_PARTS" \
  --output_dir "$RESULT_DIR/flanking_gene_annotations_splitted" \
  --prefix flanking_gene_annotations \
  2>&1 | tee "$LOG_DIR/step4_s2_flanking_report.log"

echo "Step4.2: Generate target gene annotation report"
python "$SCRIPTS_DIR/step4_s2_generate_report.py" \
  --input "$RESULT_DIR/${STEP4_OUTPUT_PREFIX}_result.tsv" \
  --num_parts "$STEP4_REPORT_PARTS" \
  --output_dir "$RESULT_DIR/target_gene_annotations_splitted" \
  --prefix target_gene_annotations \
  2>&1 | tee "$LOG_DIR/step4_s2_target_report.log"

echo "Step4 completed"
