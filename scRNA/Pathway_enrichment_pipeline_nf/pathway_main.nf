// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Define base directory and subdirectories
def outdir = "${PWD}/MAST_pathway_${params.cond1}vs${params.cond2}"
def ULM_tables = "${outdir}/ULM_tables"
def ULM_pathway_analysis_heatmap = "${outdir}/plots/ULM_pathway_analysis_heatmap"

// Process to ensure all directories exist
process createDirectories {
    executor = 'local'
    clusterOptions = '--ntasks=1 --mem=1Gb --time=00:05:00'
    output:
    val(true)
    script:
    """
    mkdir -p "${outdir}"
    mkdir -p "${ULM_tables}"
    mkdir -p "${ULM_pathway_analysis_heatmap}"
    """
}

// Process to run the ULM calculation
process ULM_calcuation {
    executor = 'slurm'
    clusterOptions = '--ntasks=1 --mem=45Gb --time=24:00:00'
    def python_path_ulm = "${PWD}/scripts/ULM_MSigDB_analysis_nf.py"
    publishDir ULM_tables, mode: 'copy' 
    input:
    val(dummy)
    output:
    val(true)
    script:
    """
    python ${python_path_ulm} --input_dir_1 "${params.inputdir}" --outdir_3 "${ULM_tables}" --cond1 ${params.cond1} --cond2 ${params.cond2} --species "${params.species}"
    """
}

// Process to polish the ULM tables
process POLISH_ULM_TABLES {
    executor = 'slurm'
    clusterOptions = '--ntasks=1 --mem=10Gb --time=02:00:00'
    def rscript_polish_path = "${PWD}/scripts/polish_ULM_tables_nf.R"
    publishDir ULM_tables, mode: 'copy', overwrite: true, pattern: "combined_*"
    input:
    val(dummy)
    val ulm_tables_dir
    val cond1
    val cond2
    path gene_annotation_file
    output:
    file("combined_*")
    script:
    """
    Rscript ${rscript_polish_path} --ulm_tables_dir "${ulm_tables_dir}" --csv_path_base "${params.inputdir}"  --gene_annotation "${gene_annotation_file}" --cond1 "${cond1}" --cond2 "${cond2}"
    """
}

// Process to run the ULM plotting
process ULM_plotting {
    executor = 'slurm'
    clusterOptions = '--ntasks=1 --mem=10Gb --time=02:00:00'
    def rscript_5_path = "${PWD}/scripts/combined_heatmap_ULM_nf_test.R"
    publishDir ULM_pathway_analysis_heatmap, mode: 'copy', overwrite: true, pattern: "*.pdf"
    input:
    file(csvs)
    output:
    file("*.pdf")
    script:
    """
    Rscript ${rscript_5_path} --input_dir_6 "${ULM_tables}" --outdir_3 "${ULM_pathway_analysis_heatmap}"
    """
}

workflow {
    // Step 1: Create directories
    create_done = createDirectories()
    
    // Step 2: ULM calculation, waits for directories
    ulm_done = ULM_calcuation(create_done)
    
    // Step 3: Polish tables, waits for ULM calculation
    def gene_annotation_path = (params.species == 'human') ? 
        "/data/iu627335/nf-tutorial/MSigDB_human_database.csv" : 
        "/data/iu627335/nf-tutorial/Msg_mouse_database.csv"
    
    gene_annotation_file = file(gene_annotation_path)
    polish_outputs = POLISH_ULM_TABLES(ulm_done, ULM_tables, params.cond1, params.cond2, gene_annotation_file)
    
    // Step 4: ULM plotting, waits for polish outputs
    ULM_plotting(polish_outputs)
}

