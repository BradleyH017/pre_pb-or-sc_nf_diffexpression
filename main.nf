#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.


// Modules to include.
include {
    wf__differential_expression;
} from "./modules/differential_expression.nf"

// Set default parameters.
params.output_dir           = "nf-differential_condition"
params.help                 = false
params.anndata_cell_label = [value: 'cluster']
// Default parameters for differential expression
params.differential_expression = [
    run_process: false
]

// Read a list of target variable names (one per line) from a text file.
def read_target_variable_list(file_path) {
    return file(file_path)
        .readLines()
        .collect { it.trim() }
        .findAll { it != "" && !it.startsWith("#") }
}


// Expand any model whose string fields contain placeholder into one model
// per target_variable, substituting placeholder for that target variable.
// Models that do not contain placeholder are passed through unchanged.
def expand_models_with_target_variables(models, target_variables, placeholder) {
    def expanded = []
    models.each { model ->
        def uses_placeholder = model.any { key, value ->
            value instanceof String && value.contains(placeholder)
        }
        if (uses_placeholder) {
            target_variables.each { target_variable ->
                def new_model = model.collectEntries { key, value ->
                    [key, value instanceof String ? value.replace(placeholder, target_variable) : value]
                }
                expanded.add(new_model)
            }
        } else {
            expanded.add(model)
        }
    }
    return expanded
}


// Read a tab-separated, headerless target_variable -> formula-terms matrix.
// Column 2 may be wrapped in double quotes; these are stripped. Returns a
// list of [target_variable, formula_terms] pairs.
def read_target_x_model_matrix(file_path) {
    return file(file_path)
        .readLines()
        .collect { it.trim() }
        .findAll { it != "" && !it.startsWith("#") }
        .collect { line ->
            def cols = line.split("\t", 2)
            def target_variable = cols[0].trim()
            def formula_terms = cols[1].trim().replaceAll('^"|"$', "")
            return [target_variable, formula_terms]
        }
}


// Like expand_models_with_target_variables, but each target variable gets
// its own formula (taken from target_x_model_rows) instead of substituting
// placeholder into the model's shared formula field.
def expand_models_with_target_x_model_matrix(models, target_x_model_rows, placeholder) {
    def expanded = []
    models.each { model ->
        def uses_placeholder = model.any { key, value ->
            value instanceof String && value.contains(placeholder)
        }
        if (uses_placeholder) {
            target_x_model_rows.each { row ->
                def target_variable = row[0]
                def formula_terms = row[1]
                def new_model = model.collectEntries { key, value ->
                    [key, value instanceof String ? value.replace(placeholder, target_variable) : value]
                }
                new_model.formula = "~ ${formula_terms}"
                expanded.add(new_model)
            }
        } else {
            expanded.add(model)
        }
    }
    return expanded
}


// Define the help messsage.
def help_message() {
    log.info """
    ============================================================================
     single cell differential condition ~ v${VERSION}
    ============================================================================

    Runs basic single cell preprocessing

    Usage:
    nextflow run main.nf -profile <local|lsf> -params-file params.yaml [options]

    Mandatory arguments:
        --file_anndata      Anndata file with cell type labels.

    Optional arguments:
        --output_dir        Directory name to save results to. (Defaults to
                            '${params.output_dir}')

        -params-file        YAML file containing analysis parameters. See
                            example in example_runtime_setup/params.yml.

    Profiles:
        lsf                 lsf cluster execution
    """.stripIndent()
}


// Boot message - either help message or the parameters supplied.
if (params.help){
    help_message()
    exit 0
} else {
    log.info """
    ============================================================================
     single cell differential expression ~ v${VERSION}
    ============================================================================
    file_anndata                  : ${params.file_anndata}
    output_dir (output folder)    : ${params.output_dir}
    """.stripIndent()
}


// Initalize Channels.
// anndata = Channel
//     .fromPath(params.file_anndata)


// Run the workflow.
workflow {
    main:
        // Run differential expression analysis
        if (params.differential_expression.run_process) {
            // If target_x_model_matrix is provided, it takes precedence:
            // expand any model containing target_variable_placeholder into
            // one model per row, each with its own formula. Otherwise, if
            // target_variable_list_file is provided, expand using the
            // shared formula in models.value, one model per listed target
            // variable.
            de_models = params.differential_expression.models
            placeholder = params.differential_expression.target_variable_placeholder.value
            target_x_model_matrix = params.differential_expression.target_x_model_matrix?.value
            target_variable_list_file = params.differential_expression.target_variable_list_file?.value
            if (target_x_model_matrix) {
                target_x_model_rows = read_target_x_model_matrix(target_x_model_matrix)
                de_models = de_models + [
                    value: expand_models_with_target_x_model_matrix(
                        de_models.value,
                        target_x_model_rows,
                        placeholder
                    )
                ]
            } else if (target_variable_list_file) {
                target_variables = read_target_variable_list(target_variable_list_file)
                de_models = de_models + [
                    value: expand_models_with_target_variables(
                        de_models.value,
                        target_variables,
                        placeholder
                    )
                ]
            }

            wf__differential_expression(
                params.output_dir,
                params.file_anndata,
                params.anndata_cell_label.value,
                params.experiment_key_column.value,
                de_models,
                params.differential_expression.de_merge_config,
                params.differential_expression.de_plot_config,
                params.differential_expression.fgsea_config,
                params.differential_expression.run_plots
            )
        }
    // NOTE: One could do publishing in the workflow like so, however
    //       that will not allow one to build the directory structure
    //       depending on the input data call. Therefore, we use publishDir
    //       within a process.
    // publish:
    //     merge_samples.out.anndata to: "${params.output_dir}",
    //         mode: "copy",
    //         overwrite: "true"
}


workflow.onComplete {
    // executed after workflow finishes
    // ------------------------------------------------------------------------
    log.info """\n
    ----------------------------------------------------------------------------
     pipeline execution summary
    ----------------------------------------------------------------------------
    Completed         : ${workflow.complete}
    Duration          : ${workflow.duration}
    Success           : ${workflow.success}
    Work directory    : ${workflow.workDir}
    Exit status       : ${workflow.exitStatus}
    """.stripIndent()
}
