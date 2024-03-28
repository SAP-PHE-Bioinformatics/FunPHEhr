/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowFunphehr.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db, ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check krakendb
if(! params.skip_kraken2){
    if(params.kraken2db){
        kraken2db = file(params.kraken2db)
    } else {
        exit 1, "Missing Kraken2 DB arg"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { MEDAKA                    } from '../modules/local/medaka'
include { KRAKEN2_DB_PREPARATION    } from '../modules/local/kraken2_db_preparation'
include { BUSCO       } from '../modules/local/busco/main'
include { BUSCO as BUSCO_AUTO  } from '../modules/local/busco/main'
include { ITSX                                  } from '../modules/local/itsx'
include { HELIXER                               } from '../modules/local/helixer'
include { AGAT_EXTRACTSEQUENCES as GFF2PROTEIN                           } from '../modules/local/agat_extract_protein'
include { AGAT_MANAGEFUNCTIONALANNOTATION as MERGE_FUNC_ANNOTATION                 } from '../modules/local/agat_merge_func_annotation'
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SEQKIT_STATS                         } from '../modules/nf-core/seqkit/stats/main'
include { NANOPLOT                              } from '../modules/nf-core/nanoplot/main'
include { PORECHOP_PORECHOP                     } from '../modules/nf-core/porechop/porechop/main'
include { MINIMAP2_ALIGN                        } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_CONSENSUS  } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_POLISH     } from '../modules/nf-core/minimap2/align/main'
include { MINIASM                               } from '../modules/nf-core/miniasm/main'
include { FLYE                                  } from '../modules/nf-core/flye/main'
include { RACON                                 } from '../modules/nf-core/racon/main'
include { SAMTOOLS_SORT                         } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                        } from '../modules/nf-core/samtools/index/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_READS      } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_ASSEMBLY   } from '../modules/nf-core/kraken2/kraken2/main'
include { QUAST                                 } from '../modules/nf-core/quast/main'
include { GUNZIP                                } from '../modules/nf-core/gunzip/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { INTERPROSCAN                          } from '../modules/nf-core/interproscan/main'  
include { BLAST_BLASTP                          } from '../modules/nf-core/blast/blastp/main'                                                      
include { CHOPPER                               } from '../modules/nf-core/chopper/main'
//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
//include { struc_func_annotation               } from '../subworkflows/local/struc_func_annotation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow FUNPHEHR {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    def criteria = multiMapCriteria {
        meta, longFastq ->
            longreads: longFastq   != 'NA' ? tuple(meta, longFastq)         : null
    }
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    Channel
        .fromSamplesheet('input')
        .multiMap (criteria)
        .set { ch_input }
    // reconfigure channels
    ch_input
        .longreads
        .filter{ it != null }
        .set { ch_longreads } 
    //
    // MODULE: Nanoplot, quality check for nanopore reads and Quality/Length Plots
    //
    //ch_input.view()
    CHOPPER (
        ch_longreads
    )
    ch_versions= ch_versions.mix(CHOPPER.out.versions.ifEmpty(null))    

    ch_porechop_log_multiqc = Channel.empty()
    
    PORECHOP_PORECHOP (
        CHOPPER.out.fastq
        )
        ch_porechop_log_multiqc = PORECHOP_PORECHOP.out.log
        ch_versions = ch_versions.mix( PORECHOP_PORECHOP.out.versions.ifEmpty(null) )
    
    ch_assembly = Channel.empty()
    
    ch_for_assembly = PORECHOP_PORECHOP.out.reads
    

    
    NANOPLOT (
        PORECHOP_PORECHOP.out.reads
    )
    ch_nanoplot_txt_multiqc = NANOPLOT.out.txt
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.ifEmpty(null))


    
    // update ch_for_assmebly so if meta.id contains NEG exclude from channel_for_assembly
    ch_for_assembly = ch_for_assembly.filter{ meta, lr -> !meta.id.contains("NEG") }

    //
    // ASSEMBLY: Miniasm
    //
    ch_assembly = Channel.empty()

    //
    // MODULE: Miniasm, genome assembly, long reads
    if ( params.assembler == 'miniasm' ) {
        MINIMAP2_ALIGN (
            ch_for_assembly.map{ meta,lr -> tuple(meta,lr) },
            [],
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.ifEmpty(null))

        ch_for_assembly
            .join(MINIMAP2_ALIGN.out.paf)
            .map { meta, lr, paf-> tuple(meta, lr, paf) }
            .set { ch_for_miniasm }

        MINIASM (
            ch_for_miniasm
        )
        ch_versions = ch_versions.mix(MINIASM.out.versions.ifEmpty(null))

        MINIMAP2_CONSENSUS (
            ch_for_assembly.map{ meta,sr,lr -> tuple(meta,lr) },
            MINIASM.out.assembly.map { meta, assembly -> assembly },
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_CONSENSUS.out.versions.ifEmpty(null))

        ch_for_assembly
            .join(MINIASM.out.assembly)
            .join(MINIMAP2_CONSENSUS.out.paf)
            .map { meta, lr, assembly, paf -> tuple(meta, lr, assembly, paf) }
            .set{ ch_for_racon }

        RACON (
            ch_for_racon
        )
        ch_assembly = ch_assembly.mix( RACON.out.improved_assembly.dump(tag: 'miniasm') )
        ch_versions = ch_versions.mix( RACON.out.versions.ifEmpty(null) )
    }

    //
    // MODULE: flye, genome assembly of long reads. inputs long read channel ch_for_assembly also needs mode "--nano-raw"
    //
    if( params.assembler == 'flye' ){
        FLYE (
            // take only meta and lr from meta, [], lr
            ch_for_assembly.map{ meta,lr -> tuple(meta,lr) },
            '--nano-hq'
        )

        ch_assembly = ch_assembly.mix( FLYE.out.fasta.dump(tag: 'flye') )
        ch_versions = ch_versions.mix( FLYE.out.versions.ifEmpty(null) )
    }
    //
    // MODULE: Medaka, polishes assembly - should take either miniasm or flye consensus sequence
    //
    if ( !params.skip_polish && params.polish_method == 'medaka' ) {
        ch_for_assembly
            .join( ch_assembly )
            .map { meta, lr, assembly -> tuple(meta, lr, assembly) }
            .set { ch_for_medaka }

        MEDAKA ( ch_for_medaka.dump(tag: 'into_medaka') )
        ch_versions = ch_versions.mix(MEDAKA.out.versions.ifEmpty(null))
    }

    
    //
    // MODULE: QUAST, assembly QC
    //
    ch_assembly
        .collect{ it[1] }
        .map { consensus_collect -> tuple([id: "report"], consensus_collect) }
        .set { ch_to_quast }

    QUAST (
        ch_to_quast,
        [[:],[]],
        [[:],[]]
    )
    ch_quast_multiqc = QUAST.out.tsv
    ch_versions      = ch_versions.mix(QUAST.out.versions.ifEmpty(null))

    // Check assemblies that require further processing for gene annotation
    ch_assembly
        .branch{ meta, fasta ->
            gzip: fasta.name.endsWith('.gz')
            skip: true
        }
        .set{ ch_assembly_for_gunzip }
    

    ch_assembly 
        .map{ meta, fasta -> tuple(meta, fasta) }
        .set{ ch_assembly_for_busco }
    //ch_assembly_for_busco.view()

    BUSCO (
        ch_assembly_for_busco,
        'genome',
        'fungi_odb10',
        params.busco_db_path
    )
    ch_busco_fungi_multiqc = BUSCO.out.batch_summary

       // MODULE: Kraken2, QC for sample purity
    //
    ch_kraken_reads_multiqc  = Channel.empty()
    ch_kraken2db = Channel.empty()
    if ( !params.skip_kraken2 ) {
    // if path for directory of kraken2db is valid/exists use it to create ch_kraken2db and use the last directory as db name, else run kraken2_db_preparation with "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz"
        if (params.kraken2db){
            // create channel of info, db where db is the path of directory and info is the last directory name eg /scratch/kraken/bigDB will be  bigDB, /scratch/kraken/bigDB
            // channel should look like tuple val("${pathname}"), path("db")
            ch_kraken2db = Channel.fromPath(params.kraken2db, type: 'dir', checkIfExists: true)
                .map{ db -> tuple(db.name, db) }
            
        } else {
            KRAKEN2_DB_PREPARATION (
                'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz'
            )
            ch_kraken2db = KRAKEN2_DB_PREPARATION.out.kraken2db
        }
       
    

        KRAKEN2_READS(
            PORECHOP_PORECHOP.out.reads
            .map { meta, reads -> 
                info = [:]
                info.id = meta.id
                info.single_end = true
                [info, reads]
            }
            .dump(tag: 'kraken2_reads'),
            ch_kraken2db.map{ info, db -> db }.dump(tag: 'kraken2_db_preparation').first(),
            false,
            false
        )
        ch_kraken_reads_multiqc = KRAKEN2_READS.out.report
        ch_versions = ch_versions.mix(KRAKEN2_READS.out.versions.ifEmpty(null))
    }

    // BUSCO_AUTO (
    //     ch_assembly_for_busco,
    //     'genome',
    //     'auto',
    //     params.busco_db_path
    // )
    // ch_busco_auto_summary = BUSCO_AUTO.out.batch_summary
    // ch_versions      = ch_versions.mix(BUSCO.out.versions.ifEmpty(null))
   

    ch_assembly
        .map{ meta, fasta -> tuple(meta, fasta) }
        .set{ ch_assembly_for_itsx }
    
    //ch_assembly_for_itsx.view()
    ITSX (
        ch_assembly_for_itsx
    )
    ch_versions = ch_versions.mix(ITSX.out.versions.ifEmpty(null))

    // //need to blast the ITSX sequences to the UNITE database
    // BLAST_BLASTN
    


    //ch_assembly.view()
    ch_kraken_assembly_multiqc = Channel.empty()
    KRAKEN2_ASSEMBLY(
            ch_assembly
            .map {meta, fasta ->
                info=[:]
                info.id = meta.id
                info.single_end = true
                [info, fasta]
            }
            .dump(tag: 'kraken2_assembly'),
            ch_kraken2db.map{ info, db -> db }.dump(tag: 'kraken2_db_preparation').first(),
            false,
            false
        )
        ch_kraken_assembly_multiqc = KRAKEN2_ASSEMBLY.out.report

    // annotation if annotation=true call subworkflow annotation
    if ( !params.skip_annotation ) {

        ch_assembly
            .map{ meta, fasta -> tuple(meta, fasta) }
            .set{ ch_assembly_for_annotation }

    if (params.cuda) {
     HELIXER (
        ch_assembly_for_annotation
    )

    ch_versions = ch_versions.mix(HELIXER.out.versions.ifEmpty(null))

     GFF2PROTEIN( 
        HELIXER.out.gff3, 
        ch_assembly_for_annotation
    )
    ch_versions = ch_versions.mix(GFF2PROTEIN.out.versions.ifEmpty(null))

    BLAST_BLASTP (
        GFF2PROTEIN.out.proteins.splitFasta(by: 1000, file: true),
        params.blast_db_directory,
        "6"
    )

    INTERPROSCAN (
        GFF2PROTEIN.out.proteins.splitFasta(by: 1000, file: true),
        params.interproscan_db_directory
    )
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions.ifEmpty(null))

    MERGE_FUNC_ANNOTATION (
        HELIXER.out.gff3,
        INTERPROSCAN.out.tsv,
        BLAST_BLASTP.out.tsv,
        params.blast_db_directory
    )
    ch_versions = ch_versions.mix(MERGE_FUNC_ANNOTATION.out.versions.ifEmpty(null))
    }
    else {
        FUNANNOTATE_MASK(
            ch_assembly_for_annotation
        )
    
    ch_versions = ch_versions.mix(FUNANNOTATE_MASK.out.versions.ifEmpty(null))

    FUNANNOTATE_PREDICT(
        FUNANNOTATE_MASK.out.masked_fasta
    )


    }
    }
    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc){
        workflow_summary        = WorkflowFunphehr.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary     = Channel.value(workflow_summary)
        methods_description     = WorkflowFunphehr.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description  = Channel.value(methods_description)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(ch_kraken_assembly_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_kraken_reads_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_quast_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_nanoplot_txt_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_porechop_log_multiqc.collect{it[1]}.ifEmpty([]))
        //ch_multiqc_files = ch_multiqc_files.mix(ch_busco_auto_summary.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_busco_fungi_multiqc.collect{it[1]}.ifEmpty([]))
        //ch_multiqc_files = ch_multiqc_files.mix(ch_kraken_assembly_multiqc.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
