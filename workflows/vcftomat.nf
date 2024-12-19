/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { GATK4_GENOTYPEGVCFS    } from '../modules/nf-core/gatk4/genotypegvcfs/main'
include { BCFTOOLS_CONCAT        } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_REHEADER      } from '../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_MERGE         } from '../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX            } from '../modules/nf-core/tabix/tabix/main'
include { VCF2MAT                } from '../modules/local/vcf2mat/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_vcftomat_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCFTOMAT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta
    fai
    dict

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // add index to non-indexed VCFs
    //
    (ch_has_index, ch_has_no_index) = ch_samplesheet
        .map{ it ->
            def name = it[1][0].baseName
                name = name
                    .replaceFirst(/\.g\.vcf$/, "")
                    .replaceFirst(/\.genome\.vcf$/, "")
                    .replaceFirst(/\.genome\.g\.vcf$/, "")
                    .replaceFirst(/\.g$/, "")
                    .replaceFirst(/\.genome$/, "")
                    .replaceFirst(/\.vcf$/, "")
            [ it[0] + [ name:name ], it[1] ]
        }
        .branch{
                has_index: !it[0].to_index
                to_index: it[0].to_index
        }

    TABIX_TABIX( ch_has_no_index )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    ch_indexed = ch_has_no_index.join(
        TABIX_TABIX.out.tbi
            .map{ it -> [ it[0], [it[1]] ] }
        ).map { meta, vcf, tbi -> [ meta, [ vcf[0], tbi[0] ] ] }

    // Join both channels back together
    ch_vcf_tbi = ch_has_index.mix(ch_indexed)

    //
    // Convert gvcfs to vcfs
    //
    (ch_gvcf, ch_normal_vcf) = ch_vcf_tbi.branch {
            gvcf: it[0].gvcf
            vcf: !it[0].gvcf
        }

    GATK4_GENOTYPEGVCFS(
        ch_gvcf.map{ it -> [ it[0], it[1][0], it[1][1], [], [] ] },
        fasta.map{ it -> [ [ id:it.baseName ], it ] },
        fai.map{ it -> [ [ id:it.baseName ], it ] },
        dict.map{ it -> [ [ id:it.baseName ], it ] },
        [[],[]], // dbsnp
        [[],[]] // dbsnp_tbi
    )

    ch_vcf_index = GATK4_GENOTYPEGVCFS.out.vcf
            .join(GATK4_GENOTYPEGVCFS.out.tbi)
            .map { meta, vcf, tbi -> [ meta, [ vcf, tbi ] ] }

    ch_vcf = ch_normal_vcf.mix(ch_vcf_index)

    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    //
    // Concatenate converted VCFs if the entries for "id" and "label" are the same
    //
    (ch_single_vcf, ch_multiple_vcf) = ch_vcf
        .map { meta, files ->
            // Assuming files is a list of all VCF and TBI files
            def vcfs = files.findAll { it.name.endsWith('.vcf.gz') }
            def tbis = files.findAll { it.name.endsWith('.vcf.gz.tbi') }
            [ [meta.id, meta.label], meta, vcfs, tbis]
        }
        .groupTuple(by: 0)
        .map { id_label, metas, vcfs, tbis ->
            def meta = metas[0]
            def vcf_count = vcfs.flatten().size()
            meta.single_vcf = (vcf_count == 1)
            [meta, vcfs.flatten(), tbis.flatten()]
        }.branch {
            single: it[0].single_vcf
            multiple: !it[0].single_vcf
        }

    BCFTOOLS_CONCAT( ch_multiple_vcf )

    ch_vcf_index = BCFTOOLS_CONCAT.out.vcf
            .join(BCFTOOLS_CONCAT.out.tbi)

    ch_vcf_concat = ch_single_vcf.mix(ch_vcf_index)
                .map { meta, vcf, tbi
                -> [ meta.findAll { it.key != 'name' }, [ vcf, tbi ] ] }

    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    if (params.rename) {
        //
        // Rename samples in vcf with the label
        //
        BCFTOOLS_REHEADER(
            ch_vcf_concat.map{ it -> [ it[0], it[1][0], [], [] ] },
            [[],[]]
        )

        ch_vcf_index_rh = BCFTOOLS_REHEADER.out.vcf
                .join(BCFTOOLS_REHEADER.out.index)
                .map { meta, vcf, tbi -> [ meta, [ vcf, tbi ] ] }

        ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)
    } else {
        ch_vcf_index_rh = ch_vcf_concat
    }


    //
    // Merge multiple VCFs per sample (patient) with BCFTOOLS_MERGE
    //

    // Bring all vcfs from one sample into a channel
    // Branch based on the number of VCFs per sample
    (ch_single_id, ch_multiple_id) = ch_vcf_index_rh
        .map { meta, files ->
            // Assuming files is a list of all VCF and TBI files
            def vcfs = files.findAll { it.name.endsWith('.vcf.gz') }
            def tbis = files.findAll { it.name.endsWith('.vcf.gz.tbi') }
            [meta.id, meta, vcfs, tbis]
        }
        .groupTuple(by: 0)
        .map { id, metas, vcfs, tbis ->
            def meta = metas[0]  // Take the first meta, they should all be the same for a given ID
            def vcf_count = vcfs.flatten().size()
            meta.single_id = (vcf_count == 1)
            [meta, vcfs.flatten(), tbis.flatten()]
        }.branch {
            single: it[0].single_id
            multiple: !it[0].single_id
        }

    // Run BCFTOOLS_MERGE only on samples with multiple VCFs
    BCFTOOLS_MERGE(
        ch_multiple_id,
        [[],[]], // fasta reference only needed for gvcf
        [[],[]], // fasta.fai reference only needed for gvcf
        [[],[]] // bed
    )

    // Merge the results back into a single channel
    ch_merged_vcfs = ch_single_id.mix(BCFTOOLS_MERGE.out.vcf)

    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    //
    // Convert VCFs to Count Matrices
    //
    VCF2MAT(
        ch_merged_vcfs.map{ it -> [it[0], it[1]] },
    )

    ch_versions = ch_versions.mix(VCF2MAT.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'vcftomat_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    csv            = VCF2MAT.out.csv             // channel: *.csv
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
