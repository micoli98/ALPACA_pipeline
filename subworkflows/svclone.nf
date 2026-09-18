// SVclone subworkflow: annotates, counts, filters and clusters structural
// variants by cancer cell fraction (per sample).

include { GET_INSERT_STD } from "../modules/GetInsertStd/main.nf"
include { CONFIG } from "../modules/Config/main.nf"
include { FILTER_SV } from "../modules/FilterSV/main.nf"
include { ANNOTATION } from "../modules/Annotation/main.nf"
include { COUNT } from "../modules/Count/main.nf"
include { PREPARE_FILTER_INPUT } from "../modules/PrepareFilterInput/main.nf"
include { FILTER } from "../modules/Filter/main.nf"
include { CLUSTERING } from "../modules/Clustering/main.nf"

workflow SVCLONE_SUBWORKFLOW {
    take:
    samples_info_ch  // value channel: path to sample_info TSV (columns: sample, patient)
    batch_info_ch    // value channel: path to batch_info TSV (also used for path_to_trees/patient-derived snv_pyclone path and the snv_vcf column)
    coverage_info    // value channel: path to coverage_info TSV

    main:
    /////// 0: Prepare configuration ///////
    // Exclude low purity samples
    high_purity_samples = channel.fromPath(params.pp)
        .splitCsv(sep: '\t', header: true)
        .filter{ row -> row.purity.isDouble() && (row.purity as Double) > params.min_purity }
        .map{ row -> tuple(row.sample, true) }

    // Get insert size standard deviation (provisionally with bamtools)
    sample_info_insert = samples_info_ch
        .splitCsv(sep:'\t', header: true, quote: '"')
        .map{ row -> tuple(row.sample, row.patient) }
        .join(high_purity_samples, by: 0)
        .map{ sample, patient, is_high_purity ->
            tuple(
                sample,
                patient,
                "${params.bam_dir}${sample}.bam"
            )
        }

    GET_INSERT_STD(sample_info_insert)

    // Join with the actual coverage info
    coverage_sample = coverage_info
        .splitCsv(sep:'\t', header: true)
        .map{ row-> tuple(row.sample, row.meanCoverage, row.MEAN_READ_LENGTH, row.MEAN_INSERT_SIZE)}

    sample_info_coverage = GET_INSERT_STD.out
        .join(coverage_sample, by: 0)
        .map{ sample, patient, insertStd, meanCov, readLen, insertSize ->
            tuple(
                sample,
                patient,
                Math.round(meanCov as Double) as Integer,
                Math.round(readLen as Double) as Integer,
                Math.round(insertSize as Double) as Integer,
                Math.round(insertStd as Double) as Integer
            )
        }

    // Produce the configuration files for each sample
    CONFIG(sample_info_coverage)

    /////// 1: Annotation ///////
    sample_info_sv = sample_info_coverage
        .map{ sample, patient, meanCov, readLen, insertSize, insertStd ->
            tuple(
                sample,
                patient,
                "${params.purple_dir}${patient}/${sample}.purple.sv.vcf.gz"
            )
        }

    // Filter only the breakpoints from the SV VCF
    FILTER_SV(sample_info_sv)

    // Join config and filtered SVs
    sample_info_anno = FILTER_SV.out
        .join(CONFIG.out, by: 0)
        .join(sample_info_insert, by:0)
        .map{ sample, patient, sv_vcf, config, patient2, bam ->
            tuple(sample, patient, bam,
            "${params.bam_dir}${sample}.bai",
            sv_vcf, config)
        }

    ANNOTATION(sample_info_anno)

    sample_info_count = sample_info_anno
        .join(ANNOTATION.out, by: 0)
        .map{ sample, patient, bam, bai, sv_vcf, config, svin, readinfo ->
            tuple(sample, patient, bam, bai, svin, config, readinfo)}

    /////// 2: Count ///////
    COUNT(sample_info_count)

    /////// 3: Filter ///////
    // Prepare SNV input for the filtering step
    snv_sample = batch_info_ch
        .splitCsv(sep:'\t', header: true)
        .map{ row ->
            def version_dir = row.file_cf.replace(row.path_to_trees, '').tokenize('/')[0]
            tuple(row.patient, "${row.path_to_trees}${version_dir}/${version_dir}_formatted_filtered.csv", row.snv_vcf)
        }

    prepare_input = sample_info_insert
        .map{ sample, patient, bam -> tuple(patient, sample)}
        .combine(snv_sample, by: 0)
        .map{ patient, sample, snv_pyclone, snv_vcf ->
            tuple(
                sample,
                patient,
                params.seg,
                params.pp,
                snv_pyclone,
                snv_vcf
            )
        }

    PREPARE_FILTER_INPUT(prepare_input)

    // Join the prepared inputs
    filter_sample = sample_info_count
        .map{sample, patient, bam, bai, svin, config, readinfo ->
            tuple(sample, patient, config, readinfo)
        }
        .join(COUNT.out, by: 0)
        .join(PREPARE_FILTER_INPUT.out, by:0)

    FILTER(filter_sample)

    /////// 4: Cluster ///////
    cluster_sample = filter_sample
        .map{sample, patient, config, readinfo, svinfo, snv, cnv, pp ->
            tuple(sample, patient, config, readinfo)
        }
        .join(FILTER.out, by: 0)

    CLUSTERING(cluster_sample)

    // CLUSTERING only emits (sample, files) — reattach patient so downstream
    // joins/groupings by patient are possible.
    clustering_out = cluster_sample
        .map{ sample, patient, config, readinfo, svs, snvs, pp -> tuple(sample, patient) }
        .join(CLUSTERING.out, by: 0)
        .map{ sample, patient, files -> tuple(patient, sample, files) }

    emit:
    clustering_out
}
