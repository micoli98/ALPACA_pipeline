process CONFIG {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", pattern: "svclone_config.ini", mode: "copy"

    input: 
    tuple val(sample),
        val(patient),
        val(meanCov), 
        val(readLen), 
        val(insertSize),
        val(insertStd)

    output:
    tuple val(sample), path("svclone_config.ini")

    script:
    def chroms = (1..22).collect{ "chr${it}" }.join(',') + ',chrX,chrY'
    """
    wget -O svclone_config_template.ini https://raw.githubusercontent.com/mcmero/SVclone/master/svclone_config.ini
    
    cat svclone_config_template.ini | \\
        sed 's/read_len: -1/read_len: ${readLen}/' | \\
        sed 's/insert_mean: -1/insert_mean: ${insertSize}/' | \\
        sed 's/insert_std: -1/insert_std: ${insertStd}/' | \\
        sed 's/mean_cov: 50/mean_cov: ${meanCov}/' | \\
        sed 's/max_cn: 10/max_cn: 100/' | \\
        sed 's/min_dep: 8/min_dep: 4/' | \\
        sed 's/filter_chroms: False/filter_chroms: True/' | \\
        sed 's/support_adjust_factor: 0/support_adjust_factor: 0.2/' | \\
        sed 's/repeat: 5/repeat: 10/' | \\
        sed 's/clus_limit: 6/clus_limit: 10/' | \\
        sed 's/male: True/male: False/' | \\
        sed 's/chroms: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y/chroms: ${chroms}/' \\
        > svclone_config.ini
    """
}