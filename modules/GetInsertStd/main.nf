process GET_INSERT_STD {
    tag "${sample}"

    input: 
    tuple val(sample),
        val(patient),
        val(bam)

    output:
    tuple val(sample), val(patient), stdout

    script:
    """
    samtools view -f 2 -F 3844 ${bam} \
    | awk '\$9 > 0 {print \$9}' \
    | head -n 1000000 \
    | awk '
        {n++; x+=\$1; xx+=\$1*\$1}
        END {
        mean=x/n
        sd=sqrt((xx - x*x/n)/(n-1))
        print sd
        }'
    """
}