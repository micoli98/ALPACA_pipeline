process PREPARE_FILTER_INPUT {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}", mode: "copy"
    conda params.svclone_env
    
    input: 
    tuple val(sample),
        val(patient),
        path(cnv),
        path(pp),
        path(snv_pyclone),
        path(snv_vcf)

    output:
    tuple val(sample),
        path("${sample}_snvs_for_svclone.vcf"),
        path("${sample}_ascat.csv"),
        path("${sample}_pp.tsv")

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(VariantAnnotation)
        library(GenomicRanges)
        library(tidyverse)
    })
    source("${params.conversion_script}")

    sample_name   <- '${sample}'
    pyclone_snv   <- '${snv_pyclone}'
    snv_vcf       <- '${snv_vcf}'
    cnv_tsv       <- '${cnv}'
    purity_ploidy <- '${pp}'

    ##### Helper functions #####

    # mutation_key() creates a standardized mutation identifier used internally for
    # matching the PyClone SNV table to the VCF.
    #
    # Example input: chrom = "chr1", pos = 12345  ->  "1:12345"
    #
    # This deliberately removes the "chr" prefix, so chr1 and 1 match each other.
    mutation_key <- function(chrom, pos) {
        chrom <- sub("^chr", "", trimws(as.character(chrom)), ignore.case = TRUE)
        paste(chrom, as.integer(pos), sep = ":")
    }

    # Extraction of the ALT read depth for the tumour sample from the VCF
    # genotype field AD.
    #
    # In many VCFs, AD means allelic depth. For a biallelic variant:
    #   AD = 20,5   -> means: REF reads = 20 ALT reads = 5
    #
    # The script uses this to remove variants where the tumour sample has zero ALT
    # reads, so making the extracted SNVs sample specific
    ad_alt_depth <- function(ad, sample_name) {
        # If the VCF does not contain AD at all, we cannot check ALT depth.
        if (is.null(ad)) stop("VCF FORMAT field AD is missing")

        # VariantAnnotation may store AD as a 3D array:
        #   dimension 1 = variants
        #   dimension 2 = samples
        #   dimension 3 = alleles
        # In this structure, allele 1 is REF and allele 2 is ALT after VCF expansion.
        if (length(dim(ad)) == 3) {
            return(as.integer(ad[, sample_name, 2]))
        }

        # In some VCFs, VariantAnnotation may store AD in a matrix-like structure,
        # where each cell contains the vector of allele depths for one variant/sample.
        if (is.matrix(ad) || is.data.frame(ad)) {
            x <- ad[, sample_name]

            return(vapply(x, function(z) {
            z <- suppressWarnings(as.integer(unlist(z)))
            if (length(z) < 2 || is.na(z[2])) NA_integer_ else z[2]
            }, integer(1)))
        }

        # If AD has an unexpected structure, stop clearly instead of silently doing
        # something wrong.
        stop("Do not know how to extract AD from this VCF object")
    }

    ###############################################################################
    # 1. Basic input checks
    ###############################################################################

    # Stop immediately if any required input file is missing.
    for (path in c(pyclone_snv, snv_vcf, cnv_tsv, purity_ploidy)) {
        if (!file.exists(path)) stop("Missing input file: ", path)
    }

    ###############################################################################
    # 2. Create the purity/ploidy file required by SVclone
    ###############################################################################

    # Read purity/ploidy input table.
    pp_in <- read_tsv(purity_ploidy, show_col_types = F) |>
        # If the table contains multiple samples, keep only this sample.
        # If there is no sample column, we assume the file contains only one row/sample.
        filter(sample == sample_name)

    # If no row remains, the requested sample was not found.
    if (nrow(pp_in) == 0) stop("No purity/ploidy row found for sample ", sample_name)

    # SVclone expects a small file with sample, purity, and ploidy.
    # We take the first matching row if more than one row exists.
    write_tsv(
        data.frame(
            sample = sample_name,
            purity = pp_in\$purity[1],
            ploidy = pp_in\$ploidy[1]
        ),
        file.path(paste0(sample_name, "_pp.tsv")))

    ###############################################################################
    # 3. Create the ASCAT-like CNV file required by SVclone
    ###############################################################################

    # Read CNV segment table.
    cnv <- read_tsv(cnv_tsv, show_col_types = F) |>
        # Column names expected in the CNV table.
        dplyr::select("sample", "chromosome", "start", "end", "copyNumber", "minorAlleleCopyNumber") |>
        # If the CNV table contains multiple samples, keep only this sample.
        filter(sample == sample_name)
    if (nrow(cnv) == 0) stop("No CNV rows found for sample ", sample_name)

    # Build the ASCAT-like file.
    # The output has no column names because the original Python script wrote rows
    # directly without a header.
    #
    # Columns are:
    #   1) segment id
    #   2) chromosome
    #   3) start
    #   4) end
    #   5) normal total copy number, fixed to 2
    #   6) normal minor copy number, fixed to 1
    #   7) tumour total copy number
    #   8) tumour minor copy number
    ascat <- data.frame(
        id = seq_len(nrow(cnv)),
        chromosome = cnv\$chromosome,
        start = cnv\$start,
        end = cnv\$end,
        normal_total = 2,
        normal_minor = 1,
        tumour_total = cnv\$copyNumber,
        tumour_minor = cnv\$minorAlleleCopyNumber,
        check.names = FALSE
        )

    write_csv(ascat, file.path(paste0(sample_name, "_ascat.csv")),, col_names = F)

    ###############################################################################
    # 4. Read the PyClone SNV table and create the SNV whitelist
    ###############################################################################

    # Read the PyClone SNV table.
    snv_table <- read_tsv(pyclone_snv, show_col_types = F)

    # The only column we need from this table is mutation_id.
    # It is used as a whitelist: only VCF variants matching these mutation IDs are
    # kept for the SVclone input VCF.
    if (!("mutation_id" %in% colnames(snv_table))) {
        stop("pyclone_snv.tsv must contain mutation_id column")
    }

    # Parse mutation_id values.
    matches <- regexec(
        "^(chr)?([^:]+):([0-9]+)",
        snv_table\$mutation_id,
        ignore.case = TRUE
    )

    hits <- regmatches(snv_table\$mutation_id, matches)

    # Convert all valid mutation IDs to the same internal key format used for VCF
    # variants: chrom:pos:REF:ALT, with chr removed and alleles upper-cased.
    selected <- unique(vapply(
        hits[lengths(hits) > 0],
        function(x) mutation_key(x[3], x[4]),
        character(1)
    ))

    ###############################################################################
    # 5. Read the VCF with VariantAnnotation
    ###############################################################################

    vcf <- readVcf(snv_vcf)

    # Some source VCFs (e.g. the v4.10 CADD annotation batch) have INFO header
    # lines with unquoted, multi-word Description values (spec-noncompliant).
    # VariantAnnotation misparses these on read (Type absorbs the free text,
    # Description becomes NA) and re-emits a malformed, duplicate-key INFO
    # line on write, which later breaks SVclone's stricter VCF parser. None of
    # these fields (PolyPhenCat/PolyPhenVal/SIFTcat/SIFTval, ...) are used
    # below, so just drop any INFO field whose header failed to parse.
    info_header <- info(header(vcf))
    malformed <- which(is.na(info_header\$Description))
    if (length(malformed) > 0) {
        bad_ids <- rownames(info_header)[malformed]
        info(header(vcf)) <- info_header[-malformed, , drop = FALSE]
        info(vcf) <- info(vcf)[, !(colnames(info(vcf)) %in% bad_ids), drop = FALSE]
    }

    # Name conversion
    conv_table <- read.table("${params.conversion_table}", sep="\t", header = T)
    updated_cols <- bind_rows(lapply(colnames(vcf), function(id) sample_name_converter(id, conv_table)))\$new_name
    colnames(vcf) <- updated_cols
    print(updated_cols)
    
    # Check that the VCF contains the tumour sample.
    # The output VCF will contain only these two samples.
    if (!(sample_name %in% colnames(vcf))) {
        stop("Sample '", sample_name, "' not found in ", snv_vcf)
    }

    normal_sample <- colnames(vcf)[grepl("BDNA|BBC|ME", colnames(vcf))]
    if ((length(normal_sample) != 1)) {
        stop("Normal sample not found or multiple matches in ", snv_vcf)
    }

    # Keep only the normal and tumour columns.
    # This removes unrelated samples from the output VCF.
    vcf <- vcf[, c(normal_sample, sample_name)]

    ###############################################################################
    # 6. Split multiallelic VCF rows into one ALT allele per row
    ###############################################################################

    # A VCF row can contain several alternate alleles, for example:
    #   REF=A ALT=G,T
    #
    # SVclone input is easier to handle as biallelic rows, for example:
    #   REF=A ALT=G
    #   REF=A ALT=T
    #
    # Before expanding the VCF, we store which ALT allele each expanded row will
    # represent.
    #
    # Example:
    #   original row 1 has ALT = C       -> alt index 1
    #   original row 2 has ALT = G,T     -> alt indexes 1 and 2
    #   original row 3 has ALT = A,C,G   -> alt indexes 1, 2, and 3
    #
    # sequence(c(1, 2, 3)) gives:
    #   1, 1, 2, 1, 2, 3
    #
    # These indexes are needed later to recode GT correctly.
    alt_index <- sequence(S4Vectors::elementNROWS(alt(vcf)))

    # expand() turns multiallelic VCF rows into one row per ALT allele.
    vcf <- VariantAnnotation::expand(vcf, row.names = FALSE)

    ###############################################################################
    # 7. Create mutation keys for VCF rows and match to PyClone whitelist
    ###############################################################################

    # rowRanges(vcf) contains genomic coordinates for each VCF row.
    rr <- rowRanges(vcf)

    # Create a mutation key for every expanded VCF row.
    # This creates the same format as the PyClone whitelist:
    #   chrom:pos:REF:ALT
    vcf_keys <- mutation_key(as.character(seqnames(rr)), start(rr))

    # selected_in_vcf is TRUE for VCF rows that appear in pyclone_snv.tsv.
    selected_in_vcf <- vcf_keys %in% selected
    print(paste0("Matched ", length(selected %in% vcf_keys), " out of ", length(selected), 
    " SNVs from pyclone to VCF"))

    # The FILTER column of the VCF indicates whether a variant passed the caller's
    # filters. This script treats PASS and . as acceptable.
    pass_record <- filt(vcf) %in% c("PASS", ".")

    keep <- selected_in_vcf & pass_record

    ###############################################################################
    # 8. Remove variants with zero tumour ALT depth
    ###############################################################################

    # Extract tumour ALT depth from the AD field.
    # After expand(), each row has only one ALT allele, so AD should correspond to:
    #   REF depth, ALT depth
    alt_depth <- ad_alt_depth(geno(vcf)\$AD, sample_name)

    # Keep only variants with positive tumour ALT depth.
    keep <- keep & !is.na(alt_depth) & alt_depth > 0

    # Apply the final VCF row filter.
    vcf <- vcf[keep, ]

    # Apply the same filter to alt_index so it remains aligned with vcf rows.
    alt_index <- alt_index[keep]

    ###############################################################################
    # 9. Recode GT after multiallelic expansion
    ###############################################################################

    # This block fixes genotype values after splitting multiallelic variants.
    #
    # Example before expansion:
    #   REF=A ALT=G,T GT=0/2
    #
    # If we keep the T allele, the expanded row becomes:
    #   REF=A ALT=T
    #
    # In a biallelic row, the kept ALT allele should be coded as 1, not 2.
    # Therefore GT=0/2 should become GT=0/1.
    #
    # Any other ALT allele that is not the current kept ALT is recoded to 0.
    # For example, if the current row is ALT=T and the original genotype was 1/2:
    #   original ALT 1 = G
    #   original ALT 2 = T
    #   new GT = 0/1
    # because G is no longer represented in this row.
    if (nrow(vcf) > 0 && !is.null(geno(vcf)\$GT)) {
    recode_gt <- function(gt, alt_i) {
        # Leave missing or empty genotypes unchanged.
        if (is.na(gt) || gt == "") return(gt)

        # Preserve phasing separator if GT uses | instead of /.
        sep <- if (grepl("\\\\|", gt)) "|" else "/"

        # Split genotype into allele numbers.
        alleles <- strsplit(gt, "[/|]")[[1]]

        # Allele coding in the original VCF:
        #   0 = REF
        #   1 = first ALT
        #   2 = second ALT
        #   . = missing
        #
        # In the expanded biallelic VCF row:
        #   0 = REF
        #   1 = the ALT allele represented by this row
        #   . = missing

        # Any ALT allele that is not the current ALT becomes REF-like 0 in this
        # biallelic representation.
        alleles[alleles != "." & alleles != "0" & alleles != as.character(alt_i)] <- "0"

        # The current ALT allele becomes 1.
        alleles[alleles == as.character(alt_i)] <- "1"

        # Join the genotype back together using the original separator.
        paste(alleles, collapse = sep)
    }

    # Recode GT for the normal sample.
    geno(vcf)\$GT[, normal_sample] <- mapply(
        recode_gt,
        geno(vcf)\$GT[, normal_sample],
        alt_index,
        USE.NAMES = FALSE
    )

    # Recode GT for the tumour sample.
    geno(vcf)\$GT[, sample_name] <- mapply(
        recode_gt,
        geno(vcf)\$GT[, sample_name],
        alt_index,
        USE.NAMES = FALSE
    )
    }

    ###############################################################################
    # 10. Write the filtered VCF
    ###############################################################################

    # Rename the normal column to exactly "normal" in the output.
    # Keep the tumour column as the sample name.
    colnames(vcf) <- c("normal", sample_name)

    out_vcf <- file.path(paste0(sample_name, "_snvs_for_svclone.vcf"))
    writeVcf(vcf, out_vcf, index = FALSE)
    """
}