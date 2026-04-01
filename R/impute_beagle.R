#' Helper for writing Beagle VCFs
#' @export
writevcf_beagle <- function(vcf, filepath, vcfversion = "4.2", genomereference = "GRCh38") {
    header <- paste0(
        "##fileformat=VCFv", vcfversion, "\n",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "##reference=", genomereference, "\n"
    )
    cat(header, file = filepath)
    data.table::fwrite(vcf, file = filepath, sep = "\t", append = TRUE, col.names = TRUE, quote = FALSE)
}

#' Convert intermediate Battenberg format to Beagle VCF
#' @export
convert_impute_input_to_beagle_vcf <- function(impute_input_data, chrom) {
    chr_vcf <- if (chrom == "23") "X" else as.character(chrom)
    coln <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMP001")

    vcf <- data.frame(
        CHROM = rep(chr_vcf, nrow(impute_input_data)),
        POS = impute_input_data$X3,
        ID = rep(".", nrow(impute_input_data)),
        REF = impute_input_data$X4,
        ALT = impute_input_data$X5,
        QUAL = rep(".", nrow(impute_input_data)),
        FILTER = rep("PASS", nrow(impute_input_data)),
        INFO = rep(".", nrow(impute_input_data)),
        FORMAT = rep("GT", nrow(impute_input_data)),
        GT = paste(impute_input_data$X6, impute_input_data$X7, impute_input_data$X8, sep = "-"),
        stringsAsFactors = FALSE
    )
    vcf$GT[vcf$GT == "1-0-0"] <- "0/0"
    vcf$GT[vcf$GT == "0-1-0"] <- "0/1"
    vcf$GT[vcf$GT == "0-0-1"] <- "1/1"
    vcf <- vcf[vcf$GT != "0-0-0", ]
    colnames(vcf) <- coln
    return(vcf)
}

#' Generate Beagle input directly from allele counts
#' @export
generate_beagle_input_from_counts <- function(chrom, tumour_allele_counts_file, normal_allele_counts_file,
                                              output_file, reference_info_file = NA, is_male = NA,
                                              problem_loci_file = NA, heterozygous_filter = 0.1,
                                              beagleref_dir = NA) {
    # Try to find a reference legend
    known_SNPs <- NULL
    if (!is.na(reference_info_file) && file.exists(reference_info_file)) {
        impute_info <- parse_imputeinfofile(reference_info_file, is_male, chrom = chrom)
        if (nrow(impute_info) > 0) {
            log_info("Reading legend from {impute_info$impute_legend}")
            known_SNPs <- vroom::vroom(unlist(impute_info$impute_legend), delim = " ", col_types = "ciccc", show_col_types = FALSE)
            data.table::setDT(known_SNPs)
        }
    }

    # If no reference info provided, attempt discovery in beagleref_dir (expecting LEGEND-style files or subsetting reference VCF)
    # Actually, if we don't have a legend, we'll try to use the matched normal loci themselves as the "legend" if no reference is specified.
    # But for a high-quality Beagle run, we really want that legend.
    if (is.null(known_SNPs)) {
        log_warning("No reference legend found for chr {chrom}. Using all loci from normal allele counts.")
        # This might be slow if the allele counts file is huge, but it's a fallback.
    }

    log_info("Reading normal allele counts from {normal_allele_counts_file}")
    snp_normal <- read_alleleFrequencies(normal_allele_counts_file)

    # Filter problem SNPs
    if (!is.na(problem_loci_file) && problem_loci_file != "NA" && file.exists(problem_loci_file)) {
        problem_snps_raw <- data.table::fread(problem_loci_file, header = TRUE, sep = "\t", data.table = FALSE)
        problem_positions <- problem_snps_raw$Pos[problem_snps_raw$Chr == chrom]
        snp_normal <- snp_normal[!(snp_normal$POS %in% problem_positions), ]
    }

    if (!is.null(known_SNPs)) {
        common_pos <- intersect(known_SNPs$position, snp_normal$POS)
        valid_known_snps <- known_SNPs[match(common_pos, known_SNPs$position), ]
        found_normal_data <- snp_normal[match(common_pos, snp_normal$POS), ]

        # Define base columns (A=3, C=4, G=5, T=6 in our table)
        bases <- c("A", "C", "G", "T")
        ref_base_idx <- match(valid_known_snps$a0, bases)
        alt_base_idx <- match(valid_known_snps$a1, bases)

        # Extract counts
        normal_counts_matrix <- as.matrix(found_normal_data[, 3:6, with = FALSE])
        ref_counts <- as.numeric(vapply(seq_along(ref_base_idx), function(i) {
            if (is.na(ref_base_idx[i])) {
                return(0)
            }
            normal_counts_matrix[i, ref_base_idx[i]]
        }, numeric(1)))
        alt_counts <- as.numeric(vapply(seq_along(alt_base_idx), function(i) {
            if (is.na(alt_base_idx[i])) {
                return(0)
            }
            normal_counts_matrix[i, alt_base_idx[i]]
        }, numeric(1)))

        total_counts <- ref_counts + alt_counts
        keep_mask <- total_counts > 0

        vcf <- data.table::data.table(
            "#CHROM" = if (chrom == "23") "X" else as.character(chrom),
            POS = valid_known_snps$position[keep_mask],
            ID = valid_known_snps$id[keep_mask],
            REF = valid_known_snps$a0[keep_mask],
            ALT = valid_known_snps$a1[keep_mask],
            QUAL = ".",
            FILTER = "PASS",
            INFO = ".",
            FORMAT = "GT"
        )

        bafs <- alt_counts[keep_mask] / total_counts[keep_mask]
        gt <- rep("0/1", length(bafs))
        gt[bafs <= heterozygous_filter] <- "0/0"
        gt[bafs >= (1.0 - heterozygous_filter)] <- "1/1"
        vcf$SAMP001 <- gt
    } else {
        # No legend: Infer REF/ALT from counts (largest count is Ref, second largest is Alt)
        # This is sub-optimal but works for pre-phasing.
        log_info("Inferring alleles from counts for chr {chrom}")
        # (Simplified logic for now: only use the top 2 bases)
        # Actually, legacy Battenberg ALWAYS requires a legend or it fails elsewhere.
        log_failure("A reference legend is currently required to generate Beagle input. Please provide a reference_info_file.")
    }

    log_info("Writing {nrow(vcf)} SNPs to {output_file}")
    writevcf_beagle(vcf, output_file)
}

#' Convert Beagle VCF to IMPUTE format
#' @export
convert_beagle_to_impute <- function(beagle_file, output_file) {
    if (!file.exists(beagle_file)) log_failure("Beagle VCF file not found: {beagle_file}")

    vcf <- tryCatch(
        {
            data.table::fread(beagle_file, skip = "#CHROM", header = TRUE)
        },
        error = function(e) {
            if (file.info(beagle_file)$size < 500) {
                return(data.table::data.table())
            }
            stop(e)
        }
    )

    if (nrow(vcf) == 0) {
        log_info("Beagle VCF is empty. Writing empty output.")
        data.table::fwrite(data.table::data.table(), file = output_file, sep = " ", col.names = FALSE)
        return(NULL)
    }

    gt_data <- vcf[[10]]
    gt_only <- data.table::tstrsplit(gt_data, ":")[[1]]
    haplo <- data.table::tstrsplit(gt_only, "[|/]")

    impute_dt <- data.table::data.table(
        V1 = "---",
        V2 = vcf[["ID"]],
        V3 = as.integer(as.numeric(vcf[["POS"]])),
        V4 = vcf[["REF"]],
        V5 = vcf[["ALT"]],
        V6 = haplo[[1]],
        V7 = haplo[[2]]
    )
    data.table::fwrite(impute_dt, file = output_file, sep = " ", col.names = FALSE, quote = FALSE)
}

#' Split a VCF into p and q arms
#' @export
split_and_writevcf_by_arm <- function(vcf, chrom, pathP, pathQ, coord_file) {
    centromere_split <- load_centromere_splits(coord_file)
    lookup_chrom <- if (chrom == "X") "23" else as.character(chrom)
    if (!(lookup_chrom %in% names(centromere_split))) {
        log_warning("Chromosome '{chrom}' not found in centromere table. Phasing as single unit.")
        writevcf_beagle(vcf, pathP)
        return(invisible(NULL))
    }
    split_point <- centromere_split[[lookup_chrom]]
    vcf_p <- vcf[as.numeric(vcf$POS) <= split_point]
    vcf_q <- vcf[as.numeric(vcf$POS) > split_point]
    if (nrow(vcf_p) > 0) writevcf_beagle(vcf_p, pathP)
    if (nrow(vcf_q) > 0) writevcf_beagle(vcf_q, pathQ)
}

#' Merge Beagle output from p and q arms back into IMPUTE format
#' @export
writebeagle_as_impute_arms <- function(vcfP = NULL, vcfQ = NULL, outfile) {
    read_vcf <- function(path) {
        if (!is.null(path) && file.exists(path)) {
            return(data.table::fread(path, skip = "#CHROM", header = TRUE))
        }
        return(NULL)
    }
    outP <- read_vcf(vcfP)
    outQ <- read_vcf(vcfQ)
    if (is.null(outP) && is.null(outQ)) log_failure("Neither p-arm nor q-arm Beagle output found.")
    combined <- data.table::rbindlist(list(outP, outQ), use.names = TRUE)
    gt_data <- combined[[10]]
    haplo <- data.table::tstrsplit(gt_data, "[|/]")
    impute_dt <- data.table::data.table(
        V1 = "---", V2 = combined$ID, V3 = combined$POS, V4 = combined$REF, V5 = combined$ALT,
        V6 = haplo[[1]], V7 = haplo[[2]]
    )
    data.table::fwrite(impute_dt, file = outfile, sep = " ", col.names = FALSE, quote = FALSE)
}

#' Run Beagle 5 internal phasing
#' @export
run_beagle_internal <- function(chrom, samplename, beagle_in, out_prefix,
                                beaglejar, beagleref_dir,
                                threads_per_chromosome = 1) {
    norm_c <- gsub("chr", "", as.character(chrom), ignore.case = TRUE)

    # Discover reference VCF
    ref_vcf <- NA
    if (!is.na(beagleref_dir) && dir.exists(beagleref_dir)) {
        ref_pats <- c(paste0("chr", chrom, ".*vcf.gz"), paste0("chr", norm_c, ".*vcf.gz"))
        for (p in ref_pats) {
            matches <- list.files(beagleref_dir, pattern = p, full.names = TRUE)
            if (length(matches) > 0) {
                ref_vcf <- matches[1]
                break
            }
        }
    }
    if (is.na(ref_vcf)) {
        log_failure("Running Beagle internal requires a reference VCF. Could not find one for chr {chrom} in {beagleref_dir}")
    }

    beagle_cmd <- sprintf(
        "java -jar %s gt=%s out=%s nthreads=%d impute=false",
        beaglejar, beagle_in, out_prefix, threads_per_chromosome
    )
    if (!is.na(ref_vcf)) beagle_cmd <- paste0(beagle_cmd, " ref=", ref_vcf)

    log_info("Executing Beagle: {beagle_cmd}")
    system(beagle_cmd)

    # Return the expected output file path
    vcf_out <- paste0(out_prefix, ".vcf.gz")
    if (!file.exists(vcf_out)) vcf_out <- paste0(out_prefix, ".vcf")
    return(vcf_out)
}
