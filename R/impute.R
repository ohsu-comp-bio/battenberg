# Phasing Dispatcher for Battenberg
# This file handles the high-level orchestration of haplotyping/phasing.

#' @param phasing_results_dir Directory containing the phasing output files
#' @author sd11, maxime.tarabichi, jdemeul
#' @export
run_haplotyping <- function(
  chrom, tumourname, normalname,
  ismale, problemloci,
  phasing_results_dir, min_normal_depth, chrom_names,
  reference_info_file = NA,
  externalhaplotypeprefix = NA,
  beagle_input_dir = NA,
  allele_frequencies_dir = NA,
  chrom_coord_file = NA,
  beaglejar = NA,
  beagleref_dir = NA,
  phasing_engine = "impute2",
  threads_per_chromosome = 1
) {
  # determine if we are using beagle based on engine flag or provided jar
  usebeagle <- (phasing_engine == "beagle") || (!is.na(beaglejar) && file.exists(beaglejar))

  # 1. DISCOVER OR GENERATE HAPLOTYPES
  local_haplo <- paste0(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt")

  if (file.exists(local_haplo)) {
    haplotype_file <- local_haplo
    log_info("Using existing local haplotype file: {haplotype_file}")
  } else if (usebeagle) {
    # BEAGLE FLOW
    if (!is.na(beaglejar) && file.exists(beaglejar)) {
      # Running Beagle Internal
      if (is.na(beagleref_dir) || !dir.exists(beagleref_dir)) {
        log_failure("Running Beagle internally requires a reference directory: beagleref_dir")
      }
      log_info("Running Beagle Phasing for Chromosome {chrom}")
      beagle_in <- paste0("beagle_in_chr", chrom, ".vcf")

      # Generate input for Beagle
      find_ac_file <- function(dir, sample, chrom) {
        options <- c(
          file.path(dir, paste0(sample, "_alleleFrequencies_chr", chrom, ".txt")),
          file.path(dir, paste0(sample, "_alleleFrequencies_chr", gsub("chr", "", as.character(chrom), ignore.case = TRUE), ".txt")),
          file.path(dir, paste0(sample, "_alleleFrequencies_", gsub("chr", "", as.character(chrom), ignore.case = TRUE), ".txt"))
        )
        for (f in options) {
          if (file.exists(f)) {
            return(f)
          }
        }
        return(NULL)
      }
      t_file <- find_ac_file(allele_frequencies_dir, tumourname, chrom)
      n_file <- find_ac_file(allele_frequencies_dir, normalname, chrom)
      if (is.null(t_file) || is.null(n_file)) log_failure("Could not find allele counts for phasing.")

      generate_beagle_input_from_counts(
        chrom = chrom, tumour_allele_counts_file = t_file, normal_allele_counts_file = n_file,
        output_file = beagle_in, reference_info_file = reference_info_file,
        is_male = ismale, problem_loci_file = problemloci, beagleref_dir = beagleref_dir
      )

      out_prefix <- paste0(tumourname, "_beagle_output_chr", chrom)
      vcf_out <- run_beagle_internal(
        chrom, tumourname, beagle_in, out_prefix, beaglejar, beagleref_dir,
        threads_per_chromosome = threads_per_chromosome
      )
      if (file.exists(beagle_in)) file.remove(beagle_in)

      convert_beagle_to_impute(vcf_out, local_haplo)
      haplotype_file <- local_haplo
    } else {
      # Beagle Discovery Flow (Using pre-calculated Beagle)
      beagle_search_dir <- if (!is.na(beagle_input_dir)) beagle_input_dir else phasing_results_dir
      beagle_vcf_p <- file.path(beagle_search_dir, paste0(tumourname, "_beagle5_output_chr", chrom, "_P.vcf.gz"))
      beagle_vcf_q <- file.path(beagle_search_dir, paste0(tumourname, "_beagle5_output_chr", chrom, "_Q.vcf.gz"))

      if (file.exists(beagle_vcf_p) || file.exists(beagle_vcf_q)) {
        writebeagle_as_impute_arms(vcfP = beagle_vcf_p, vcfQ = beagle_vcf_q, outfile = local_haplo)
      } else {
        # Single pattern discovery
        patterns <- c(
          paste0(tumourname, "_beagle5_output_chr", chrom, ".txt.vcf.gz"),
          paste0(tumourname, "_beagle_output_chr", chrom, ".vcf.gz")
        )
        found_vcf <- NA
        for (p in patterns) {
          tmp <- file.path(beagle_search_dir, p)
          if (file.exists(tmp)) {
            found_vcf <- tmp
            break
          }
        }
        if (is.na(found_vcf)) log_failure("Could not find pre-calculated Beagle VCF for {tumourname} chr {chrom} in {beagle_search_dir}")
        convert_beagle_to_impute(found_vcf, local_haplo)
      }
      haplotype_file <- local_haplo
    }
  } else {
    # IMPUTE2 / DIRECT DISCOVERY FLOW
    if (!is.na(phasing_results_dir)) {
      haplotype_file <- file.path(phasing_results_dir, local_haplo)
    } else {
      haplotype_file <- local_haplo
    }
    if (!file.exists(haplotype_file)) log_failure("No haplotype file found for {tumourname} chr {chrom} and no pre-phased results provided.")
  }

  # 2. TRANSFORM HAPLOTYPES INTO BAFs
  # Discovery of Allele Frequency Data
  find_ac_file <- function(dir, sample, chrom) {
    options <- c(
      file.path(dir, paste0(sample, "_alleleFrequencies_chr", chrom, ".txt")),
      file.path(dir, paste0(sample, "_alleleFrequencies_chr", gsub("chr", "", as.character(chrom), ignore.case = TRUE), ".txt")),
      file.path(dir, paste0(sample, "_alleleFrequencies_", gsub("chr", "", as.character(chrom), ignore.case = TRUE), ".txt"))
    )
    for (f in options) {
      if (file.exists(f)) {
        return(f)
      }
    }
    return(NULL)
  }
  allelefrequenciesfile <- find_ac_file(allele_frequencies_dir, tumourname, chrom)

  if (!is.null(allelefrequenciesfile) && file.exists(allelefrequenciesfile)) {
    # WGS FLOW
    if (!is.na(externalhaplotypeprefix) && file.exists(paste0(externalhaplotypeprefix, chrom, ".vcf"))) {
      # Incorporate external hapblocks
      ext_baf <- paste0(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt")
      GetChromosomeBAFs(chrom, allelefrequenciesfile, haplotype_file, tumourname, ext_baf, chrom_names, min_normal_depth)
      plot_haplotype_data(ext_baf, paste0(tumourname, "_chr", chrom, "_heterozygousData_noExt.png"), tumourname, chrom)
      input_known_haplotypes(chrom, chrom_names, haplotype_file, paste0(externalhaplotypeprefix, chrom, ".vcf"))
    }
    GetChromosomeBAFs(
      chrom, allelefrequenciesfile, haplotype_file, tumourname,
      paste0(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
      chrom_names, min_normal_depth
    )
  } else {
    # SNP6 FLOW
    GetChromosomeBAFs_SNP6(
      chrom, paste0(tumourname, "_impute_input_chr", chrom, "_withAlleleFreq.csv"),
      haplotype_file, tumourname,
      paste0(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"), chrom_names
    )
  }

  # Final Plot
  plot_haplotype_data(
    paste0(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
    paste0(tumourname, "_chr", chrom, "_heterozygousData.png"), tumourname, chrom
  )
}

#' @export
run_haplotyping_germline <- function(
  chrom, germlinename, normalname, ismale, problemloci,
  phasing_results_dir, min_normal_depth, chrom_names,
  reference_info_file = NA,
  externalhaplotypeprefix = NA,
  beagle_input_dir = NA,
  allele_frequencies_dir = NA,
  chrom_coord_file = NA,
  beaglejar = NA,
  beagleref_dir = NA,
  phasing_engine = "impute2",
  threads_per_chromosome = 8
) {
  usebeagle <- (phasing_engine == "beagle") || (!is.na(beaglejar) && file.exists(beaglejar))
  local_haplo <- paste0(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt")

  if (file.exists(local_haplo)) {
    haplotype_file <- local_haplo
  } else if (usebeagle) {
    if (!is.na(beaglejar) && file.exists(beaglejar)) {
      if (is.na(beagleref_dir) || !dir.exists(beagleref_dir)) log_failure("Internal Beagle requires beagleref_dir")
      log_info("Running internal Beagle for Germline chr {chrom}")
      beagle_in <- paste0("beagle_in_chr", chrom, ".vcf")

      find_ac_file <- function(dir, sample, chrom) {
        opts <- c(
          file.path(dir, paste0(sample, "_alleleFrequencies_chr", chrom, ".txt")),
          file.path(dir, paste0(sample, "_alleleFrequencies_chr", gsub("chr", "", as.character(chrom), ignore.case = TRUE), ".txt"))
        )
        for (f in opts) {
          if (file.exists(f)) {
            return(f)
          }
        }
        return(NULL)
      }
      ac_file <- find_ac_file(allele_frequencies_dir, germlinename, chrom)
      if (is.null(ac_file)) log_failure("No allele frequencies for Germline chr {chrom}")

      generate_beagle_input_from_counts(chrom, ac_file, ac_file, beagle_in, reference_info_file, ismale, problemloci, beagleref_dir = beagleref_dir)
      vcf_out <- run_beagle_internal(
        chrom, germlinename, beagle_in, paste0(germlinename, "_beagle_output_chr", chrom),
        beaglejar, beagleref_dir,
        threads_per_chromosome = threads_per_chromosome
      )
      if (file.exists(beagle_in)) file.remove(beagle_in)
      convert_beagle_to_impute(vcf_out, local_haplo)
      haplotype_file <- local_haplo
    } else {
      beagle_search_dir <- if (!is.na(beagle_input_dir)) beagle_input_dir else phasing_results_dir
      patterns <- c(
        paste0(germlinename, "_beagle_output_chr", chrom, ".vcf.gz"),
        paste0(germlinename, "_beagle5_output_chr", chrom, ".txt.vcf.gz")
      )
      found_vcf <- NA
      for (p in patterns) {
        tmp <- file.path(beagle_search_dir, p)
        if (file.exists(tmp)) {
          found_vcf <- tmp
          break
        }
      }
      if (is.na(found_vcf)) log_failure("Could not find pre-phased Beagle VCF for Germline")
      convert_beagle_to_impute(found_vcf, local_haplo)
      haplotype_file <- local_haplo
    }
  } else {
    haplotype_file <- if (!is.na(phasing_results_dir)) file.path(phasing_results_dir, local_haplo) else local_haplo
    if (!file.exists(haplotype_file)) log_failure("Expected haplotype file for germline missing: {haplotype_file}")
  }

  # 2. TRANSFORM HAPLOTYPES INTO BAFs (Restore missing logic for germline)
  find_ac_file <- function(dir, sample, chrom) {
    opts <- c(
      file.path(dir, paste0(sample, "_alleleFrequencies_chr", chrom, ".txt")),
      file.path(dir, paste0(sample, "_alleleFrequencies_chr", gsub("chr", "", as.character(chrom), ignore.case = TRUE), ".txt"))
    )
    for (f in opts) {
      if (file.exists(f)) {
        return(f)
      }
    }
    return(NULL)
  }
  allelefrequenciesfile <- find_ac_file(allele_frequencies_dir, germlinename, chrom)

  if (!is.null(allelefrequenciesfile) && file.exists(allelefrequenciesfile)) {
    if (!is.na(externalhaplotypeprefix) && file.exists(paste0(externalhaplotypeprefix, chrom, ".vcf"))) {
      ext_baf <- paste0(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt")
      GetChromosomeBAFs(chrom, allelefrequenciesfile, haplotype_file, germlinename, ext_baf, chrom_names, min_normal_depth)
      plot_haplotype_data(ext_baf, paste0(germlinename, "_chr", chrom, "_heterozygousData_noExt.png"), germlinename, chrom)
      input_known_haplotypes(chrom, chrom_names, haplotype_file, paste0(externalhaplotypeprefix, chrom, ".vcf"))
    }
    GetChromosomeBAFs(
      chrom, allelefrequenciesfile, haplotype_file, germlinename,
      paste0(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
      chrom_names, min_normal_depth
    )
  } else {
    log_failure("Germline calling requires WGS allele counts.")
  }

  plot_haplotype_data(
    paste0(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
    paste0(germlinename, "_chr", chrom, "_heterozygousData.png"), germlinename, chrom
  )
}
