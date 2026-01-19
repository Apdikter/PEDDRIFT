peddrift <- function(pedigree, allele_freq = NULL, nrep = 1000, rseed = 0, prune = TRUE, allall = FALSE) {

  # Load required packages
  require(dplyr)
  require(tidyr)
  require(ggplot2)

  # Set random seed
  if (rseed != 0) set.seed(rseed)

  # Validate input
  required_cols <- c("tag", "sire", "dam", "line")
  if (!all(required_cols %in% names(pedigree))) {
    stop("Pedigree must contain columns: tag, sire, dam, line")
  }

  has_locus <- "locus" %in% names(pedigree)

  # Check for locus data
  if (!has_locus) {
    stop("Must provide 'locus' column with genotype data")
  }

  # Detect genotype format and convert if needed
  if (has_locus) {
    # Get non-missing genotypes
    non_missing_geno <- pedigree$locus[!is.na(pedigree$locus)]

    if (length(non_missing_geno) == 0) {
      stop("No genotypes found in locus column")
    }


    # Check first few genotypes to detect format
    sample_geno <- head(non_missing_geno, 100)

    # Check if numeric (0/1/2 format)
    is_numeric <- all(sample_geno %in% c(0, 1, 2, "0", "1", "2", NA))

    # Check if character with 2 letters (AB format)
    is_letter <- is.character(sample_geno) && all(nchar(sample_geno) == 2)

    if (is_numeric) {
      cat("Detected 0/1/2 genotype format. Converting to AB format...\n")
      pedigree <- pedigree %>%
        mutate(locus = case_when(
          is.na(locus) ~ NA_character_,
          locus %in% c(0, "0") ~ "AA",
          locus %in% c(1, "1") ~ "AB",
          locus %in% c(2, "2") ~ "BB",
          TRUE ~ NA_character_
        ))
    } else if (is_letter) {
      cat("Detected two-letter genotype format (e.g., 'AA', 'AB', 'BB')\n")

    } else {
      stop("Genotype format not recognized. Must be either:\n",
           "  - Numeric: 0, 1, 2 (for homozygous ref, heterozygous, homozygous alt)\n",
           "  - Character: Two letters like 'AA', 'AB', 'BB'")
    }

    # Auto-assign test = 1 for animals with genotypes
    if (!("test" %in% names(pedigree))) {
      cat("Creating 'test' column: assigning test=1 to all animals with genotypes\n")
      pedigree <- pedigree %>%
        mutate(test = ifelse(!is.na(locus), 1, 0))
    } else {
      # Update test column for animals with genotypes
      cat("Updating 'test' column: setting test=1 for all animals with genotypes\n")
      pedigree <- pedigree %>%
        mutate(test = ifelse(!is.na(locus), 1, test))
    }

    cat(sprintf("Number of genotyped animals: %d\n\n", sum(pedigree$test == 1)))
  }

  # Extract locus information if available
  if (has_locus) {
    genotyped <- pedigree %>% filter(test == 1)

    # Create allele dataset
    alleles <- genotyped %>%
      mutate(
        allele1 = substr(locus, 1, 1),
        allele2 = substr(locus, 2, 2)
      ) %>%
      select(line, allele1, allele2) %>%
      pivot_longer(cols = c(allele1, allele2),
                   names_to = "allele_pos",
                   values_to = "allele") %>%
      filter(!is.na(allele))

    # Calculate contingency table chi-square
    cont_table <- table(alleles$allele, alleles$line)
    chi_test <- chisq.test(cont_table)
    obs_chisq <- chi_test$statistic
    obs_pchi <- chi_test$p.value

    cat("\n=== PEDDRIFT: Simulation of Genetic Drift ===\n")
    cat("Locus analyzed\n\n")
  } else {
    obs_chisq <- 0
    obs_pchi <- NA
    obs_diff <- 0
    obs_t <- 0
    cat("\n=== PEDDRIFT: Simulation of Genetic Drift ===\n")
    cat("No locus specified - test statistics set to zero\n\n")
  }

  # Count number of selection lines
  nsel <- pedigree %>%
    filter(test == 1) %>%
    pull(line) %>%
    unique() %>%
    length()

  # Check number of selection lines is equal to 2
  if(nsel != 2) stop("Input data must have 2 selection lines.")

  # Set up allele frequencies
  if (is.null(allele_freq)) {
    if (!has_locus) {
      stop("If no locus specified, must provide allele_freq")
    }

    # Calculate mean of line frequencies (equal weight per line) for initial line frequencies if allele_freq is not specified.
    allele_freq <- alleles %>%
      group_by(line, allele) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(line) %>%
      mutate(percent = 100 * n / sum(n)) %>%
      group_by(allele) %>%
      summarise(freq = sum(percent) / (nsel * 100), .groups = "drop")
  }

  # Print initial (founder) allele frequencies
  cat("Initial allele frequencies:\n")
  print(allele_freq, row.names = FALSE)
  cat("Sum:", sum(allele_freq$freq), "\n\n")

  nalls <- nrow(allele_freq)

  # For 2x2 case, calculate additional statistics
  if (nsel == 2 && nalls == 2 && has_locus) {
    freq_table <- alleles %>%
      group_by(line, allele) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = allele, values_from = count, values_fill = 0)

    line_totals <- rowSums(freq_table[, -1])
    p1 <- freq_table[1, 2] / line_totals[1]
    p2 <- freq_table[2, 2] / line_totals[2]

    obs_diff <- abs(as.numeric(p1 - p2))

    # Calculate t-statistic
    base_freq <- allele_freq$freq[1]
    sed <- sqrt(base_freq * (1 - base_freq) * sum(1 / line_totals))
    obs_t <- obs_diff / sed
  } else {
    obs_diff <- 0
    obs_t <- 0
  }

  # Print observed statistics
  cat("Size and test of original data:\n")
  cat(sprintf("Number of lines:             %d\n", nsel))
  cat(sprintf("Number of alleles:           %d\n", nalls))
  if (nsel == 2 && nalls == 2) {
    cat(sprintf("Allele frequency difference: %.4f\n", obs_diff))
    cat(sprintf("t test for difference:       %.4f\n", obs_t))
  }
  cat(sprintf("Contingency chi-square:      %.4f\n\n", obs_chisq))

  # Run simulations if requested
  if (nrep > 0) {
    cat(sprintf("Running %d simulation replicates...\n", nrep))

    # Prune pedigree if requested
    if (prune) {
      nindiv_before <- nrow(pedigree)
      need <- pedigree$test

      # Work backwards to mark all ancestors of genotyped animals (e.g., test = 1)
      for (i in nrow(pedigree):1) {
        if (need[i] == 1) {
          if (!is.na(pedigree$sire[i])) {
            sire_pos <- which(pedigree$tag == pedigree$sire[i])
            need[sire_pos] <- 1
          }
          if (!is.na(pedigree$dam[i])) {
            dam_pos <- which(pedigree$tag == pedigree$dam[i])
            need[dam_pos] <- 1
          }
        }
      }

      ped <- pedigree[need == 1, ]
      cat(sprintf("Number of animals before pruning: %d\n", nindiv_before))
      cat(sprintf("Number of animals after pruning:  %d\n", nrow(ped)))
    }








    # Create lookup tables for parents
    tag_to_idx <- setNames(1:nrow(ped), ped$tag)

   sire_idx <- sapply(ped$sire, function(s) {
      if (is.na(s)) return(NA)
      tag_to_idx[as.character(s)]
    })

    dam_idx <- sapply(ped$dam, function(d) {
      if (is.na(d)) return(NA)
      tag_to_idx[as.character(d)]
    })








    # Prepare allele frequency cumulative distribution
    cumfreq <- cumsum(allele_freq$freq)

    # Initialize storage for results
    sim_chisq <- numeric(nrep)
    sim_diff <- numeric(nrep)
    sim_t <- numeric(nrep)

    # Run simulations
    for (rep in 1:nrep) {
      # Initialize genotypes
      geno <- matrix(NA, nrow = nrow(ped), ncol = 2)

      # Simulate genotypes through pedigree
      for (i in 1:nrow(ped)) {
        # Allele 1 (from sire)
        if (is.na(sire_idx[i])) {
          # Sample from base population
          geno[i, 1] <- sum(cumfreq < runif(1)) + 1
        } else {
          # Inherit from sire
          geno[i, 1] <- geno[sire_idx[i], sample(1:2, 1)]
        }

        # Allele 2 (from dam)
        if (is.na(dam_idx[i])) {
          # Sample from base population
          geno[i, 2] <- sum(cumfreq < runif(1)) + 1
        } else {
          # Inherit from dam
          geno[i, 2] <- geno[dam_idx[i], sample(1:2, 1)]
        }
      }

      # Count alleles in genotyped animals by line
      sampact <- matrix(0, nrow = nsel, ncol = nalls)

      genotyped_idx <- which(ped$test == 1)
      for (i in genotyped_idx) {
        line_num <- ped$line[i]
        sampact[line_num, geno[i, 1]] <- sampact[line_num, geno[i, 1]] + 1
        sampact[line_num, geno[i, 2]] <- sampact[line_num, geno[i, 2]] + 1
      }

      # Check if all alleles present
      allemarg <- colSums(sampact)
      keep_result <- if (allall) all(allemarg > 0) else TRUE

      if (keep_result) {
        # Calculate chi-square
        linemarg <- rowSums(sampact)
        ntest <- sum(sampact)
        expect <- outer(linemarg, allemarg) / ntest
        sim_chisq[rep] <- sum((sampact - expect)^2 / expect)

        # For 2x2, calculate additional statistics
        if (nsel == 2 && nalls == 2) {
          p1_sim <- sampact[1, 1] / linemarg[1]
          p2_sim <- sampact[2, 1] / linemarg[2]
          diff_sim <- p1_sim - p2_sim

          pnull <- (sampact[1, 1] / linemarg[1] + sampact[2, 1] / linemarg[2]) / 2
          var_sim <- pnull * (1 - pnull) * sum(1 / linemarg)

          sim_diff[rep] <- diff_sim
          sim_t[rep] <- abs(diff_sim / sqrt(var_sim))
        }
      } else {
        sim_chisq[rep] <- NA
        sim_diff[rep] <- NA
        sim_t[rep] <- NA
      }

      if (rep %% 100 == 0) cat(sprintf("  Completed %d/%d replicates\n", rep, nrep))
    }

    # Remove NAs if allall=TRUE
    if (allall) {
      valid <- !is.na(sim_chisq)
      sim_chisq <- sim_chisq[valid]
      sim_diff <- sim_diff[valid]
      sim_t <- sim_t[valid]
    }

    # Calculate p-values
    pval_chisq <- mean(sim_chisq >= obs_chisq)

    cat("\n=== Results ===\n")
    cat(sprintf("P-value (contingency chi-square): %.4f\n", pval_chisq))

    if (nsel == 2 && nalls == 2) {
      pval_diff <- mean(abs(sim_diff) >= obs_diff)
      pval_t <- mean(sim_t >= obs_t)
      cat(sprintf("P-value (allele freq difference): %.4f\n", pval_diff))
      cat(sprintf("P-value (t-test):                 %.4f\n", pval_t))
    }

    # Create plots
    # Chi-square histogram with theoretical distribution
    df <- (nsel - 1) * (nalls - 1)

    p1 <- ggplot(data.frame(x2 = sim_chisq), aes(x = x2)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30,
                     fill = "skyblue", color = "black", alpha = 0.7) +
      stat_function(fun = dchisq, args = list(df = df),
                    color = "red", linewidth = 1) +
      geom_vline(xintercept = obs_chisq, color = "blue",
                 linewidth = 1, linetype = "dashed") +
      annotate("text", x = obs_chisq, y = Inf,
               label = "Observed", vjust = 1.5, hjust = -0.1, angle = 90) +
      labs(title = "Null Distribution of Chi-Square Statistic",
           subtitle = sprintf("Red line: theoretical chi-square(df=%d)", df),
           x = "Chi-square statistic", y = "Density") +
      theme_minimal()

    print(p1)

    # For 2x2, additional plots
    if (nsel == 2 && nalls == 2) {
      p2 <- ggplot(data.frame(diff = sim_diff), aes(x = diff)) +
        geom_histogram(bins = 30, fill = "lightgreen",
                       color = "black", alpha = 0.7) +
        geom_vline(xintercept = c(-obs_diff, obs_diff),
                   color = "blue", linewidth = 1, linetype = "dashed") +
        labs(title = "Null Distribution of Allele Frequency Difference",
             x = "Allele frequency difference", y = "Count") +
        theme_minimal()

      print(p2)

      p3 <- ggplot(data.frame(t = sim_t), aes(x = t)) +
        geom_histogram(bins = 30, fill = "lightyellow",
                       color = "black", alpha = 0.7) +
        geom_vline(xintercept = obs_t, color = "blue",
                   linewidth = 1, linetype = "dashed") +
        labs(title = "Null Distribution of t-statistic",
             x = "t-statistic", y = "Count") +
        theme_minimal()

      print(p3)
    }

    # Return results
    results <- list(
      observed = list(
        chisq = obs_chisq,
        diff = obs_diff,
        t = obs_t
      ),
      pvalues = list(
        chisq = pval_chisq,
        diff = if (nsel == 2 && nalls == 2) pval_diff else NA,
        t = if (nsel == 2 && nalls == 2) pval_t else NA
      ),
      simulated = data.frame(
        chisq = sim_chisq,
        diff = sim_diff,
        t = sim_t
      ),
      summary = list(
        n_lines = nsel,
        n_alleles = nalls,
        n_replicates = length(sim_chisq),
        allele_freq = allele_freq
      )
    )

    return(invisible(results))

  } else {
    cat("No simulations run (nrep = 0)\n")
    return(NULL)
  }
}
