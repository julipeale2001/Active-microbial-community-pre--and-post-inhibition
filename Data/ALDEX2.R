# Calcute and plot ALDEX2
# Roey Angel
# roey.angel@bc.cas.cz
# Feb. 2020



# Functions ---------------------------------------------------------------
DropRareSpecies <-
  function(Ps_obj = Ps_obj_subset,
           prevalence = 0.1) {
    prevdf <- apply(
      X = otu_table(Ps_obj),
      MARGIN = ifelse(taxa_are_rows(Ps_obj), yes = 1, no = 2),
      FUN = function(x) {
        sum(x > 0)
      }
    )
    # Add taxonomy and total read counts to this data.frame
    prevdf <- data.frame(Prevalence = prevdf,
                         TotalAbundance = taxa_sums(Ps_obj),
                         tax_table(Ps_obj))
    
    # Define prevalence threshold as 0.X of total samples
    prevalenceThreshold <- prevalence * nsamples(Ps_obj)
    prevalenceThreshold
    
    # Execute prevalence filter, using `prune_taxa()` function
    prevdf_phylum_filt <-
      subset(prevdf,
             Phylum %in% get_taxa_unique(Ps_obj, "Phylum"))
    keepTaxa <-
      row.names(prevdf_phylum_filt)[(prevdf_phylum_filt$Prevalence >= prevalenceThreshold)]
    Ps_obj_small <- prune_taxa(keepTaxa, Ps_obj)
    sample_data(Ps_obj_small)$Lib.size <-
      rowSums(otu_table(Ps_obj_small))
    print(Ps_obj)
    print(Ps_obj_small)
    return(Ps_obj_small)
  }

calc_ALDEx2 <-
  function(physeq_obj = Ps_obj_filt_subset,
           vars2test = "Spill.Treatment",
           rare_phyla = Rare_phyla,
           sig_level = 0.05,
           LFC = 0.322,
           ...) {
    require(phyloseq)
    require(speedyseq)
    require(ALDEx2)
    physeq_obj <- filter_taxa(physeq_obj, function(x)
      sum(x) > 0, TRUE)
    # physeq_obj <- prune_taxa(sig_taxa, physeq_obj) # remove taxa not significant under the full model
    data2test <- (otu_table(physeq_obj))
    comparison <- as.character(get_variable(physeq_obj, vars2test))
    ALDEx <- aldex.clr(
      data2test,
      comparison,
      mc.samples = 128,
      denom = "iqlr",
      # iqlr for slight assymetry in composition
      verbose = TRUE,
      useMC = TRUE
    )
    ALDEx_tt <-
      aldex.ttest(ALDEx, paired.test = FALSE, verbose = TRUE) # for two conditions
    ALDEx_effect <- aldex.effect(
      ALDEx,
      include.sample.summary = TRUE,
      verbose = TRUE,
      useMC = TRUE
    ) # estimate effect sizes
    ALDEx2plot <-
    prep_AlDEx_data(ALDEx_tt,
                    ALDEx_effect,
                    physeq_obj,
                    sig_level,
                    LFC,
                    Taxa_rank,
                    rare_phyla)
    return(ALDEx2plot)
  }

prep_AlDEx_data <-
  function(ALDEx_tt,
           ALDEx_effect,
           physeq_obj = Ps_obj_filt_subset,
           sig_level,
           LFC,
           Taxa_rank,
           rare_phyla,
           ...) {
    ALDEx2plot <- data.frame(ALDEx_tt, ALDEx_effect) # merge results
    # group dataframe by OTU, calculate median rel. abundance
    physeq_obj %>%
      transform_sample_counts(., function(x)
        x / sum(x) * 100) %>%
      psmelt() %>%
      group_by(OTU) %>%
      # filter(OTU %in% sig_taxa) %>%
      summarise(baseMean = mean(Abundance)) ->
      baseMean
    
    ALDEx2plot$OTU <- rownames(ALDEx2plot)
    ALDEx2plot %<>% left_join(., baseMean, by = "OTU") # add mean abundance to results table
    ALDEx2plot %<>% cbind(., tax_table(physeq_obj)[taxa_names(physeq_obj) %in% ALDEx2plot$OTU, ], stringsAsFactors = FALSE) # add taxonomy
    # change their name to "Rare"
    
    ALDEx2plot[ALDEx2plot$Phylum %in% rare_phyla, ]$Phylum <-
      'Rare' # rare_phyla is calcuted for the taxa box plots
    ALDEx2plot$Significance <-
      factor("Fail", levels = c("Fail", "Pass")) # define significance factor
    ALDEx2plot$Significance[ALDEx2plot$wi.eBH < sig_level &
                              !is.na(ALDEx2plot$wi.eBH) &
                              abs(ALDEx2plot$effect) > LFC] <-
      "Pass"
    # ALDEx2plot$Significance <- as.factor(sapply(ALDEx2plot$wi.eBH, function(x) if (is.na(x) | x > 0.05) {x <- "Fail"} else {x <- "Pass"}))
    # Rank by taxa abundance
    ALDEx2plot$Phylum %<>%
      factor(., levels = Taxa_rank$Phylum) %>%  # Taxa_rank is calcuted for the taxa box plots
      fct_relevel(., "Rare", after = Inf)
    return(ALDEx2plot)
  }

plot_ALDEx_tax <-
  function(ALDEx2plot,
           OTU_labels = FALSE,
           Taxa = "Phylum",
           Y_val = "effect",
           sig_level = 0.1) {
    pos <- position_jitter(width = 0.3, seed = 1)
    p <-
      ggplot(ALDEx2plot) +
      geom_point(
        aes_string(
          x = Taxa,
          y = Y_val,
          colour = "Significance",
          size = "baseMean"
        ),
        position = pos,
        alpha = 2 / 3,
        stroke = 0
      ) +
      xlab("") +
      ylab(expression(paste("Effect size (lo", g[2], " fold change)"))) +
      # ylab("Fold change") +
      labs(colour = paste("Significance at \n p <", sig_level),
           size = "Mean count (%)") +
      theme_grey(base_size = 18,  base_family = f_name) +
      theme(axis.text.x = element_text(
        angle = 45.0,
        vjust = 1,
        hjust = 1
      )) +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      scale_colour_manual(
        values = c(
          ggpomological:::pomological_base[[7]],
          ggpomological:::pomological_palette[[1]]
        )
      ) +
      scale_size_continuous(range = c(1, 5), breaks = c(1, 2.5, 5, 10))
    
    if (OTU_labels) {
      p <- p + geom_label_repel(
        aes_string(x = Taxa, y = Y_val),
        size = 6,
        label = sub("Seq_([0-9]+)", "\\1", ALDEx2plot[ALDEx2plot$Significance == "Pass", "OTU"]),
        position = pos,
        data = ALDEx2plot[ALDEx2plot$Significance == "Pass", ],
        # nudge_x = 0.4,
        colour = "#4a4a4a",
        label.size = NA,
        alpha = 0.75,
        # fontface = 'bold',
        box.padding = 0.80,
        point.padding = 0.5
        
      )
    }
    return(p)
  }

plot_single_ASV <- function(physeq_obj = Ps_obj_filt_New_Oil_Control, vars2test = "Spill.Treatment", OTU_name = "Seq_8"){
  # adapted from https://link.springer.com/protocol/10.1007%2F978-1-4939-8728-3_10
  require("phyloseq")
  require("ggplot2")
  require("ggpomological")
  # Get Proportion transform
  otuRA <- transform_sample_counts(physeq_obj, function(x) x/sum(x))
  gtab <- tibble(
    Count = get_sample(physeq_obj, i = OTU_name),
    Proportion = get_sample(otuRA, i = OTU_name),
    SAMPLE_ID = unlist(get_variable(physeq_obj, "Joint.sample.name")),
    TREATMENT = get_variable(physeq_obj, vars2test)
  )
  suppressWarnings({
    mgtab <- gather(
      gtab,
      "Count", "Proportion",
      key = "Type",
      value = "Abundance",
    )
  })
  # detect if pairs are both zero
  mgtab %<>% group_by(SAMPLE_ID) %>% mutate(., BothZero := all(Abundance == 0))
  
  # Create a dummy min-value for display
  mgtab %<>% group_by(Type) %>% mutate(Zero = min(Abundance[(Abundance > 0.0)], na.rm = TRUE)/10)
  mgtab %<>% group_by(SAMPLE_ID) %>% mutate(Abundance = replace(Abundance, Abundance == 0.0, Zero[Abundance == 0.0]))

  
  pointSize <- 3
  p <- ggplot(data = mgtab,
             mapping = aes(
               x = TREATMENT, 
               y = Abundance,
               color = TREATMENT,
               shape = TREATMENT)) + 
    facet_wrap(~Type, scales = "free_y") +
    # Not both zero
    geom_point(
      data = filter(mgtab, BothZero == FALSE),
      size = pointSize, 
      alpha = 0.8) +
    # Both Zero
    geom_point(
      data = filter(mgtab, BothZero == TRUE),
      size = pointSize, 
      alpha = 0.8,
      position = position_jitter(width = 0.2, height = 0)) + 
    geom_path(
      data = filter(mgtab, BothZero == FALSE),
      mapping = aes(group = SAMPLE_ID), 
      color = "darkgray", 
      size = 0.25,
      position = position_jitter(width = 0, height = 0.001)) + 
    geom_text(mapping = aes(label = SAMPLE_ID),
              data = filter(mgtab, TREATMENT == levels(gtab$TREATMENT)[1] & Abundance > Zero),
              color = "black",
              size = 2,
              nudge_x = -0.15) +
    scale_colour_manual(values = pom2) +
    # scale_y_sqrt() +
    scale_y_log10() +
    theme_bw() + 
    theme(text = element_text(size = f_size),
          legend.position = "none") +
    scale_size_continuous(range = c(2, 5))
    # ggtitle(paste("Abundance plot for OTU", OTU))
  return(p)
}

plot_top_ASVs <-
  function(physeq_obj = Ps_obj_filt_subset_New_Oil_Control,
           ALDEx_obj = ALDEx2plot_New_Oil_Control,
           vars2test = "Spill.Treatment",
           rank_by = effect,
           Ntop = 12) {
    # adapted from https://link.springer.com/protocol/10.1007%2F978-1-4939-8728-3_10
    require("phyloseq")
    require("ggplot2")
    require("ggpomological")
    
    # Rank OTUs according to rank_by
    quo_rank_by = enquo(rank_by)
    ALDEx_obj %>%
      filter(Significance == "Pass") %>%
      select(OTU, Phylum, baseMean, effect) %>%
      arrange(desc(abs(!!quo_rank_by))) -> OTU_rank
    if (nrow(OTU_rank) < Ntop)
      Ntop <- nrow(OTU_rank)
    OTU_rank %>% .[1:Ntop, ] %>% pull(OTU) -> OTU_names
    
    # Get rel. abund. transform
    otuRA <-
      transform_sample_counts(physeq_obj, function(x)
        x / sum(x) * 100)
    
    if (is.na(OTU_names[1]) | is.null(OTU_names[1])) {
      message("No significant differentially abundant OTUs to display")
    } else {
      if (length(OTU_names) == 1) {
        gtab <- tibble(
          OTU =  as_factor(OTU_names),
          # Count = gather(as_tibble(get_sample(physeq_obj, i = OTU_names)))$value,
          Rel_abundance = gather(as_tibble(get_sample(otuRA, i = OTU_names)))$value,
          SAMPLE_ID = as_factor(rep(unlist(
            get_variable(otuRA, "Joint.sample.name")
          ), Ntop)),
          TREATMENT = as_factor(rep(
            get_variable(otuRA, vars2test), Ntop
          ))
        )
      } else {
        gtab <- tibble(
          OTU =  as_factor(gather(as_tibble(
            get_sample(otuRA, i =
                         OTU_names)
          ))$key),
          # Count = gather(as_tibble(get_sample(otuRA, i = OTU_names)))$value,
          Rel_abundance = gather(as_tibble(get_sample(otuRA, i = OTU_names)))$value,
          SAMPLE_ID = as_factor(rep(unlist(
            get_variable(otuRA, "Joint.sample.name")
          ), Ntop)),
          TREATMENT = as_factor(rep(
            get_variable(otuRA, vars2test), Ntop
          ))
        )
      }
      
      gtab %<>% group_by(SAMPLE_ID, OTU) %>% mutate(., BothZero := all(Rel_abundance == 0))
      
      # Create a dummy min-value for display
      gtab %<>% group_by(OTU) %>% mutate(Zero = min(Rel_abundance[(Rel_abundance > 0.0)], na.rm = TRUE) /
                                           10)
      gtab %<>% group_by(SAMPLE_ID, OTU) %>% mutate(Rel_abundance = replace(Rel_abundance, Rel_abundance == 0.0, Zero[Rel_abundance == 0.0]))
      gtab$SAMPLE_ID
      
      pointSize <- 3
      p <- ggplot(
        data = gtab,
        mapping = aes(
          x = TREATMENT,
          y = Rel_abundance,
          color = TREATMENT,
          shape = TREATMENT
        )
      ) +
        facet_wrap(~ OTU, scales = "free_y") +
        # Not both zero
        geom_point(
          data = filter(gtab, BothZero == FALSE),
          size = pointSize,
          alpha = 0.8
        ) +
        # Both Zero
        geom_point(
          data = filter(gtab, BothZero == TRUE),
          size = pointSize,
          alpha = 0.8,
          position = position_jitter(width = 0.2, height = 0)
        ) +
        geom_path(
          data = filter(gtab, BothZero == FALSE),
          mapping = aes(group = SAMPLE_ID),
          color = "darkgray",
          size = 0.25,
          position = position_jitter(width = 0, height = 0.001)
        ) +
        geom_text(
          mapping = aes(label = SAMPLE_ID),
          data = filter(
            gtab,
            TREATMENT == levels(gtab$TREATMENT)[1] &
              Rel_abundance > Zero
          ),
          color = "black",
          size = 2,
          nudge_x = -0.15
        ) +
        scale_colour_manual(values = pom2) +
        # scale_y_sqrt() +
        scale_y_log10() +
        theme_bw() +
        theme(text = element_text(size = f_size),
              legend.position = "none") +
        scale_size_continuous(range = c(2, 5)) +
        ylab("Abundance (%)")
      # ggtitle(paste("Abundance plot for OTU", OTU))
      return(p)
    }
  }


# Import data, set parameters --------------------------------------------------------------------

Ps_obj <-
  import_biom(...) # if you have your OTU table in biom format
significance <- 0.05

data2test <- t(otu_table(Ps_obj))

ALDEx_comparisons <- list() # here we store the results
ALDEx_comparisons$Comparison <-
  "Control vs. Treatment" # Pairwise comparison string

# subset your phyloseq object to what you're comparing
Ps_obj %>%
  subset_samples(Treatment == "Control" |
                   # "Treatment" is the column in your sample data that contains the relevant labels
                   Treatment == "Treatment") ->
  Ps_obj_subset_pairwise

#  Remove species with prevalence < 10%
Ps_obj_subset_pairwise_s <-
  DropRareSpecies(Ps_obj = Ps_obj_subset_pairwise, prevalence = 0.1)


# Mark rare species -------------------------------------------------------

Ps_obj_subset_glom <- tax_glom(Ps_obj_subset_pairwise_s,
                               "Phylum",
                               NArm = TRUE)
Ps_obj_subset_glom_rel <-
  transform_sample_counts(Ps_obj_subset_glom, function(x)
    x / sum(x))
Ps_obj_subset_glom_rel_DF <- psmelt(Ps_obj_subset_glom_rel)
Ps_obj_subset_glom_rel_DF$Phylum %<>% as.character()

# group dataframe by Phylum, calculate median rel. abundance
Ps_obj_subset_glom_rel_DF %>%
  group_by(Phylum) %>%
  summarise(median = median(Abundance)) ->
  medians

# find Phyla whose rel. abund. is less than 0.5%
Rare_phyla <- medians[medians$median <= 0.005,]$Phylum

# change their name to "Rare"
Ps_obj_subset_glom_rel_DF[Ps_obj_subset_glom_rel_DF$Phylum %in% Rare_phyla,]$Phylum <-
  'Rare'

# re-group
Ps_obj_subset_glom_rel_DF %>%
  group_by(Phylum) %>%
  summarise(Abundance = sum(Abundance)) %>%
  arrange(desc(Abundance)) -> Taxa_rank


# Run ALDEX2 --------------------------------------------------------------

ALDEx2plot_pairwise <- calc_ALDEx2(
  physeq_obj = Ps_obj_subset_pairwise_s,
  vars2test = "Treatment",
  rare_phyla = Rare_phyla,
  sig_level = significance,
  LFC = 0
)
ALDEx_comparisons$Results <-
  ALDEx2plot_pairwise # store results
ALDEx_comparisons$Results$Var1 <-
  str_split(ALDEx_comparisons$Comparison, " vs. ", simplify = TRUE)[1]
ALDEx_comparisons$Results$Var2 <-
  str_split(ALDEx_comparisons$Comparison, " vs. ", simplify = TRUE)[2]
ALDEX_summary <-
  tibble(Label = c(paste0(
    "⬆",
    sum(
      ALDEx_comparisons$Results$effect > 0 &
        ALDEx_comparisons$Results$Significance == "Pass"
    ),
    " ⬇",
    sum(
      ALDEx_comparisons$Results$effect < 0 &
        ALDEx_comparisons$Results$Significance == "Pass"
    ),
    " (",
    nrow(ALDEx_comparisons$Results),
    ")"
  )))

ALDEx2plot_pairwise %>%
  filter(Significance == "Pass") %>%
  select(OTU, baseMean, effect, Phylum, Class, Order, Family, Genus) %>%
  arrange(desc(abs(effect))) ->
  ALDEx2plot_pairwise_results

write.csv(ALDEx2plot_pairwise_results,
          file = paste0("Aldex",
                        "_",
                        paste0(comparison_string, collapse = "_"),
                        ".csv"))

# Plot ALDEX plot
p1 <-
  plot_ALDEx_tax(ALDEx2plot_pairwise,
                 OTU_labels = TRUE,
                 sig_level = significance) +
  # ggtitle(ALDEx_comparisons$Comparisons[j]) +
  geom_text(
    data    = ALDEX_summary,
    mapping = aes(x = Inf, y = Inf, label = Label),
    hjust   = 1.1,
    vjust   = 1.6
  )
p1 <- p1 + labs(title = ALDEx_comparisons$Comparisons[j])
print(p1)

# Plot OTU plots
# GGPlotOTU(Ps_obj_subset_pairwise_s, vars2test = "Spill.Treatment", "Seq_8")
p2 <- plot_top_ASVs(
  Ps_obj = Ps_obj_subset_pairwise_s,
  vars2test = "Spill.Treatment",
  ALDEx_obj = ALDEx2plot_pairwise,
  rank_by = effect,
  Ntop = 12
)
print(p2)

filter(ALDEx_comparisons$Results[[1]], Significance == "Pass")$OTU[common_taxa <-
                                                                     filter(ALDEx_comparisons$Results[[1]], Significance == "Pass")$OTU %in% filter(ALDEx_comparisons$Results[[2]], Significance == "Pass")$OTU]
