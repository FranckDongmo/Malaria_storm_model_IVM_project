
# Fonction pour générer les tableaux LaTeX
generate_latex_tables <- function(IVM_field_dependancy = 0, IVM_Pregnancy = 0,
                                  Prev_mos10 = 1) {
  
  # Charger les données du scénario optimal
  load_file = paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy,
                     "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData")
  load(load_file)
  
  if (length(dim(GainPrev_array_oral)) == 3) {
    GainPrev2_array_oral = GainPrev_array_oral[,,1]
    GainPrev2_array_LAIF = GainPrev_array_LAIF[,,,1]
  } else {
    GainPrev2_array_oral = GainPrev_array_oral
    GainPrev2_array_LAIF = GainPrev_array_LAIF
  }
  
  # Table 1: Differences R-gain pour plot_formulations_comparison()
  ##########################
  # Générer tableau LaTeX 
  ##########################
  VectPropIVM_target = c(0.5, 0.7, 0.9)
  VectNumber_of_cycle = c(1, 2, 3, 4)
  
  format_value <- function(val) {
    sprintf("%.2f", val)
  }
  
  # Créer le tableau complet pour tous les PropIVM
  LatexTable = "\\begin{table}[h!]\n"
  LatexTable = paste0(LatexTable, "\\scriptsize\n")
  LatexTable = paste0(LatexTable, "\\centering\n")
  LatexTable = paste0(LatexTable, "\\caption{Performance differences (\\%) between formulations for different TP coverage levels}\n")
  LatexTable = paste0(LatexTable, "\\label{tab:comparisons_complete}\n")
  LatexTable = paste0(LatexTable, "\\begin{tabular}{l|ccc}\n")
  LatexTable = paste0(LatexTable, "\\toprule\n")
  
  # Section 1: Oral vs Oral
  LatexTable = paste0(LatexTable, "\\multicolumn{4}{c}{\\textbf{Oral vs Oral}} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  LatexTable = paste0(LatexTable, "\\textbf{Comparison} & \\textbf{50\\%} & \\textbf{70\\%} & \\textbf{90\\%} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  
  # Comparaisons Oral vs Oral
  oral_comparisons = list(
    list(name="Δ(RII[S],BOH)", idx1=2, idx2=1),
    list(name="Δ(RII[K],BOH)", idx1=3, idx2=1),
    list(name="Δ(RII[K],RII[S])", idx1=3, idx2=2)
  )
  
  for (comp in oral_comparisons) {
    # Formater pour LaTeX
    comp_latex = gsub("Δ\\(", "$\\\\Delta$(", comp$name)
    comp_latex = gsub("RII\\[S\\]", "RII$_{S}$", comp_latex)
    comp_latex = gsub("RII\\[K\\]", "RII$_{K}$", comp_latex)
    comp_latex = paste0(comp_latex, ")")
    
    row = comp_latex
    
    # Valeurs pour chaque PropIVM
    for (PropIVM_target in VectPropIVM_target) {
      IndexPropIVM = which.min(abs(VectPropIVM - PropIVM_target))
      
      val = (GainPrev2_array_oral[comp$idx1, IndexPropIVM] - 
               GainPrev2_array_oral[comp$idx2, IndexPropIVM]) * 100
      
      row = paste0(row, " & ", format_value(val))
    }
    
    LatexTable = paste0(LatexTable, row, " \\\\ \n")
  }
  
  # Section 2: Injectable vs Oral
  LatexTable = paste0(LatexTable, "\\midrule\n")
  LatexTable = paste0(LatexTable, "\\end{tabular}\n")
  LatexTable = paste0(LatexTable, "\\\\[0.5em]\n")
  
  # Nouveau tableau pour Injectable vs Oral
  LatexTable = paste0(LatexTable, "\\begin{tabular}{l|ccc|ccc|ccc|ccc}\n")
  LatexTable = paste0(LatexTable, "\\toprule\n")
  LatexTable = paste0(LatexTable, "\\multicolumn{13}{c}{\\textbf{Injectable vs Oral}} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  LatexTable = paste0(LatexTable, "\\multirow{2}{*}{\\textbf{Comparison}} & \\multicolumn{3}{c|}{\\textbf{$n_c=1$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=2$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=3$}} & \\multicolumn{3}{c}{\\textbf{$n_c=4$}} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\cmidrule{2-13}\n")
  LatexTable = paste0(LatexTable, " & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  
  # Comparaisons Injectable vs Oral
  laif_strategies = c("0.6-vk5", "1.0", "0.6-kis")
  oral_strategies = c("BOH", "RII[S]", "RII[K]")
  
  for (laif_idx in 1:3) {
    for (oral_idx in 1:3) {
      laif_name = laif_strategies[laif_idx]
      oral_name = oral_strategies[oral_idx]
      
      # Formater les noms
      if (laif_name == "1.0") {
        laif_display = "1.0"
      } else {
        laif_display = laif_name
      }
      
      if (oral_name == "BOH") {
        oral_display = "BOH"
      } else if (oral_name == "RII[S]") {
        oral_display = "RII$_{S}$"
      } else {
        oral_display = "RII$_{K}$"
      }
      
      comp_latex = paste0("$\\Delta$(", laif_display, ",", oral_display, ")")
      row = comp_latex
      
      # Valeurs pour chaque cycle et chaque PropIVM
      for (nc in VectNumber_of_cycle) {
        for (PropIVM_target in VectPropIVM_target) {
          IndexPropIVM = which.min(abs(VectPropIVM - PropIVM_target))
          
          val = (GainPrev2_array_LAIF[laif_idx, IndexPropIVM, nc] - 
                   GainPrev2_array_oral[oral_idx, IndexPropIVM]) * 100
          
          row = paste0(row, " & ", format_value(val))
        }
      }
      
      LatexTable = paste0(LatexTable, row, " \\\\ \n")
    }
  }
  
  # Section 3: Injectable vs Injectable
  LatexTable = paste0(LatexTable, "\\midrule\n")
  LatexTable = paste0(LatexTable, "\\end{tabular}\n")
  LatexTable = paste0(LatexTable, "\\\\[0.5em]\n")
  
  # Nouveau tableau pour Injectable vs Injectable
  LatexTable = paste0(LatexTable, "\\begin{tabular}{l|ccc|ccc|ccc|ccc}\n")
  LatexTable = paste0(LatexTable, "\\toprule\n")
  LatexTable = paste0(LatexTable, "\\multicolumn{13}{c}{\\textbf{Injectable vs Injectable}} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  LatexTable = paste0(LatexTable, "\\multirow{2}{*}{\\textbf{Comparison}} & \\multicolumn{3}{c|}{\\textbf{$n_c=1$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=2$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=3$}} & \\multicolumn{3}{c}{\\textbf{$n_c=4$}} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\cmidrule{2-13}\n")
  LatexTable = paste0(LatexTable, " & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  
  # Comparaisons Injectable vs Injectable
  laif_comparisons = list(
    list(name="Δ(1.0,0.6-vk5)", idx1=2, idx2=1),
    list(name="Δ(1.0,0.6-kis)", idx1=2, idx2=3),
    list(name="Δ(0.6-kis,0.6-vk5)", idx1=3, idx2=1)
  )
  
  for (comp_info in laif_comparisons) {
    comp_latex = gsub("Δ\\(", "$\\\\Delta$(", comp_info$name)
    comp_latex = paste0(comp_latex, ")")
    
    row = comp_latex
    
    # Valeurs pour chaque cycle et chaque PropIVM
    for (nc in VectNumber_of_cycle) {
      for (PropIVM_target in VectPropIVM_target) {
        IndexPropIVM = which.min(abs(VectPropIVM - PropIVM_target))
        
        val = (GainPrev2_array_LAIF[comp_info$idx1, IndexPropIVM, nc] - 
                 GainPrev2_array_LAIF[comp_info$idx2, IndexPropIVM, nc]) * 100
        
        row = paste0(row, " & ", format_value(val))
      }
    }
    
    LatexTable = paste0(LatexTable, row, " \\\\ \n")
  }
  
  # Finaliser le dernier tableau
  LatexTable = paste0(LatexTable, "\\bottomrule\n")
  LatexTable = paste0(LatexTable, "\\end{tabular}\n")
  
  # Notes
  LatexTable = paste0(LatexTable, "\\\\[0.5em]\n")
  LatexTable = paste0(LatexTable, "\\begin{tablenotes}\n")
  LatexTable = paste0(LatexTable, "\\scriptsize\n")
  LatexTable = paste0(LatexTable, "\\item \\textbf{Note:} Coverage percentages refer to TP coverage levels (50\\%, 70\\%, 90\\%). ")
  LatexTable = paste0(LatexTable, "$n_c$ indicates number of cycles for injectable formulations. ")
  LatexTable = paste0(LatexTable, "Values show percentage point differences in performance.\n")
  LatexTable = paste0(LatexTable, "\\end{tablenotes}\n")
  LatexTable = paste0(LatexTable, "\\end{table}")
  
  filename1 = paste0("Comparisons_Delta_PairF_table_field", IVM_field_dependancy,
                     "_Preg", IVM_Pregnancy, "_Mos", Name_Prev_mos(Prev_mos10), ".txt")
  writeLines(LatexTable, filename1)
  
  # Table 2: GainPrev Values pour plot_gainprev_barplot()
  ##########################
  # Générer tableaux LaTeX avec valeurs GainPrev
  ##########################
  format_value2 <- function(val) {
    sprintf("%.2f", val * 100)  # Multiplier par 100 pour avoir des pourcentages
  }
  
  # Créer le tableau complet
  LatexTable2 = "\\begin{table}[h!]\n"
  LatexTable2 = paste0(LatexTable2, "\\scriptsize\n")
  LatexTable2 = paste0(LatexTable2, "\\centering\n")
  LatexTable2 = paste0(LatexTable2, "\\caption{GainPrev values (\\%) for different formulations and TP coverage levels}\n")
  LatexTable2 = paste0(LatexTable2, "\\label{tab:gainprev_values}\n")
  
  # TABLEAU 1: FORMULATIONS ORALES
  LatexTable2 = paste0(LatexTable2, "\\begin{tabular}{l|ccc}\n")
  LatexTable2 = paste0(LatexTable2, "\\toprule\n")
  LatexTable2 = paste0(LatexTable2, "\\multicolumn{4}{c}{\\textbf{Oral Formulations}} \\\\ \n")
  LatexTable2 = paste0(LatexTable2, "\\midrule\n")
  LatexTable2 = paste0(LatexTable2, "\\textbf{Formulation} & \\textbf{50\\%} & \\textbf{70\\%} & \\textbf{90\\%} \\\\ \n")
  LatexTable2 = paste0(LatexTable2, "\\midrule\n")
  
  oral_names = c("BOH", "RII$_{S}$", "RII$_{K}$")
  
  for (oral_idx in 1:3) {
    row = oral_names[oral_idx]
    
    # Valeurs pour chaque PropIVM
    for (PropIVM_target in VectPropIVM_target) {
      IndexPropIVM = which.min(abs(VectPropIVM - PropIVM_target))
      
      val = GainPrev2_array_oral[oral_idx, IndexPropIVM]
      row = paste0(row, " & ", format_value2(val))
    }
    
    LatexTable2 = paste0(LatexTable2, row, " \\\\ \n")
  }
  
  LatexTable2 = paste0(LatexTable2, "\\bottomrule\n")
  LatexTable2 = paste0(LatexTable2, "\\end{tabular}\n")
  LatexTable2 = paste0(LatexTable2, "\\\\[1em]\n")
  
  # TABLEAU 2: FORMULATIONS INJECTABLES
  LatexTable2 = paste0(LatexTable2, "\\begin{tabular}{l|ccc|ccc|ccc|ccc}\n")
  LatexTable2 = paste0(LatexTable2, "\\toprule\n")
  LatexTable2 = paste0(LatexTable2, "\\multicolumn{13}{c}{\\textbf{Injectable Formulations}} \\\\ \n")
  LatexTable2 = paste0(LatexTable2, "\\midrule\n")
  LatexTable2 = paste0(LatexTable2, "\\multirow{2}{*}{\\textbf{Formulation}} & \\multicolumn{3}{c|}{\\textbf{$n_c=1$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=2$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=3$}} & \\multicolumn{3}{c}{\\textbf{$n_c=4$}} \\\\ \n")
  LatexTable2 = paste0(LatexTable2, "\\cmidrule{2-13}\n")
  LatexTable2 = paste0(LatexTable2, " & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} \\\\ \n")
  LatexTable2 = paste0(LatexTable2, "\\midrule\n")
  
  # Noms des formulations injectables
  laif_names = c("0.6-vk5", "1.0", "0.6-kis")
  
  for (laif_idx in 1:3) {
    row = laif_names[laif_idx]
    
    # Valeurs pour chaque cycle et chaque PropIVM
    for (nc in VectNumber_of_cycle) {
      for (PropIVM_target in VectPropIVM_target) {
        IndexPropIVM = which.min(abs(VectPropIVM - PropIVM_target))
        
        val = GainPrev2_array_LAIF[laif_idx, IndexPropIVM, nc]
        row = paste0(row, " & ", format_value2(val))
      }
    }
    
    LatexTable2 = paste0(LatexTable2, row, " \\\\ \n")
  }
  
  LatexTable2 = paste0(LatexTable2, "\\bottomrule\n")
  LatexTable2 = paste0(LatexTable2, "\\end{tabular}\n")
  
  # Notes 
  LatexTable2 = paste0(LatexTable2, "\\\\[0.5em]\n")
  LatexTable2 = paste0(LatexTable2, "\\begin{tablenotes}\n")
  LatexTable2 = paste0(LatexTable2, "\\scriptsize\n")
  LatexTable2 = paste0(LatexTable2, "\\item \\textbf{Note:} Coverage percentages refer to TP coverage levels (50\\%, 70\\%, 90\\%). ")
  LatexTable2 = paste0(LatexTable2, "$n_c$ indicates number of cycles for injectable formulations. ")
  LatexTable2 = paste0(LatexTable2, "Values represent GainPrev performance metrics in percentage.\n")
  LatexTable2 = paste0(LatexTable2, "\\end{tablenotes}\n")
  LatexTable2 = paste0(LatexTable2, "\\end{table}")
  
  filename2 = paste0("GainPrev_Values_table_field", IVM_field_dependancy,
                     "_Preg", IVM_Pregnancy, "_Mos", Name_Prev_mos(Prev_mos10), ".txt")
  writeLines(LatexTable2, filename2)
  
  cat("Tables LaTeX générées:\n", filename1, "\n", filename2, "\n")
}

# Utilisation
generate_latex_tables(IVM_field_dependancy = 0, IVM_Pregnancy = 0, Prev_mos10 = 1)

# Fonction pour générer les tableaux de comparaison (opt_vs_base, pregnant vs chilbearing etc.)
generate_comparison_latex_table <- function(comparison_type, IVM_Pregnancy_ref = 0,
                                            Prev_mos10_ref = 1,
                                            IVM_field_dependancy_ref = 0) {
  
  VectPropIVM_target = c(0.5, 0.7, 0.9)
  VectNumber_of_cycle = c(1, 2, 3, 4)
  
  format_value <- function(val) {
    sprintf("%.2f", val)
  }
  
  # Charger les données selon le type de comparaison
  if (comparison_type == "opt_vs_base") {
    # Charger optimal
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy_ref, "_Pregnant", IVM_Pregnancy_ref,
                "Prev_mos", Name_Prev_mos(Prev_mos10_ref), ".RData"))
    GainPrevopt_array_oral = GainPrev_array_oral
    GainPrevopt_array_LAIF = GainPrev_array_LAIF
    VectPropIVM_opt = VectPropIVM
    
    # Charger base
    load(paste0("ResultatsPrev_IhBase_field_dep_", IVM_field_dependancy_ref, "_Pregnant", IVM_Pregnancy_ref,
                "Prev_mos", Name_Prev_mos(Prev_mos10_ref), ".RData"))
    GainPrevBase_array_oral = GainPrev_array_oral
    GainPrevBase_array_LAIF = GainPrev_array_LAIF
    VectPropIVM = VectPropIVM_opt
    
    comparison_label = "Optimal vs Baseline"
    
  } else if (comparison_type == "field_dep") {
    # Comparer IVM_field_dependancy 0 vs 1
    load(paste0("ResultatsPrev_Ihopt_field_dep_0_Pregnant", IVM_Pregnancy_ref,
                "Prev_mos", Name_Prev_mos(Prev_mos10_ref), ".RData"))
    GainPrev1_array_oral = GainPrev_array_oral
    GainPrev1_array_LAIF = GainPrev_array_LAIF
    VectPropIVM_1 = VectPropIVM
    
    load(paste0("ResultatsPrev_Ihopt_field_dep_1_Pregnant", IVM_Pregnancy_ref,
                "Prev_mos", Name_Prev_mos(Prev_mos10_ref), ".RData"))
    GainPrev2_array_oral = GainPrev_array_oral
    GainPrev2_array_LAIF = GainPrev_array_LAIF
    VectPropIVM = VectPropIVM_1
    
    comparison_label = "Field dependency (0 vs 1)"
    
  } else if (comparison_type == "prev_mos") {
    # Comparer Prev_mos10 0 vs 1
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy_ref, "_Pregnant", IVM_Pregnancy_ref,
                "Prev_mos10", ".RData"))
    GainPrev1_array_oral = GainPrev_array_oral
    GainPrev1_array_LAIF = GainPrev_array_LAIF
    VectPropIVM_1 = VectPropIVM
    
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy_ref, "_Pregnant", IVM_Pregnancy_ref,
                "Prev_mos5", ".RData"))
    GainPrev2_array_oral = GainPrev_array_oral
    GainPrev2_array_LAIF = GainPrev_array_LAIF
    VectPropIVM = VectPropIVM_1
    
    comparison_label = "Mosquito prevalence (5% vs 10%)"
    
  } else if (comparison_type == "pregnancy") {
    # Comparer IVM_Pregnancy 0 vs 1
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy_ref, "_Pregnant1",
                "Prev_mos", Name_Prev_mos(Prev_mos10_ref), ".RData"))
    GainPrev1_array_oral = GainPrev_array_oral
    GainPrev1_array_LAIF = GainPrev_array_LAIF
    VectPropIVM_1 = VectPropIVM
    
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy_ref, "_Pregnant0",
                "Prev_mos", Name_Prev_mos(Prev_mos10_ref), ".RData"))
    GainPrev2_array_oral = GainPrev_array_oral
    GainPrev2_array_LAIF = GainPrev_array_LAIF
    VectPropIVM = VectPropIVM_1
    
    comparison_label = "Pregnancy (with vs without)"
    
  } else if (comparison_type == "deltoptim_vs_15") {
    # Comparer DeltaOptim vs DeltaOptim_15
    load(paste0("ResultatsPrev_Ihopt_field_dep_", IVM_field_dependancy_ref, "_Pregnant", IVM_Pregnancy_ref,
                "Prev_mos", Name_Prev_mos(Prev_mos10_ref), ".RData"))
    GainPrev1_array_oral = GainPrev_array_oral
    GainPrev1_array_LAIF = GainPrev_array_LAIF
    GainPrev2_array_oral = GainPrev_array_oral_15
    GainPrev2_array_LAIF = GainPrev_array_LAIF_15
    
    comparison_label = "All ages vs Under 5"
  }
  
  # Créer le tableau LaTeX
  LatexTable = "\\begin{table}[h!]\n"
  LatexTable = paste0(LatexTable, "\\scriptsize\n")
  LatexTable = paste0(LatexTable, "\\centering\n")
  LatexTable = paste0(LatexTable, "\\caption{Performance comparison: ", comparison_label, "}\n")
  LatexTable = paste0(LatexTable, "\\label{tab:comparison_", comparison_type, "}\n")
  
  # Section 1: Formulations orales
  LatexTable = paste0(LatexTable, "\\begin{tabular}{l|ccc}\n")
  LatexTable = paste0(LatexTable, "\\toprule\n")
  LatexTable = paste0(LatexTable, "\\multicolumn{4}{c}{\\textbf{Oral Formulations}} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  LatexTable = paste0(LatexTable, "\\textbf{Formulation} & \\textbf{50\\%} & \\textbf{70\\%} & \\textbf{90\\%} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  
  oral_names = c("BOH", "RII$_{S}$", "RII$_{K}$")
  
  for (oral_idx in 1:3) {
    row = oral_names[oral_idx]
    
    for (PropIVM_target in VectPropIVM_target) {
      IndexPropIVM = which.min(abs(VectPropIVM - PropIVM_target))
      
      if (comparison_type == "opt_vs_base") {
        if (length(dim(GainPrevopt_array_oral)) == 2) {
          val = (GainPrevopt_array_oral[oral_idx, IndexPropIVM] - 
                   GainPrevBase_array_oral[oral_idx, IndexPropIVM]) * 100
        } else {
          val = (GainPrevopt_array_oral[oral_idx, IndexPropIVM, 1] - 
                   GainPrevBase_array_oral[oral_idx, IndexPropIVM, 1]) * 100
        }
      } else {
        if (length(dim(GainPrev1_array_oral)) == 2) {
          val = (GainPrev1_array_oral[oral_idx, IndexPropIVM] - 
                   GainPrev2_array_oral[oral_idx, IndexPropIVM]) * 100
        } else {
          val = (GainPrev1_array_oral[oral_idx, IndexPropIVM, 1] - 
                   GainPrev2_array_oral[oral_idx, IndexPropIVM, 1]) * 100
        }
      }
      
      row = paste0(row, " & ", format_value(val))
    }
    
    LatexTable = paste0(LatexTable, row, " \\\\ \n")
  }
  
  LatexTable = paste0(LatexTable, "\\bottomrule\n")
  LatexTable = paste0(LatexTable, "\\end{tabular}\n")
  LatexTable = paste0(LatexTable, "\\\\[1em]\n")
  
  # Section 2: Formulations injectables
  LatexTable = paste0(LatexTable, "\\begin{tabular}{l|ccc|ccc|ccc|ccc}\n")
  LatexTable = paste0(LatexTable, "\\toprule\n")
  LatexTable = paste0(LatexTable, "\\multicolumn{13}{c}{\\textbf{Injectable Formulations}} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  LatexTable = paste0(LatexTable, "\\multirow{2}{*}{\\textbf{Formulation}} & \\multicolumn{3}{c|}{\\textbf{$n_c=1$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=2$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=3$}} & \\multicolumn{3}{c}{\\textbf{$n_c=4$}} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\cmidrule{2-13}\n")
  LatexTable = paste0(LatexTable, " & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} \\\\ \n")
  LatexTable = paste0(LatexTable, "\\midrule\n")
  
  laif_names = c("0.6-vk5", "1.0", "0.6-kis")
  
  for (laif_idx in 1:3) {
    row = laif_names[laif_idx]
    
    for (nc in VectNumber_of_cycle) {
      for (PropIVM_target in VectPropIVM_target) {
        IndexPropIVM = which.min(abs(VectPropIVM - PropIVM_target))
        
        if (comparison_type == "opt_vs_base") {
          if (length(dim(GainPrevopt_array_LAIF)) == 3) {
            val = (GainPrevopt_array_LAIF[laif_idx, IndexPropIVM, nc] - 
                     GainPrevBase_array_LAIF[laif_idx, IndexPropIVM, nc]) * 100
          } else {
            val = (GainPrevopt_array_LAIF[laif_idx, IndexPropIVM, nc, 1] - 
                     GainPrevBase_array_LAIF[laif_idx, IndexPropIVM, nc, 1]) * 100
          }
        } else {
          if (length(dim(GainPrev1_array_LAIF)) == 3) {
            val = (GainPrev1_array_LAIF[laif_idx, IndexPropIVM, nc] - 
                     GainPrev2_array_LAIF[laif_idx, IndexPropIVM, nc]) * 100
          } else {
            val = (GainPrev1_array_LAIF[laif_idx, IndexPropIVM, nc, 1] - 
                     GainPrev2_array_LAIF[laif_idx, IndexPropIVM, nc, 1]) * 100
          }
        }
        
        row = paste0(row, " & ", format_value(val))
      }
    }
    
    LatexTable = paste0(LatexTable, row, " \\\\ \n")
  }
  
  LatexTable = paste0(LatexTable, "\\bottomrule\n")
  LatexTable = paste0(LatexTable, "\\end{tabular}\n")
  
  # Notes
  LatexTable = paste0(LatexTable, "\\\\[0.5em]\n")
  LatexTable = paste0(LatexTable, "\\begin{tablenotes}\n")
  LatexTable = paste0(LatexTable, "\\scriptsize\n")
  LatexTable = paste0(LatexTable, "\\item \\textbf{Note:} Coverage percentages refer to TP coverage levels (50\\%, 70\\%, 90\\%). ")
  LatexTable = paste0(LatexTable, "$n_c$ indicates number of cycles for injectable formulations. ")
  
  if (comparison_type == "opt_vs_base") {
    LatexTable = paste0(LatexTable, "Values show difference between optimal and baseline scenarios (percentage points).\n")
  } else if (comparison_type == "field_dep") {
    LatexTable = paste0(LatexTable, "Values show difference between field dependency scenarios (percentage points).\n")
  } else if (comparison_type == "prev_mos") {
    LatexTable = paste0(LatexTable, "Values show difference between mosquito prevalence levels (percentage points).\n")
  } else if (comparison_type == "pregnancy") {
    LatexTable = paste0(LatexTable, "Values show difference with and without pregnancy inclusion (percentage points).\n")
  } else if (comparison_type == "deltoptim_vs_15") {
    LatexTable = paste0(LatexTable, "Values show difference between all ages and under-5 metrics (percentage points).\n")
  }
  
  LatexTable = paste0(LatexTable, "\\end{tablenotes}\n")
  LatexTable = paste0(LatexTable, "\\end{table}")
  
  # Sauvegarder le fichier
  filename = paste0("Comparison_", comparison_type, "_table.txt")
  writeLines(LatexTable, filename)
  
  cat("Tableau LaTeX généré:", filename, "\n")
}

# Utilisation
 generate_comparison_latex_table("opt_vs_base", IVM_Pregnancy_ref = 0, Prev_mos10_ref = 1, IVM_field_dependancy_ref = 0)
 generate_comparison_latex_table("field_dep", IVM_Pregnancy_ref = 0, Prev_mos10_ref = 1)
 generate_comparison_latex_table("prev_mos", IVM_Pregnancy_ref = 0, IVM_field_dependancy_ref = 0)
 generate_comparison_latex_table("pregnancy", Prev_mos10_ref = 1, IVM_field_dependancy_ref = 0)
 generate_comparison_latex_table("deltoptim_vs_15", IVM_Pregnancy_ref = 0, Prev_mos10_ref = 1, IVM_field_dependancy_ref = 0)
 
 
 # Fonction pour générer le tableau LaTeX des données saisonnières
 # Fonction pour générer le tableau LaTeX des données saisonnières (temps 250 uniquement)
 generate_seasonality_table_250 <- function(IVM_field_dependancy = 0, IVM_Pregnancy = 0,
                                            Prev_mos10 = 1) {
   
   # Charger les données saisonnières
   load_file = paste0("ResultatsPrev_Seasonal_field_dep_", IVM_field_dependancy,
                      "_Pregnant", IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".RData")
   load(load_file)
   
   # Paramètres
   VectPropIVM_target = c(0.5, 0.7, 0.9)
   t_idx = 2  # Index pour t_begin_Camp = 250 (position 2 dans c(220, 250, 280))
   VectNumber_of_cycle = c(1, 2, 3, 4)
   
   # Fonction de formatage
   format_value <- function(val) {
     sprintf("%.2f", val * 100)  # Convertir en pourcentage
   }
   
   # Créer le tableau LaTeX
   LatexTable = "\\begin{table}[h!]\n"
   LatexTable = paste0(LatexTable, "\\scriptsize\n")
   LatexTable = paste0(LatexTable, "\\centering\n")
   LatexTable = paste0(LatexTable, "\\caption{Seasonal GainPrev values (\\%) - Campaign start: Day 250}\n")
   LatexTable = paste0(LatexTable, "\\label{tab:seasonal_gainprev_250}\n")
   
   # TABLEAU 1: FORMULATIONS ORALES
   LatexTable = paste0(LatexTable, "\\begin{tabular}{l|ccc}\n")
   LatexTable = paste0(LatexTable, "\\toprule\n")
   LatexTable = paste0(LatexTable, "\\multicolumn{4}{c}{\\textbf{Oral Formulations}} \\\\ \n")
   LatexTable = paste0(LatexTable, "\\midrule\n")
   LatexTable = paste0(LatexTable, "\\textbf{Formulation} & \\textbf{50\\%} & \\textbf{70\\%} & \\textbf{90\\%} \\\\ \n")
   LatexTable = paste0(LatexTable, "\\midrule\n")
   
   # Noms des formulations orales
   oral_names = c("BOH", "RII$_{S}$", "RII$_{K}$")
   
   for (oral_idx in 1:3) {
     row = oral_names[oral_idx]
     
     # Valeurs pour chaque PropIVM
     for (prop_idx in 1:length(VectPropIVM_target)) {
       val = GainPrev_array_oral_Seas[oral_idx, prop_idx, t_idx]
       row = paste0(row, " & ", format_value(val))
     }
     
     LatexTable = paste0(LatexTable, row, " \\\\ \n")
   }
   
   LatexTable = paste0(LatexTable, "\\bottomrule\n")
   LatexTable = paste0(LatexTable, "\\end{tabular}\n")
   LatexTable = paste0(LatexTable, "\\\\[1em]\n")
   
   # TABLEAU 2: FORMULATIONS INJECTABLES
   LatexTable = paste0(LatexTable, "\\begin{tabular}{l|ccc|ccc|ccc|ccc}\n")
   LatexTable = paste0(LatexTable, "\\toprule\n")
   LatexTable = paste0(LatexTable, "\\multicolumn{13}{c}{\\textbf{Injectable Formulations}} \\\\ \n")
   LatexTable = paste0(LatexTable, "\\midrule\n")
   LatexTable = paste0(LatexTable, "\\multirow{2}{*}{\\textbf{Formulation}} & \\multicolumn{3}{c|}{\\textbf{$n_c=1$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=2$}} & \\multicolumn{3}{c|}{\\textbf{$n_c=3$}} & \\multicolumn{3}{c}{\\textbf{$n_c=4$}} \\\\ \n")
   LatexTable = paste0(LatexTable, "\\cmidrule{2-13}\n")
   LatexTable = paste0(LatexTable, " & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} & \\textbf{50} & \\textbf{70} & \\textbf{90} \\\\ \n")
   LatexTable = paste0(LatexTable, "\\midrule\n")
   
   # Noms des formulations injectables
   laif_names = c("0.6-vk5", "1.0", "0.6-kis")
   
   for (laif_idx in 1:3) {
     row = laif_names[laif_idx]
     
     # Valeurs pour chaque cycle et chaque PropIVM
     for (nc in VectNumber_of_cycle) {
       for (prop_idx in 1:length(VectPropIVM_target)) {
         val = GainPrev_array_LAIF_Seas[laif_idx, prop_idx, nc, t_idx]
         row = paste0(row, " & ", format_value(val))
       }
     }
     
     LatexTable = paste0(LatexTable, row, " \\\\ \n")
   }
   
   LatexTable = paste0(LatexTable, "\\bottomrule\n")
   LatexTable = paste0(LatexTable, "\\end{tabular}\n")
   
   # Notes explicatives
   LatexTable = paste0(LatexTable, "\\\\[0.5em]\n")
   LatexTable = paste0(LatexTable, "\\begin{tablenotes}\n")
   LatexTable = paste0(LatexTable, "\\scriptsize\n")
   LatexTable = paste0(LatexTable, "\\item \\textbf{Note:} Coverage percentages refer to TP coverage levels (50\\%, 70\\%, 90\\%). ")
   LatexTable = paste0(LatexTable, "$n_c$ indicates number of cycles for injectable formulations. ")
   LatexTable = paste0(LatexTable, "Campaign start at day 250 (mid-season). ")
   LatexTable = paste0(LatexTable, "Values represent GainPrev performance metrics in percentage.\n")
   LatexTable = paste0(LatexTable, "\\end{tablenotes}\n")
   LatexTable = paste0(LatexTable, "\\end{table}")
   
   # Sauvegarder le fichier
   filename = paste0("Seasonal_GainPrev_t250_table_field", IVM_field_dependancy,
                     "_Preg", IVM_Pregnancy, "_Mos", Name_Prev_mos(Prev_mos10), ".txt")
   writeLines(LatexTable, filename)
   
   cat("Tableau LaTeX saisonnier (t=250) généré: ", filename, "\n")
 }
 
 # Utilisation
 generate_seasonality_table_250(IVM_field_dependancy = 0, IVM_Pregnancy = 0, Prev_mos10 = 1)

