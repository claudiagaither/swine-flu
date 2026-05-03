## Functions for US_flu script

## part three: continuous trait diffusion ----
## ── 3. Helper: build a single clade XML ─────────────────────────────────────
## v2: enforced monophyly + uncorrelated log-normal relaxed clock (UCLD)
##     following Ebola/WNV phylogeography approach (Dudas & Rambaut 2014,
##     Pybus et al. 2012, Lemey et al. 2010)
build_clade_xml <- function(clade_name, clade_tips, tip_meta_df,
                            taxon_dates, seq_seqs,
                            loc_lat, loc_lon,
                            chain_length, log_every, save_every) {
  
  ## clean file-safe name
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", clade_name)
  file_stem <- paste0("1990_", safe_name)
  
  ## match tips to template (some tree tips may not be in XML if naming differs)
  keep <- clade_tips[clade_tips %in% names(taxon_dates) &
                       clade_tips %in% names(seq_seqs)]
  
  ## look up state for each tip
  tip_states <- tip_meta_df$state[match(keep, tip_meta_df$sequence_name)]
  tip_lat    <- loc_lat[tip_states]
  tip_lon    <- loc_lon[tip_states]
  
  ## BEAST continuous phylogeography requires EVERY taxon to have a location.
  ## Drop any tips without coordinates (no state, or state not in locations table).
  has_coords <- !is.na(tip_lat) & !is.na(tip_lon)
  n_dropped  <- sum(!has_coords)
  if (n_dropped > 0) {
    dropped_states <- tip_states[!has_coords]
    message("  Dropping ", n_dropped, " tips without coordinates: ",
            paste(unique(na.omit(dropped_states)), collapse = ", "),
            if (any(is.na(dropped_states))) paste0(" + ", sum(is.na(dropped_states)), " with NA state"))
  }
  keep       <- keep[has_coords]
  tip_states <- tip_states[has_coords]
  tip_lat    <- tip_lat[has_coords]
  tip_lon    <- tip_lon[has_coords]
  
  if (length(keep) < 3) {
    message("  Skipping clade '", clade_name, "': only ", length(keep),
            " tips with coordinates (need >=3).")
    return(invisible(NULL))
  }
  
  n_with <- length(keep)
  
  ## nchar from first sequence
  nchar_aln <- nchar(seq_seqs[keep[1]])
  
  ## ── taxa block ──
  ## Every taxon gets lat/lon attrs -- BEAST will crash without them
  taxa_lines <- vapply(seq_along(keep), function(i) {
    tid  <- keep[i]
    dval <- taxon_dates[tid]
    
    sprintf('\t\t<taxon id="%s">\n\t\t\t<date value="%s" direction="forwards" units="years"/>\n\t\t\t<attr name="location">%s %s</attr>\n\t\t</taxon>',
            tid, dval, tip_lat[i], tip_lon[i])
  }, character(1))
  
  ## ── alignment block ──
  aln_lines <- vapply(keep, function(tid) {
    sprintf('\t\t<sequence>\n\t\t\t<taxon idref="%s"/>\n\t\t\t%s\n\t\t</sequence>',
            tid, seq_seqs[tid])
  }, character(1))
  
  ## ── assemble full XML ──
  xml_out <- paste0(
    '<?xml version="1.0" standalone="yes"?>\n',
    '\n',
    '<!-- Clade-specific BEAST XML: ', clade_name, ' -->\n',
    '<!-- Generated programmatically from 1990_USflu_tree1.xml -->\n',
    '<!-- ntax=', length(keep), '  (', n_dropped, ' tips dropped for missing coordinates) -->\n',
    '<!-- v2: enforced monophyly + UCLD relaxed clock -->\n',
    '<beast version="1.10.4">\n',
    '\n',
    
    ## ── TAXA ──
    '\t<!-- ntax=', length(keep), ' -->\n',
    '\t<taxa id="taxa">\n',
    paste(taxa_lines, collapse = "\n"), '\n',
    '\t</taxa>\n\n',
    
    ## ── MONOPHYLY: taxon set for the clade ──
    '\t<!-- Enforce monophyly: all taxa belong to same clade -->\n',
    '\t<taxa id="clade">\n',
    paste(vapply(keep, function(tid) {
      sprintf('\t\t<taxon idref="%s"/>', tid)
    }, character(1)), collapse = "\n"), '\n',
    '\t</taxa>\n\n',
    
    ## ── ALIGNMENT ──
    '\t<!-- ntax=', length(keep), ' nchar=', nchar_aln, ' -->\n',
    '\t<alignment id="alignment" dataType="nucleotide">\n',
    paste(aln_lines, collapse = "\n"), '\n',
    '\t</alignment>\n\n',
    
    ## ── PATTERNS ──
    '\t<patterns id="patterns" from="1" strip="false">\n',
    '\t\t<alignment idref="alignment"/>\n',
    '\t</patterns>\n\n',
    
    ## ── COALESCENT: exponential growth ──
    '\t<exponentialGrowth id="exponential" units="years">\n',
    '\t\t<populationSize>\n',
    '\t\t\t<parameter id="exponential.popSize" value="1.0" lower="0.0"/>\n',
    '\t\t</populationSize>\n',
    '\t\t<growthRate>\n',
    '\t\t\t<parameter id="exponential.growthRate" value="0.0"/>\n',
    '\t\t</growthRate>\n',
    '\t</exponentialGrowth>\n\n',
    
    ## ── STARTING TREE ──
    '\t<coalescentSimulator id="startingTree">\n',
    '\t\t<taxa idref="taxa"/>\n',
    '\t\t<exponentialGrowth idref="exponential"/>\n',
    '\t</coalescentSimulator>\n\n',
    
    ## ── TREE MODEL ──
    '\t<treeModel id="treeModel">\n',
    '\t\t<coalescentTree idref="startingTree"/>\n',
    '\t\t<rootHeight>\n',
    '\t\t\t<parameter id="treeModel.rootHeight"/>\n',
    '\t\t</rootHeight>\n',
    '\t\t<nodeHeights internalNodes="true">\n',
    '\t\t\t<parameter id="treeModel.internalNodeHeights"/>\n',
    '\t\t</nodeHeights>\n',
    '\t\t<nodeHeights internalNodes="true" rootNode="true">\n',
    '\t\t\t<parameter id="treeModel.allInternalNodeHeights"/>\n',
    '\t\t</nodeHeights>\n',
    '\t</treeModel>\n\n',
    
    '\t<treeLengthStatistic id="treeLength">\n',
    '\t\t<treeModel idref="treeModel"/>\n',
    '\t</treeLengthStatistic>\n\n',
    
    '\t<tmrcaStatistic id="age(root)" absolute="true">\n',
    '\t\t<treeModel idref="treeModel"/>\n',
    '\t</tmrcaStatistic>\n\n',
    
    ## ── MONOPHYLY STATISTIC ──
    '\t<!-- Monophyly statistic: constrains clade to be monophyletic -->\n',
    '\t<tmrcaStatistic id="tmrca(clade)" monophyletic="true">\n',
    '\t\t<mrca>\n',
    '\t\t\t<taxa idref="clade"/>\n',
    '\t\t</mrca>\n',
    '\t\t<treeModel idref="treeModel"/>\n',
    '\t</tmrcaStatistic>\n\n',
    
    ## ── COALESCENT LIKELIHOOD ──
    '\t<coalescentLikelihood id="coalescent">\n',
    '\t\t<model>\n',
    '\t\t\t<exponentialGrowth idref="exponential"/>\n',
    '\t\t</model>\n',
    '\t\t<populationTree>\n',
    '\t\t\t<treeModel idref="treeModel"/>\n',
    '\t\t</populationTree>\n',
    '\t</coalescentLikelihood>\n\n',
    
    ## ── CLOCK: uncorrelated log-normal relaxed clock (UCLD) ──
    '\t<!-- Uncorrelated log-normal relaxed clock (Drummond et al. 2006) -->\n',
    '\t<discretizedBranchRates id="branchRates">\n',
    '\t\t<treeModel idref="treeModel"/>\n',
    '\t\t<distribution>\n',
    '\t\t\t<logNormalDistributionModel meanInRealSpace="true">\n',
    '\t\t\t\t<mean>\n',
    '\t\t\t\t\t<parameter id="ucld.mean" value="0.004" lower="0.0"/>\n',
    '\t\t\t\t</mean>\n',
    '\t\t\t\t<stdev>\n',
    '\t\t\t\t\t<parameter id="ucld.stdev" value="0.3333" lower="0.0"/>\n',
    '\t\t\t\t</stdev>\n',
    '\t\t\t</logNormalDistributionModel>\n',
    '\t\t</distribution>\n',
    '\t\t<rateCategories>\n',
    '\t\t\t<parameter id="branchRates.categories"/>\n',
    '\t\t</rateCategories>\n',
    '\t</discretizedBranchRates>\n\n',
    
    '\t<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">\n',
    '\t\t<treeModel idref="treeModel"/>\n',
    '\t\t<discretizedBranchRates idref="branchRates"/>\n',
    '\t</rateStatistic>\n\n',
    
    '\t<rateStatistic id="coefficientOfVariation" name="coefficientOfVariation"\n',
    '\t\t\tmode="coefficientOfVariation" internal="true" external="true">\n',
    '\t\t<treeModel idref="treeModel"/>\n',
    '\t\t<discretizedBranchRates idref="branchRates"/>\n',
    '\t</rateStatistic>\n\n',
    
    '\t<rateCovarianceStatistic id="covariance" name="covariance">\n',
    '\t\t<treeModel idref="treeModel"/>\n',
    '\t\t<discretizedBranchRates idref="branchRates"/>\n',
    '\t</rateCovarianceStatistic>\n\n',
    
    ## ── SUBSTITUTION MODEL: HKY ──
    '\t<HKYModel id="hky">\n',
    '\t\t<frequencies>\n',
    '\t\t\t<frequencyModel dataType="nucleotide">\n',
    '\t\t\t\t<frequencies>\n',
    '\t\t\t\t\t<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>\n',
    '\t\t\t\t</frequencies>\n',
    '\t\t\t</frequencyModel>\n',
    '\t\t</frequencies>\n',
    '\t\t<kappa>\n',
    '\t\t\t<parameter id="kappa" value="2.0" lower="0.0"/>\n',
    '\t\t</kappa>\n',
    '\t</HKYModel>\n\n',
    
    '\t<siteModel id="siteModel">\n',
    '\t\t<substitutionModel>\n',
    '\t\t\t<HKYModel idref="hky"/>\n',
    '\t\t</substitutionModel>\n',
    '\t</siteModel>\n\n',
    
    ## ── TREE LIKELIHOOD ──
    '\t<treeDataLikelihood id="treeLikelihood" useAmbiguities="false">\n',
    '\t\t<partition>\n',
    '\t\t\t<patterns idref="patterns"/>\n',
    '\t\t\t<siteModel idref="siteModel"/>\n',
    '\t\t</partition>\n',
    '\t\t<treeModel idref="treeModel"/>\n',
    '\t\t<discretizedBranchRates idref="branchRates"/>\n',
    '\t</treeDataLikelihood>\n\n',
    
    ## ── CONTINUOUS TRAIT: multivariate diffusion for lat/lon ──
    '\t<!-- Multivariate diffusion model for continuous phylogeography -->\n',
    '\t<multivariateDiffusionModel id="location.diffusionModel">\n',
    '\t\t<precisionMatrix>\n',
    '\t\t\t<matrixParameter id="location.precision">\n',
    '\t\t\t\t<parameter id="location.precision.col1" value="0.05 0.002"/>\n',
    '\t\t\t\t<parameter id="location.precision.col2" value="0.002 0.05"/>\n',
    '\t\t\t</matrixParameter>\n',
    '\t\t</precisionMatrix>\n',
    '\t</multivariateDiffusionModel>\n\n',
    
    ## ── CONTINUOUS TRAIT LIKELIHOOD ──
    '\t<multivariateTraitLikelihood id="location.traitLikelihood" traitName="location"\n',
    '\t\t\tuseTreeLength="true" scaleByTime="true"\n',
    '\t\t\treportAsSum="true" reciprocalRates="true"\n',
    '\t\t\tintegrateInternalTraits="true">\n',
    '\t\t<multivariateDiffusionModel idref="location.diffusionModel"/>\n',
    '\t\t<treeModel idref="treeModel"/>\n',
    '\t\t<traitParameter>\n',
    '\t\t\t<parameter id="leaf.location"/>\n',
    '\t\t</traitParameter>\n',
    '\t\t<conjugateRootPrior>\n',
    '\t\t\t<meanParameter>\n',
    '\t\t\t\t<parameter value="0.0 0.0"/>\n',
    '\t\t\t</meanParameter>\n',
    '\t\t\t<priorSampleSize>\n',
    '\t\t\t\t<parameter value="0.001"/>\n',
    '\t\t\t</priorSampleSize>\n',
    '\t\t</conjugateRootPrior>\n',
    '\t</multivariateTraitLikelihood>\n\n',
    
    ## ── OPERATORS ──
    '\t<operators id="operators" optimizationSchedule="default">\n',
    
    ## substitution model operators
    '\t\t<scaleOperator scaleFactor="0.75" weight="1">\n',
    '\t\t\t<parameter idref="kappa"/>\n',
    '\t\t</scaleOperator>\n',
    '\t\t<deltaExchange delta="0.01" weight="1">\n',
    '\t\t\t<parameter idref="frequencies"/>\n',
    '\t\t</deltaExchange>\n',
    
    ## UCLD operators
    '\t\t<!-- UCLD relaxed clock operators -->\n',
    '\t\t<scaleOperator scaleFactor="0.75" weight="3">\n',
    '\t\t\t<parameter idref="ucld.mean"/>\n',
    '\t\t</scaleOperator>\n',
    '\t\t<scaleOperator scaleFactor="0.75" weight="3">\n',
    '\t\t\t<parameter idref="ucld.stdev"/>\n',
    '\t\t</scaleOperator>\n',
    
    ## up-down: scale tree heights up while scaling ucld.mean down (critical for mixing)
    '\t\t<upDownOperator scaleFactor="0.75" weight="3">\n',
    '\t\t\t<up>\n',
    '\t\t\t\t<parameter idref="treeModel.allInternalNodeHeights"/>\n',
    '\t\t\t</up>\n',
    '\t\t\t<down>\n',
    '\t\t\t\t<parameter idref="ucld.mean"/>\n',
    '\t\t\t</down>\n',
    '\t\t</upDownOperator>\n',
    
    ## rate category operators (swap & random walk on integer categories)
    '\t\t<swapOperator size="1" weight="10" autoOptimize="false">\n',
    '\t\t\t<parameter idref="branchRates.categories"/>\n',
    '\t\t</swapOperator>\n',
    '\t\t<randomWalkIntegerOperator windowSize="1" weight="10">\n',
    '\t\t\t<parameter idref="branchRates.categories"/>\n',
    '\t\t</randomWalkIntegerOperator>\n',
    '\t\t<uniformIntegerOperator weight="10">\n',
    '\t\t\t<parameter idref="branchRates.categories"/>\n',
    '\t\t</uniformIntegerOperator>\n',
    
    ## tree topology operators
    '\t\t<!-- Tree operators -->\n',
    '\t\t<subtreeSlide size="1.0" gaussian="true" weight="30">\n',
    '\t\t\t<treeModel idref="treeModel"/>\n',
    '\t\t</subtreeSlide>\n',
    '\t\t<narrowExchange weight="30">\n',
    '\t\t\t<treeModel idref="treeModel"/>\n',
    '\t\t</narrowExchange>\n',
    '\t\t<wideExchange weight="3">\n',
    '\t\t\t<treeModel idref="treeModel"/>\n',
    '\t\t</wideExchange>\n',
    '\t\t<wilsonBalding weight="3">\n',
    '\t\t\t<treeModel idref="treeModel"/>\n',
    '\t\t</wilsonBalding>\n',
    '\t\t<scaleOperator scaleFactor="0.75" weight="3">\n',
    '\t\t\t<parameter idref="treeModel.rootHeight"/>\n',
    '\t\t</scaleOperator>\n',
    '\t\t<uniformOperator weight="30">\n',
    '\t\t\t<parameter idref="treeModel.internalNodeHeights"/>\n',
    '\t\t</uniformOperator>\n',
    
    ## coalescent operators
    '\t\t<scaleOperator scaleFactor="0.75" weight="3">\n',
    '\t\t\t<parameter idref="exponential.popSize"/>\n',
    '\t\t</scaleOperator>\n',
    '\t\t<randomWalkOperator windowSize="1.0" weight="3">\n',
    '\t\t\t<parameter idref="exponential.growthRate"/>\n',
    '\t\t</randomWalkOperator>\n',
    
    ## continuous trait operators
    '\t\t<!-- Continuous trait operators -->\n',
    '\t\t<precisionGibbsOperator weight="2">\n',
    '\t\t\t<multivariateTraitLikelihood idref="location.traitLikelihood"/>\n',
    '\t\t\t<multivariateWishartPrior id="location.precisionPrior" df="2">\n',
    '\t\t\t\t<scaleMatrix>\n',
    '\t\t\t\t\t<matrixParameter>\n',
    '\t\t\t\t\t\t<parameter value="1.0 0.0"/>\n',
    '\t\t\t\t\t\t<parameter value="0.0 1.0"/>\n',
    '\t\t\t\t\t</matrixParameter>\n',
    '\t\t\t\t</scaleMatrix>\n',
    '\t\t\t\t<data>\n',
    '\t\t\t\t\t<parameter idref="location.precision"/>\n',
    '\t\t\t\t</data>\n',
    '\t\t\t</multivariateWishartPrior>\n',
    '\t\t</precisionGibbsOperator>\n',
    '\t</operators>\n\n',
    
    ## ── MCMC ──
    '\t<mcmc id="mcmc" chainLength="', format(chain_length, scientific = FALSE), '" autoOptimize="true"\n',
    '\t\t  operatorAnalysis="', file_stem, '.ops.txt">\n',
    '\t\t<joint id="joint">\n',
    '\t\t\t<prior id="prior">\n',
    
    ## -- priors --
    '\t\t\t\t<logNormalPrior mu="1.0" sigma="1.25" offset="0.0">\n',
    '\t\t\t\t\t<parameter idref="kappa"/>\n',
    '\t\t\t\t</logNormalPrior>\n',
    '\t\t\t\t<dirichletPrior alpha="1.0" sumsTo="1.0">\n',
    '\t\t\t\t\t<parameter idref="frequencies"/>\n',
    '\t\t\t\t</dirichletPrior>\n',
    
    ## UCLD priors: CTMC reference prior on ucld.mean, exponential on ucld.stdev
    '\t\t\t\t<!-- UCLD relaxed clock priors -->\n',
    '\t\t\t\t<ctmcScalePrior>\n',
    '\t\t\t\t\t<ctmcScale>\n',
    '\t\t\t\t\t\t<parameter idref="ucld.mean"/>\n',
    '\t\t\t\t\t</ctmcScale>\n',
    '\t\t\t\t\t<treeModel idref="treeModel"/>\n',
    '\t\t\t\t</ctmcScalePrior>\n',
    '\t\t\t\t<exponentialPrior mean="0.3333" offset="0.0">\n',
    '\t\t\t\t\t<parameter idref="ucld.stdev"/>\n',
    '\t\t\t\t</exponentialPrior>\n',
    
    ## coalescent priors
    '\t\t\t\t<oneOnXPrior>\n',
    '\t\t\t\t\t<parameter idref="exponential.popSize"/>\n',
    '\t\t\t\t</oneOnXPrior>\n',
    '\t\t\t\t<laplacePrior mean="0.0" scale="1.0">\n',
    '\t\t\t\t\t<parameter idref="exponential.growthRate"/>\n',
    '\t\t\t\t</laplacePrior>\n',
    '\t\t\t\t<coalescentLikelihood idref="coalescent"/>\n',
    
    ## continuous trait prior
    '\t\t\t\t<multivariateWishartPrior idref="location.precisionPrior"/>\n',
    
    ## monophyly enforcement (hard constraint via booleanLikelihood)
    '\t\t\t\t<!-- Enforce monophyly of the clade -->\n',
    '\t\t\t\t<booleanLikelihood>\n',
    '\t\t\t\t\t<monophylyStatistic id="monophyly(clade)">\n',
    '\t\t\t\t\t\t<mrca>\n',
    '\t\t\t\t\t\t\t<taxa idref="clade"/>\n',
    '\t\t\t\t\t\t</mrca>\n',
    '\t\t\t\t\t\t<treeModel idref="treeModel"/>\n',
    '\t\t\t\t\t</monophylyStatistic>\n',
    '\t\t\t\t</booleanLikelihood>\n',
    
    '\t\t\t</prior>\n',
    '\t\t\t<likelihood id="likelihood">\n',
    '\t\t\t\t<treeDataLikelihood idref="treeLikelihood"/>\n',
    '\t\t\t\t<multivariateTraitLikelihood idref="location.traitLikelihood"/>\n',
    '\t\t\t</likelihood>\n',
    '\t\t</joint>\n',
    '\t\t<operators idref="operators"/>\n\n',
    
    ## ── SCREEN LOG ──
    '\t\t<log id="screenLog" logEvery="', log_every, '">\n',
    '\t\t\t<column label="Joint" dp="4" width="12">\n',
    '\t\t\t\t<joint idref="joint"/>\n',
    '\t\t\t</column>\n',
    '\t\t\t<column label="Prior" dp="4" width="12">\n',
    '\t\t\t\t<prior idref="prior"/>\n',
    '\t\t\t</column>\n',
    '\t\t\t<column label="Likelihood" dp="4" width="12">\n',
    '\t\t\t\t<likelihood idref="likelihood"/>\n',
    '\t\t\t</column>\n',
    '\t\t\t<column label="age(root)" sf="6" width="12">\n',
    '\t\t\t\t<tmrcaStatistic idref="age(root)"/>\n',
    '\t\t\t</column>\n',
    '\t\t\t<column label="ucld.mean" sf="6" width="12">\n',
    '\t\t\t\t<parameter idref="ucld.mean"/>\n',
    '\t\t\t</column>\n',
    '\t\t</log>\n\n',
    
    ## ── FILE LOG ──
    '\t\t<log id="fileLog" logEvery="', log_every, '" fileName="', file_stem, '.log.txt" overwrite="false">\n',
    '\t\t\t<joint idref="joint"/>\n',
    '\t\t\t<prior idref="prior"/>\n',
    '\t\t\t<likelihood idref="likelihood"/>\n',
    '\t\t\t<parameter idref="treeModel.rootHeight"/>\n',
    '\t\t\t<tmrcaStatistic idref="age(root)"/>\n',
    '\t\t\t<tmrcaStatistic idref="tmrca(clade)"/>\n',
    '\t\t\t<treeLengthStatistic idref="treeLength"/>\n',
    '\t\t\t<parameter idref="exponential.popSize"/>\n',
    '\t\t\t<parameter idref="exponential.growthRate"/>\n',
    '\t\t\t<parameter idref="kappa"/>\n',
    '\t\t\t<parameter idref="frequencies"/>\n',
    '\t\t\t<parameter idref="ucld.mean"/>\n',
    '\t\t\t<parameter idref="ucld.stdev"/>\n',
    '\t\t\t<rateStatistic idref="meanRate"/>\n',
    '\t\t\t<rateStatistic idref="coefficientOfVariation"/>\n',
    '\t\t\t<rateCovarianceStatistic idref="covariance"/>\n',
    '\t\t\t<treeDataLikelihood idref="treeLikelihood"/>\n',
    '\t\t\t<coalescentLikelihood idref="coalescent"/>\n',
    '\t\t\t<multivariateTraitLikelihood idref="location.traitLikelihood"/>\n',
    '\t\t\t<matrixParameter idref="location.precision"/>\n',
    '\t\t</log>\n\n',
    
    ## ── TREE LOG ──
    '\t\t<logTree id="treeFileLog" logEvery="', log_every, '" nexusFormat="true"\n',
    '\t\t\t\t fileName="', file_stem, '.trees.txt" sortTranslationTable="true">\n',
    '\t\t\t<treeModel idref="treeModel"/>\n',
    '\t\t\t<trait name="rate" tag="rate">\n',
    '\t\t\t\t<discretizedBranchRates idref="branchRates"/>\n',
    '\t\t\t</trait>\n',
    '\t\t\t<multivariateTraitLikelihood idref="location.traitLikelihood"/>\n',
    '\t\t\t<joint idref="joint"/>\n',
    '\t\t</logTree>\n\n',
    
    '\t</mcmc>\n\n',
    
    '\t<report>\n',
    '\t\t<property name="timer">\n',
    '\t\t\t<mcmc idref="mcmc"/>\n',
    '\t\t</property>\n',
    '\t</report>\n\n',
    
    '</beast>\n'
  )
  
  return(list(xml = xml_out, file_stem = file_stem, n_tips = length(keep),
              n_dropped = n_dropped))
}
