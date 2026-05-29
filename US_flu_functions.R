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
  tip_states <- tip_meta_df$state[base::match(keep, tip_meta_df$sequence_name)]
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


## part four: time-series models -----

# Functions for modeling the monthly prevalence of an INDIVIDUAL clade
# (e.g. "2010.1like", "1990.4.a") in a state as an ARIMAX outcome, with
# weighted climate variables as predictors (trending temps detrended), an
# AIC comparison across predictors, a (p,q) sensitivity sweep at a
# programmatically chosen d, and residual checks.
#
# Pure function definitions only -- source() this from your analysis script.
# Required packages (load them in your analysis script):
#   tidyverse, lubridate, forecast, tseries, gt, scales
#
# Expected inputs (same objects as part 3A):
#   state_prev        : year_month, state, clade, count, total, prevalence
#   climate_state_wt  : STATE_NAME, year, month, tmax..abs_humidity

# Default climate variables and which ones get detrended (the temps that
# showed a significant linear trend in part 3A).
clade_clim_vars    <- c("tmax", "tmin", "tavg", "ppt", "vap", "ws", "abs_humidity")
clade_detrend_vars <- c("tmax", "tmin", "tavg")


# predictor_columns(): names of the columns fed to models = detrended temps
# (*_dt) + the remaining climate vars used raw.

predictor_columns <- function(clim_vars = clade_clim_vars,
                              detrend_vars = clade_detrend_vars) {
  c(paste0(detrend_vars, "_dt"), base::setdiff(clim_vars, detrend_vars))
}


# build_clade_series(): monthly prevalence of ONE clade in ONE state on a
# complete monthly grid. Distinguishes true zeros (state sequenced that month
# but not this clade -> prev = 0) from missing months (no sequencing -> NA).
#   prev      = n_clade / n_total
#   emp_logit = log((n_clade + 0.5) / (n_total - n_clade + 0.5))   [0/1-safe]

build_clade_series <- function(state_prev, st, focal_clade, date_start, date_end) {
  grid <- tibble::tibble(date  = seq(date_start, date_end, by = "month"),
                         state = st)
  
  totals <- state_prev %>%
    dplyr::filter(state == st) %>%
    dplyr::distinct(year_month, total) %>%
    dplyr::transmute(date = as.Date(year_month), n_total = total)
  
  focal <- state_prev %>%
    dplyr::filter(state == st, clade == focal_clade) %>%
    dplyr::transmute(date = as.Date(year_month), n_clade = count)
  
  grid %>%
    dplyr::left_join(totals, by = "date") %>%
    dplyr::left_join(focal,  by = "date") %>%
    dplyr::mutate(
      # sequenced month but clade absent -> true zero; no sequencing -> NA
      n_clade   = dplyr::if_else(!is.na(n_total) & is.na(n_clade), 0L, n_clade),
      prev      = n_clade / n_total,
      emp_logit = log((n_clade + 0.5) / (n_total - n_clade + 0.5))) %>%
    dplyr::arrange(date)
}


# build_state_clade_data(): join the clade outcome to that state's climate and
# detrend the flagged temperature predictors (residual of value ~ time index).
# Returns the per-(state, clade) modeling frame.

build_state_clade_data <- function(state_prev, climate_state_wt, st, focal_clade,
                                   date_start, date_end,
                                   clim_vars = clade_clim_vars,
                                   detrend_vars = clade_detrend_vars) {
  out <- build_clade_series(state_prev, st, focal_clade, date_start, date_end)
  
  cl <- climate_state_wt %>%
    dplyr::mutate(date = as.Date(sprintf("%d-%02d-01", year, month))) %>%
    dplyr::rename(state = STATE_NAME) %>%
    dplyr::filter(state == st, date >= date_start, date <= date_end) %>%
    dplyr::select(date, dplyr::all_of(clim_vars))
  
  df <- out %>%
    dplyr::left_join(cl, by = "date") %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(t_index = as.numeric(difftime(date, min(date), units = "days")))
  
  for (v in detrend_vars) {
    df[[paste0(v, "_dt")]] <-
      as.numeric(resid(lm(df[[v]] ~ df$t_index, na.action = na.exclude)))
  }
  attr(df, "state") <- st
  attr(df, "clade") <- focal_clade
  df
}


# get_outcome(): pull the outcome as a ts on the chosen scale.
#   "logit" (recommended; variance-stabilized) or "prev" (raw proportion).

get_outcome <- function(df, scale = c("logit", "prev"), freq = 12) {
  scale <- match.arg(scale)
  v <- if (scale == "logit") df$emp_logit else df$prev
  ts(as.numeric(v), frequency = freq)
}


# choose_differencing(): pick d (KPSS via ndiffs) and seasonal D (nsdiffs),
# with optional manual overrides. Gaps are interpolated ONLY for the tests;
# model fits keep the real NAs. Returns d, D, test p-values, n_obs.

choose_differencing <- function(y, seasonal, freq, force_d = NULL, force_D = NULL) {
  yt   <- ts(as.numeric(y), frequency = freq)
  yt_i <- na.omit(forecast::na.interp(yt))
  
  d <- if (!is.null(force_d)) as.integer(force_d) else forecast::ndiffs(yt_i, test = "kpss")
  D <- if (!is.null(force_D)) as.integer(force_D)
  else if (seasonal && freq > 1)
    tryCatch(forecast::nsdiffs(ts(yt_i, frequency = freq)), error = function(e) 0L)
  else 0L
  
  list(d = as.integer(d), D = as.integer(D),
       kpss_p   = suppressWarnings(tseries::kpss.test(yt_i)$p.value),
       adf_p    = suppressWarnings(tseries::adf.test(yt_i)$p.value),
       n_obs    = sum(!is.na(as.numeric(y))),
       forced_d = !is.null(force_d),
       forced_D = !is.null(force_D))
}


# fit_arimax_by_predictor(): for each climate predictor (entered singly as
# xreg) plus a no-predictor baseline, let auto.arima pick the order and record
# AIC/AICc/BIC. Used to rank which climate variable best explains the clade.

fit_arimax_by_predictor <- function(df, st, focal_clade,
                                    predictors = predictor_columns(),
                                    scale = "logit", seasonal = TRUE, freq = 12) {
  y <- get_outcome(df, scale, freq)
  
  one <- function(label, xreg) {
    fit <- tryCatch(
      forecast::auto.arima(y, xreg = xreg, seasonal = seasonal,
                           stepwise = TRUE, approximation = FALSE),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ord <- forecast::arimaorder(fit)
    tibble::tibble(
      state = st, clade = focal_clade, predictor = label,
      order = sprintf("(%d,%d,%d)", ord["p"], ord["d"], ord["q"]),
      seasonal = if (seasonal && length(ord) >= 6)
        sprintf("(%d,%d,%d)[%d]", ord["P"], ord["D"], ord["Q"], freq) else "-",
      n_obs = sum(!is.na(y)),
      AIC = round(fit$aic, 2), AICc = round(fit$aicc, 2), BIC = round(fit$bic, 2))
  }
  
  dplyr::bind_rows(
    one("none (baseline)", NULL),
    purrr::map(predictors, ~ one(.x, df[[.x]]))
  ) %>%
    dplyr::mutate(delta_AIC = round(AIC - min(AIC, na.rm = TRUE), 2)) %>%
    dplyr::arrange(AIC)
}


# arimax_sensitivity_grid(): sweep (p, q) at a FIXED d and seasonal (P, D, Q),
# so AICc/ΔAICc are comparable across every cell. Returns the grid.

arimax_sensitivity_grid <- function(df, predictor, scale, seasonal, freq,
                                    p_range, q_range, d, seasonal_PQ, D) {
  y  <- get_outcome(df, scale, freq)
  xr <- df[[predictor]]
  s_order <- c(seasonal_PQ[1], D, seasonal_PQ[2])
  
  tidyr::expand_grid(p = p_range, q = q_range) %>%
    purrr::pmap_dfr(function(p, q) {
      fit <- suppressWarnings(tryCatch(
        forecast::Arima(y, order = c(p, d, q),
                        seasonal = list(order = s_order, period = freq),
                        xreg = xr, method = "ML"),
        error = function(e) NULL))
      tibble::tibble(
        p, d, q,
        n_par     = if (is.null(fit)) NA_integer_ else length(fit$coef),
        AIC       = if (is.null(fit)) NA_real_ else round(fit$aic, 2),
        AICc      = if (is.null(fit)) NA_real_ else round(fit$aicc, 2),
        BIC       = if (is.null(fit)) NA_real_ else round(fit$bic, 2),
        converged = !is.null(fit))
    }) %>%
    dplyr::mutate(delta_AICc = round(AICc - min(AICc, na.rm = TRUE), 2)) %>%
    dplyr::arrange(AICc)
}


# run_clade_arimax(): orchestrator for ONE (state, clade). Builds data, ranks
# predictors, chooses d/D (flagged via message), runs the (p,q) sweep on the
# best predictor, and returns everything in a list.

run_clade_arimax <- function(state_prev, climate_state_wt, st, focal_clade,
                             date_start, date_end,
                             scale = "logit", seasonal = TRUE,
                             clim_vars = clade_clim_vars,
                             detrend_vars = clade_detrend_vars,
                             p_range = 0:3, q_range = 0:3, seasonal_PQ = c(1, 0),
                             force_d = NULL, force_D = NULL, verbose = TRUE) {
  freq <- if (seasonal) 12 else 1
  
  df <- build_state_clade_data(state_prev, climate_state_wt, st, focal_clade,
                               date_start, date_end, clim_vars, detrend_vars)
  preds <- predictor_columns(clim_vars, detrend_vars)
  
  pred_res  <- fit_arimax_by_predictor(df, st, focal_clade, preds, scale, seasonal, freq)
  best_pred <- pred_res %>%
    dplyr::filter(predictor != "none (baseline)") %>%
    dplyr::slice_min(AIC, n = 1, with_ties = FALSE) %>%
    dplyr::pull(predictor)
  
  dsel <- choose_differencing(get_outcome(df, scale, freq), seasonal, freq, force_d, force_D)
  
  if (verbose) message(sprintf(
    paste0(">> %s ~ %s [clade outcome] | n_obs = %d\n",
           ">> Chosen d = %d (KPSS p = %.3f, ADF p = %.3f)%s | seasonal D = %d%s\n",
           ">> d and D fixed for the sweep -> AICc comparable across all cells."),
    st, focal_clade, dsel$n_obs,
    dsel$d, dsel$kpss_p, dsel$adf_p, if (dsel$forced_d) " [forced]" else "",
    dsel$D, if (dsel$forced_D) " [forced]" else ""))
  
  sens <- arimax_sensitivity_grid(df, best_pred, scale, seasonal, freq,
                                  p_range, q_range, dsel$d, seasonal_PQ, dsel$D)
  
  list(state = st, clade = focal_clade, data = df,
       predictor_results = pred_res, best_predictor = best_pred,
       differencing = dsel, sensitivity = sens,
       params = list(scale = scale, seasonal = seasonal, freq = freq,
                     seasonal_PQ = seasonal_PQ))
}


# render_predictor_table(): shaded gt of the per-predictor AIC comparison.

render_predictor_table <- function(run) {
  r <- run$predictor_results
  r %>%
    dplyr::select(predictor, order, seasonal, n_obs, AIC, delta_AIC) %>%
    gt::gt() %>%
    gt::data_color(columns = delta_AIC,
                   fn = scales::col_numeric(c("#1a9850", "#ffffbf", "#d73027"),
                                            domain = range(r$delta_AIC, na.rm = TRUE),
                                            na.color = "grey85")) %>%
    gt::tab_header(
      title = gt::md(sprintf("**Predictor comparison — %s ~ climate (%s)**",
                             run$state, run$clade)),
      subtitle = "One climate predictor per model; temps detrended. Ranked by AIC.") %>%
    gt::cols_label(predictor = "Predictor (xreg)", delta_AIC = "ΔAIC")
}


# render_sensitivity_table(): shaded gt of the (p,q) sweep at fixed d/D.

render_sensitivity_table <- function(run) {
  s <- run$sensitivity; d <- run$differencing; pq <- run$params$seasonal_PQ
  fixed_label <- sprintf("d = %d, seasonal (%d,%d,%d)[%d] — fixed",
                         d$d, pq[1], d$D, pq[2], run$params$freq)
  s %>%
    dplyr::transmute(order = sprintf("(%d,%d,%d)", p, d, q),
                     n_par, AIC, AICc, BIC, delta_AICc, converged) %>%
    gt::gt() %>%
    gt::data_color(columns = delta_AICc,
                   fn = scales::col_numeric(c("#1a9850", "#ffffbf", "#d73027"),
                                            domain = range(s$delta_AICc, na.rm = TRUE),
                                            na.color = "grey85")) %>%
    gt::tab_header(
      title = gt::md(sprintf("**ARIMA(p,d,q) sensitivity — %s ~ %s (%s)**",
                             run$state, run$best_predictor, run$clade)),
      subtitle = gt::md(sprintf("%s • KPSS p = %.3f, ADF p = %.3f • ΔAICc comparable across all rows",
                                fixed_label, d$kpss_p, d$adf_p))) %>%
    gt::cols_label(n_par = "k", delta_AICc = "ΔAICc")
}


# render_sensitivity_heatmap(): ΔAICc over the (p,q) plane at the fixed d.

render_sensitivity_heatmap <- function(run) {
  s <- run$sensitivity; d <- run$differencing; pq <- run$params$seasonal_PQ
  fixed_label <- sprintf("d = %d, seasonal (%d,%d,%d)[%d]",
                         d$d, pq[1], d$D, pq[2], run$params$freq)
  ggplot2::ggplot(s, ggplot2::aes(factor(p), factor(q), fill = delta_AICc)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(
      label = ifelse(is.na(delta_AICc), "×", sprintf("%.1f", delta_AICc))), size = 3) +
    ggplot2::scale_fill_viridis_c(direction = -1, na.value = "grey85", name = "ΔAICc") +
    ggplot2::labs(
      title = sprintf("Order sensitivity: %s ~ %s (%s)",
                      run$state, run$best_predictor, run$clade),
      subtitle = sprintf("%s • × = did not converge", fixed_label),
      x = "AR order (p)", y = "MA order (q)") +
    ggplot2::theme_minimal(base_size = 11)
}


# check_clade_model(): refit the best (min-AICc) order from the sweep and run
# residual diagnostics. Returns the fitted model invisibly.

check_clade_model <- function(run) {
  best <- run$sensitivity %>% dplyr::filter(converged) %>%
    dplyr::slice_min(AICc, n = 1, with_ties = FALSE)
  p <- run$params; d <- run$differencing
  y  <- get_outcome(run$data, p$scale, p$freq)
  xr <- run$data[[run$best_predictor]]
  fit <- forecast::Arima(
    y, order = c(best$p, best$d, best$q),
    seasonal = list(order = c(p$seasonal_PQ[1], d$D, p$seasonal_PQ[2]), period = p$freq),
    xreg = xr, method = "ML")
  print(forecast::checkresiduals(fit))
  invisible(fit)
}



# Part four (cont.): predictor significance for the selected ARIMAX models

# Consume the `run` objects from run_clade_arimax() and produce inference on
# the SELECTED climate predictor. Analysis is exploratory: the outcome is clade
# prevalence among publicly deposited sequences (convenience sample), so this
# is association within sampled sequences, not population clade prevalence.
# ARIMA coefficients use asymptotic (Wald) normal inference -> z-tests, not
# t-tests. SE = sqrt(diag(var.coef)). Outcome is logit, so exp(coef) is an OR.

# refit_best_model(): refit the min-AICc order from the sweep with a NAMED
# one-column xreg, so the climate coefficient carries the predictor's name
# (an unnamed numeric vector would be labelled "xreg"). Shared helper so the
# diagnostics and the significance table describe the exact same fitted model.

refit_best_model <- function(run) {
  best <- run$sensitivity %>%
    dplyr::filter(converged) %>%
    dplyr::slice_min(AICc, n = 1, with_ties = FALSE)
  if (nrow(best) == 0L)
    stop("No converged model in the sensitivity grid for ",
         run$state, " / ", run$clade)

  p <- run$params; d <- run$differencing
  y  <- get_outcome(run$data, p$scale, p$freq)
  xr <- matrix(run$data[[run$best_predictor]], ncol = 1,
               dimnames = list(NULL, run$best_predictor))

  fit <- forecast::Arima(
    y, order = c(best$p, best$d, best$q),
    seasonal = list(order = c(p$seasonal_PQ[1], d$D, p$seasonal_PQ[2]),
                    period = p$freq),
    xreg = xr, method = "ML")
  attr(fit, "best_order") <- best
  fit
}


# arimax_wald(): Wald table (estimate / SE / z / p / 95% CI) for ALL terms of
# a fitted Arima object. var.coef names align with coef names; we reindex SEs
# defensively in case any coefficient was held fixed.

arimax_wald <- function(fit) {
  est <- fit$coef
  se  <- sqrt(diag(fit$var.coef))[names(est)]
  z   <- est / se
  tibble::tibble(
    term      = names(est),
    estimate  = as.numeric(est),
    std_error = as.numeric(se),
    z         = as.numeric(z),
    p_value   = 2 * stats::pnorm(-abs(as.numeric(z))),
    ci_lower  = as.numeric(est - 1.96 * se),
    ci_upper  = as.numeric(est + 1.96 * se))
}


# predictor_significance(): one row of inference for the SELECTED climate
# predictor of a single run, on its best-AICc order. Adds odds-ratio columns
# when the outcome is on the logit scale.

predictor_significance <- function(run) {
  fit  <- refit_best_model(run)
  best <- attr(fit, "best_order")
  w    <- arimax_wald(fit)

  term <- run$best_predictor
  if (!term %in% w$term) term <- "xreg"          # fallback if unnamed
  row  <- w[w$term == term, ]

  base_aic <- run$predictor_results$AIC[
    run$predictor_results$predictor == "none (baseline)"]
  pred_aic <- run$predictor_results$AIC[
    run$predictor_results$predictor == run$best_predictor]

  out <- tibble::tibble(
    state        = run$state,
    clade        = run$clade,
    predictor    = run$best_predictor,
    order        = sprintf("(%d,%d,%d)", best$p, best$d, best$q),
    n_obs        = sum(!is.na(get_outcome(run$data, run$params$scale, run$params$freq))),
    estimate     = round(row$estimate, 4),
    std_error    = round(row$std_error, 4),
    z            = round(row$z, 3),
    p_value      = signif(row$p_value, 3),
    ci_lower     = round(row$ci_lower, 4),
    ci_upper     = round(row$ci_upper, 4),
    # does adding this predictor actually beat the no-climate baseline?
    delta_AIC_vs_baseline = round(pred_aic - base_aic, 2))

  if (identical(run$params$scale, "logit")) {
    out <- dplyr::mutate(out,
      OR       = round(exp(estimate), 3),
      OR_lower = round(exp(ci_lower), 3),
      OR_upper = round(exp(ci_upper), 3))
  }
  out
}


# build_predictor_summary(): stack predictor_significance() over a list of runs
# (e.g. the `runs` from pmap over your state x clade grid). Set adjust = TRUE
# to add a Benjamini-Hochberg FDR column across the whole family of tests.

build_predictor_summary <- function(runs, adjust = TRUE) {
  tbl <- purrr::map_dfr(runs, predictor_significance) %>%
    dplyr::arrange(state, clade)
  if (adjust)
    tbl <- dplyr::mutate(tbl, p_BH = signif(stats::p.adjust(p_value, "BH"), 3))
  tbl
}


# render_significance_table(): shaded gt of the selected-predictor inference,
# grouped by state. Shades the p column (green = small).

render_significance_table <- function(sig_tbl, by_state = TRUE) {
  has_or <- "OR" %in% names(sig_tbl)
  p_col  <- if ("p_BH" %in% names(sig_tbl)) "p_BH" else "p_value"

  keep <- c("clade", "predictor", "order", "n_obs",
            if (has_or) c("OR", "OR_lower", "OR_upper")
            else c("estimate", "ci_lower", "ci_upper"),
            "z", "p_value", if ("p_BH" %in% names(sig_tbl)) "p_BH",
            "delta_AIC_vs_baseline")

  g <- sig_tbl %>%
    dplyr::select(state, dplyr::all_of(keep)) %>%
    { if (by_state) gt::gt(., groupname_col = "state") else gt::gt(.) } %>%
    gt::data_color(
      columns = dplyr::all_of(p_col),
      fn = scales::col_numeric(c("#1a9850", "#ffffbf", "#d73027"),
                               domain = c(0, 1), na.color = "grey85")) %>%
    gt::tab_header(
      title = gt::md("**Selected climate predictor per state x clade**"),
      subtitle = gt::md(paste0(
        "ARIMAX xreg coefficient at the best-AICc order. ",
        "Wald (z) inference; outcome = logit prevalence",
        if (has_or) "; OR = exp(coef)." else "."))) %>%
    gt::tab_source_note(gt::md(
      "_Exploratory: outcome is clade prevalence among publicly deposited "
    )) %>%
    gt::tab_source_note(gt::md(
      "_sequences (convenience sample), and the predictor is chosen by AIC then "
    )) %>%
    gt::tab_source_note(gt::md(
      "_tested on the same series, so p-values are associational only._"))

  if (has_or)
    g <- gt::cols_label(g, OR_lower = "OR 2.5%", OR_upper = "OR 97.5%",
                        delta_AIC_vs_baseline = "dAIC vs base")
  else
    g <- gt::cols_label(g, ci_lower = "2.5%", ci_upper = "97.5%",
                        delta_AIC_vs_baseline = "dAIC vs base")
  g
}


# coef_stability(): refit the SELECTED predictor across every (p,q) cell of the
# sweep and report its coefficient, z and p in each. Answers "is this predictor
# robustly associated, or only at the order AIC happened to pick?" -- a stronger
# robustness statement than a single p-value, and free since the grid exists.

coef_stability <- function(run) {
  p <- run$params; d <- run$differencing
  y  <- get_outcome(run$data, p$scale, p$freq)
  nm <- run$best_predictor
  xr <- matrix(run$data[[nm]], ncol = 1, dimnames = list(NULL, nm))

  run$sensitivity %>%
    dplyr::filter(converged) %>%
    dplyr::select(p, d, q) %>%
    purrr::pmap_dfr(function(p_, d_, q_) {
      fit <- tryCatch(forecast::Arima(
        y, order = c(p_, d_, q_),
        seasonal = list(order = c(p$seasonal_PQ[1], d$D, p$seasonal_PQ[2]),
                        period = p$freq),
        xreg = xr, method = "ML"), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      term <- if (nm %in% names(fit$coef)) nm else "xreg"
      est  <- fit$coef[[term]]
      se   <- sqrt(diag(fit$var.coef))[[term]]
      z    <- est / se
      tibble::tibble(
        state = run$state, clade = run$clade, predictor = nm,
        order = sprintf("(%d,%d,%d)", p_, d_, q_),
        estimate = round(est, 4), z = round(z, 3),
        p_value = signif(2 * stats::pnorm(-abs(z)), 3),
        sig05 = abs(z) > 1.96)
    }) %>%
    dplyr::arrange(dplyr::desc(sig05), p_value)
}


# stability_summary(): collapse coef_stability() to one line per run -- share of
# orders where the predictor is significant and whether the sign ever flips.

stability_summary <- function(runs) {
  purrr::map_dfr(runs, function(r) {
    s <- coef_stability(r)
    tibble::tibble(
      state = r$state, clade = r$clade, predictor = r$best_predictor,
      n_orders      = nrow(s),
      pct_sig_05    = round(100 * mean(s$sig05), 1),
      sign_consistent = length(unique(sign(s$estimate))) == 1L,
      median_estimate = round(stats::median(s$estimate), 4))
  }) %>% dplyr::arrange(state, clade)
}



# reorder_within() / scale_y_reordered(): small tidytext-style helpers so the
# state rows can be ordered by effect size INDEPENDENTLY within each clade
# facet (a state has a different OR in each clade). Dependency-free.

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun, ...)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}


# forest_predictor(): faceted forest plot of the SELECTED predictor's effect,
# one panel per clade, one row per state, colored by which predictor was
# chosen. Uses OR (log x-axis, null = 1) when the table carries OR columns,
# else the raw coefficient (null = 0). Feed it build_predictor_summary().

forest_predictor <- function(sig_tbl, ncol = NULL, point_size = 3.2,
                             state_levels = NULL, base_size = 13,
                             sig_level = 0.10) {
  has_or <- "OR" %in% names(sig_tbl)
  d <- sig_tbl
  if (has_or) {
    d$est <- d$OR; d$lo <- d$OR_lower; d$hi <- d$OR_upper
    null_line <- 1; xlab <- "Odds ratio (95% CI), log scale"; null_lab <- "OR = 1"
  } else {
    d$est <- d$estimate; d$lo <- d$ci_lower; d$hi <- d$ci_upper
    null_line <- 0; xlab <- "xreg coefficient (95% CI)"; null_lab <- "coef = 0"
  }
  # one SHARED, fixed state order across both clade facets so rows line up and
  # you can read straight across the two clades. Pass state_levels to override.
  if (is.null(state_levels)) state_levels <- sort(unique(d$state))
  d$state <- factor(d$state, levels = rev(state_levels))
  # significance from the BH-adjusted p (falls back to raw p if absent).
  # Filled point if p < sig_level, hollow otherwise. NB: keyed to the
  # adjusted p, so a filled point need not have a CI clear of the null.
  pcol  <- if ("p_BH" %in% names(d)) d$p_BH else d$p_value
  p_src <- if ("p_BH" %in% names(d)) "BH p" else "p"
  d$sig <- pcol < sig_level
  lab_t <- sprintf("%s < %.2g", p_src, sig_level)
  lab_f <- sprintf("%s \u2265 %.2g", p_src, sig_level)

  p <- ggplot2::ggplot(d, ggplot2::aes(est, state, color = predictor)) +
    ggplot2::geom_vline(xintercept = null_line, linetype = "dashed",
                        color = "grey50") +
    ggplot2::annotate("text", x = null_line, y = -Inf, label = null_lab,
                      angle = 90, hjust = -0.05, vjust = -0.4,
                      size = base_size / 3.2, color = "grey45") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi),
                           orientation = "y", width = 0.25, linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(shape = sig), size = point_size,
                        stroke = 0.9, fill = "white") +
    ggplot2::facet_wrap(~ clade, ncol = ncol, scales = "free_x") +
    ggplot2::scale_color_manual(values = c("turquoise","cornflowerblue","magenta","darkblue","tomato","gold"),
                                name = "Selected predictor") +
    ggplot2::scale_shape_manual(name = "Significance",
                                values = c(`FALSE` = 21, `TRUE` = 16),
                                labels = c(`FALSE` = lab_f, `TRUE` = lab_t),
                                drop = FALSE) +
    ggplot2::labs(
      x = xlab, y = NULL,
      title = "Selected climate predictor association, by clade",
      subtitle = paste("ARIMAX xreg at best-AICc order; one row per state.",
                       "Exploratory \u2014 prevalence among sampled sequences.")) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   legend.position  = "bottom",
                   axis.text  = ggplot2::element_text(size = base_size + 2),
                   axis.title = ggplot2::element_text(size = base_size + 3),
                   strip.text = ggplot2::element_text(size = base_size + 2)) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 16)))
  if (has_or) p <- p + ggplot2::scale_x_log10()
  p
}


## part five: maps ----

make_diffusion_map <- function(tree_file,
                               subclade_name,
                               lat_col          = "location1",
                               lon_col          = "location2",
                               pad              = 4.0,
                               bw_mult_lon      = 1.2,
                               bw_mult_lat      = 1.5,
                               kde_n            = 400,
                               density_quantile = 0.50,
                               map_xlim         = c(-105, -72),
                               map_ylim         = c(24, 50),
                               states_data      = states,
                               centers_data     = state_centers,
                               save_path        = NULL,
                               save_width       = 9,
                               save_height      = 7,
                               save_dpi         = 300) {
  
  ## Read tree and pull node table
  mcc_tree  <- read.beast(tree_file)
  node_data <- as_tibble(mcc_tree)
  
  ## Extract node coordinates
  coords <- node_data %>%
    filter(!is.na(.data[[lat_col]]), !is.na(.data[[lon_col]])) %>%
    mutate(lat = as.numeric(.data[[lat_col]]),
           lon = as.numeric(.data[[lon_col]])) %>%
    dplyr::select(node, label, lat, lon, height)
  
  ## KDE bounds (data extent + padding)
  lon_range <- c(min(coords$lon) - pad, max(coords$lon) + pad)
  lat_range <- c(min(coords$lat) - pad, max(coords$lat) + pad)
  
  ## 2D kernel density estimate
  kde <- MASS::kde2d(
    x    = coords$lon,
    y    = coords$lat,
    n    = kde_n,
    lims = c(lon_range, lat_range),
    h    = c(MASS::bandwidth.nrd(coords$lon) * bw_mult_lon,
             MASS::bandwidth.nrd(coords$lat) * bw_mult_lat)
  )
  
  ## Convert KDE to data frame, keep top cells only
  kde_df <- expand.grid(lon = kde$x, lat = kde$y) %>%
    mutate(density = as.vector(kde$z)) %>%
    filter(density > quantile(density, density_quantile))
  
  ## Build plot
  p <- ggplot() +
    geom_sf(data = states_data,
            fill = "#f0ede8", color = "#d0ccc8", linewidth = 0.3) +
    geom_tile(data = kde_df,
              aes(x = lon, y = lat, fill = density, alpha = density)) +
    scale_fill_scico(palette = "hawaii", direction = -1, name = "Density") +
    scale_alpha_continuous(range = c(0.2, 0.9), guide = "none") +
    geom_sf(data = centers_data,
            color = "pink4", fill = "gold", size = 3, shape = 22) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    theme_classic(base_size = 16) +
    theme(legend.position    = "right",
          legend.key.height  = unit(1.2, "cm"),
          plot.title         = element_text(face = "bold", size = 18),
          axis.line          = element_blank(),
          axis.ticks         = element_blank(),
          axis.text          = element_blank()) +
    labs(subtitle = sprintf("KDE of %s locations from subclade MCC tree",
                            subclade_name), x = "", y = "")
  
  ## Optional save
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p,
           width = save_width, height = save_height,
           dpi = save_dpi, bg = "white")
  }
  
  return(p)
}

