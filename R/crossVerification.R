

crossVerification <- function(Mf,Mm,Mp,
                              Mexp=NULL,
                              sc_filter=NULL,
                              ploidy=2,
                              het = NULL){

  # take the markers of the progeny (Mp) and
  # calculate the probability of belonging to a given male and female
  if(any(is.null(nrow(Mf)) | is.null(nrow(Mm)) | is.null(nrow(Mp)) )){
    stop("One or more of the provided arguments is not a matrix", call. = FALSE)
  }
  if(is.null(Mexp)){ Mexp <- (Mf + Mm)/2 } # expected genotype

  score_Mf = (Mf==0|Mf==ploidy)
  score_Mm = (Mm==0|Mm==ploidy)
  score_Mexp = score_Mm + score_Mf

  #Calculate match probabilities

  #midpoint
  m = ploidy %/% 2L

  # Initialize result object
  resultMatch = Matrix::Matrix(0, nrow = nrow(Mp), ncol = ncol(Mp),
                        dimnames = dimnames(Mp), sparse = T)

  # Only score cells where all three are present
  ok = !is.na(Mf) & !is.na(Mm) & !is.na(Mp)

  #Make NA if there's missing genotype data
  resultMatch[!ok] <- NA_real_

  ## homo × homo -----

  # 0 × 0  -> offspring must be 0
  idx = ok & (Mf == 0L & Mm == 0L)
  if (any(idx)) resultMatch[idx] <- as.numeric(Mp[idx] == 0L)

  # P × P  -> offspring must be P
  idx = ok & (Mf == ploidy & Mm == ploidy)
  if (any(idx)) resultMatch[idx] = as.numeric(Mp[idx] == ploidy)

  # 0 × P or P × 0 -> offspring must be m
  idx = ok & ((Mf == 0L & Mm == ploidy) | (Mf == ploidy & Mm == 0L))
  if (any(idx)) resultMatch[idx] <- as.numeric(Mp[idx] == m)

  ## ----- Heterozygote × homozygote (either side), no double reduction -----
  ## Gamete from het draws m alleles without replacement (hypergeometric).
  ## Offspring dosage d_off = h + k, where:
  ##   h = m if homo parent is all A (dosage P), else 0 if homo is all a (dosage 0)
  ##   k ~ Hypergeom(K = i_het, N = P, draws = m), so
  ##   P(D = d_off) = [ C(i_het, k) * C(P - i_het, m - k) ] / C(P, m),  where k = d_off - h

  idx_het_hom = ok & score_Mexp == 1

  if (any(idx_het_hom)) {
    lin = which(idx_het_hom)  # indices of relevant cells

    # Identify which side is the heterozygote for each cell
    mom_is_het = (Mf[lin] > 0L & Mf[lin] < ploidy) & (Mm[lin] %in% c(0L, ploidy))

    i_het   = ifelse(mom_is_het, Mf[lin], Mm[lin])  # A dosage in heterozygous parent
    hom_dos = ifelse(mom_is_het, Mm[lin], Mf[lin])  # 0 or ploidy
    d_off   = Mp[lin]                                # offspring A dosage

    # h = contribution of homozygous parent to A in the gamete pair (m if all A, else 0)
    h = ifelse(hom_dos == ploidy, m, 0L)
    k = d_off - h

    # Valid parameter region for hypergeometric combo counts
    okk = (k >= 0L & k <= m) & (i_het > 0L & i_het < ploidy) & (hom_dos %in% c(0L, ploidy))

    p = numeric(length(lin))
    if (any(okk)) {
      denom = choose(ploidy, m)
      p[okk] = choose(i_het[okk], k[okk]) *
        choose(ploidy - i_het[okk], m - k[okk]) / denom
    }
    p[!is.finite(p)] = 0
    resultMatch[lin] = p
  }

  ## Heterozygote × heterozygote (double het)
  ## Gamete dosages are hypergeometric (0..m),
  ## offspring dosage distribution is the convolution of parental gamete distributions.

  idx_het_het = ok & score_Mexp == 0

  if (any(idx_het_het)) {
    lin  = which(idx_het_het)
    keys = paste(Mf[lin], Mm[lin], sep=":")
    ukeys = unique(keys)

    for (key in ukeys) {
      pos = lin[keys == key]
      ab  = strsplit(key, ":", fixed = TRUE)[[1]]
      a   = as.integer(ab[1])  # maternal dosage (heterozygote)
      b   = as.integer(ab[2])  # paternal dosage (heterozygote)

      # Gamete distributions for allele 'A' from each parent
      x  = 0:m
      px = dhyper(x, m = a, n = ploidy - a, k = m)  # mother
      py = dhyper(x, m = b, n = ploidy - b, k = m)  # father

      # Convolution to get offspring dosage probabilities pD
      pD = numeric(ploidy + 1L)
      for (d in 0:ploidy) {
        x_min = max(0L, d - m)
        x_max = min(m,   d)
        if (x_max >= x_min) {
          xs = x_min:x_max
          pD[d + 1L] = sum(px[xs + 1L] * py[(d - xs) + 1L])
        }
      }
      s = sum(pD)
      if (s > 0) pD = pD / s  # numerical tidy

      # Support bounds (everything outside gets probability 0)
      D_low  <- max(0L, a - m) + max(0L, b - m)
      D_high <- min(m,  a)     + min(m,  b)

      # Observed dosage
      lut = numeric(ploidy + 1L)
      if (D_high >= D_low) lut[(D_low:D_high) + 1L] = pD[(D_low:D_high) + 1L]

      # Assign for all cells sharing this (a,b) pair
      d_obs = Mp[pos]
      okd   = which(!is.na(d_obs) & d_obs >= 0L & d_obs <= ploidy)
      if (length(okd)) resultMatch[pos[okd]] = lut[d_obs[okd] + 1L]

      # Keep NA where offspring dosage is missing
      ina = which(is.na(d_obs))
      if (length(ina)) resultMatch[pos[ina]] = NA_real_
    }
  }

  # compute metrics for individuals
  # Minimum number of non-missing markers required to compute probMatch
  min_markers <- 8L

  if(is.null(sc_filter)){
    probMatch <- apply(resultMatch, 1, function(x) {
      den <- sum(!is.na(x))
      if (den >= min_markers) sum(x, na.rm = TRUE) / den else NA_real_
    })
  }else if(sc_filter == "Score2"){
    sel = (score_Mexp == 2L)
    probMatch <- apply(ifelse(sel, resultMatch, NA_real_), 1, function(x) {
      den <- sum(!is.na(x))
      if (den >= min_markers) sum(x, na.rm = TRUE) / den else NA_real_
    })
  }else if(sc_filter == "ScoreNon0"){
    sel = (score_Mexp != 0)
    probMatch <- apply(ifelse(sel, resultMatch, NA_real_), 1, function(x) {
      den <- sum(!is.na(x))
      if (den >= min_markers) sum(x, na.rm = TRUE) / den else NA_real_
    })
  }

  #Observed heterozigosity
  heteroMp <- {
    num <- rowSums(Mp == 1, na.rm = TRUE)
    den <- rowSums(!is.na(Mp))
    ifelse(den > 0, num / den, NA_real_)
  } # heterozigosity found in progeny

  #Expected heterozigosity
  # Precompute hypergeometric: gamete Pr(0) and Pr(m) for all parental dosages
  px0 = sapply(0:ploidy, function(a) dhyper(0, m=a, n=ploidy-a, k=m))
  pxm = sapply(0:ploidy, function(a) dhyper(m, m=a, n=ploidy-a, k=m))


  phet <- matrix(NA_real_, nrow(Mexp), ncol(Mexp))

  # Case masks (all N×M)
  is_hom0 = (Mf == 0L);     is_homP = (Mf == ploidy)
  js_hom0 = (Mm == 0L);     js_homP = (Mm == ploidy)
  is_hetM = (Mf > 0L & Mf < ploidy)
  is_hetF = (Mm > 0L & Mm < ploidy)

  # homo×homo (same allele) -> always homozygous
  mask_hh_same = (is_hom0 & js_hom0) | (is_homP & js_homP)
  phet[mask_hh_same] = 0

  # 0×P or P×0 -> always heterozygous
  mask_hh_oppo = (is_hom0 & js_homP) | (is_homP & js_hom0)
  phet[mask_hh_oppo] = 1

  # het×hom (either side)
  mask_het_hom <- (is_hetM & (js_hom0 | js_homP)) | (is_hetF & (is_hom0 | is_homP))
  if (any(mask_het_hom)) {
    # pick the heterozygote dosage and the homozygote dosage for each cell
    i_het   <- ifelse(is_hetM & (js_hom0 | js_homP), Mf, Mm)
    hom_dos <- ifelse(is_hetM & (js_hom0 | js_homP), Mm, Mf)
    # phomo
    idx <- which(mask_het_hom)
    i_vec <- i_het[idx] + 1L
    h_vec <- hom_dos[idx]
    phomo <- ifelse(h_vec == 0L, px0[i_vec], pxm[i_vec])
    phet[idx] <- 1 - phomo
  }

  # het×het
  mask_het_het <- (is_hetM & is_hetF)
  if (any(mask_het_het)) {
    idx <- which(mask_het_het)
    a <- Mf[idx] + 1L
    b <- Mm[idx] + 1L
    # phet
    phet[idx] <- 1 - (px0[a] * px0[b] + pxm[a] * pxm[b])
  }

  # Per individual
  mask_use <- if (exists("Mexp")) !is.na(Mexp) else !is.na(phet)
  num <- rowSums(phet * mask_use, na.rm = TRUE)
  den <- rowSums(mask_use & !is.na(phet))
  heteroMexp <- ifelse(den > 0, num / den, NA_real_)

  #Het deviation
  heteroDeviation <- abs(heteroMexp - heteroMp) # deviation from expected heterosigosity

  #Average information score
  if(is.null(sc_filter)){
    avgScore = rowMeans(score_Mexp, na.rm = T)
  }else if(sc_filter == "Score2"){
    sel = (score_Mexp == 2L)
    avgScore = rowMeans(ifelse(sel, score_Mexp, NA_real_), na.rm = TRUE)
  }else if(sc_filter == "ScoreNon0"){
    sel = (score_Mexp != 0)
    avgScore = rowMeans(ifelse(sel, score_Mexp, NA_real_), na.rm = TRUE)
  }

  #Number of highly informative markers
  if(is.null(sc_filter)){
    nScore2 = rowSums((score_Mexp == 2), na.rm = T)
  }else if(sc_filter == "Score2"){
    sel = (score_Mexp == 2L)
    nScore2 = rowSums(ifelse(sel, (score_Mexp == 2), NA_real_), na.rm = TRUE)
  }else if(sc_filter == "ScoreNon0"){
    sel = (score_Mexp != 0)
    nScore2 = rowSums(ifelse(sel, (score_Mexp == 2), NA_real_), na.rm = TRUE)
  }

  #Total number of markers
  if(is.null(sc_filter)){
    nMarkers = rep(ncol(score_Mexp),nrow(score_Mexp))
  }else if(sc_filter == "Score2"){
    nMarkers = rowSums(score_Mexp == 2L, na.rm = T)
  }else if(sc_filter == "ScoreNon0"){
    nMarkers = rowSums(score_Mexp != 0, na.rm = T)
  }

  #Parental heterozygosity filter
  homo_Mf <- (Mf == 0 | Mf == ploidy)
  n_non_missing <- rowSums(!is.na(Mf))
  n_homo <- rowSums(homo_Mf, na.rm = TRUE)
  het_pct_fem <- (pmax(n_non_missing - n_homo, 0) / pmax(n_non_missing, 1)) * 100

  homo_Mm <- (Mm == 0 | Mm == ploidy)
  n_non_missing <- rowSums(!is.na(Mm))
  n_homo <- rowSums(homo_Mm, na.rm = TRUE)
  het_pct_mal <- (pmax(n_non_missing - n_homo, 0) / pmax(n_non_missing, 1)) * 100

  parHetFilter = as.numeric(het_pct_fem <= het & het_pct_mal <= het)

  res <- data.frame(designation=rownames(Mp),probMatch,heteroMp,heteroMexp,heteroDeviation,avgScore,nScore2,nMarkers,parHetFilter)
  rownames(res) <- NULL
  # compute metrics for markers
  return(list(metricsInd=res,matchMat=resultMatch,Mprogeny=Mp, Mfemale=Mf, Mmale=Mm, Mexpected=Mexp))

}

