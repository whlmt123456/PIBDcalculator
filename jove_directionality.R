args <- commandArgs(TRUE)
if (length(args) != 6) stop("Usage: Rscript jove_directionality.R archive_dir extracts_dir manifest.tsv plink eur_bfile output_dir")
arc <- args[1]; ext <- args[2]; man <- args[3]
plink <- normalizePath(args[4], mustWork=TRUE); ref <- args[5]; out <- args[6]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(TwoSampleMR))
endo <- c(APPENDACUT="K11_APPENDACUT",APPENDACUT_COMPLIC="APPENDACUT_COMPLIC",APPENDOTH="K11_APPENDOTH")
traits <- c(saturated="ebi-a-GCST90092980",Monounsaturated="ebi-a-GCST90092928",Polyunsaturated="ebi-a-GCST90092939",
            OMEGA6="ebi-a-GCST90092933",OMEGA3="ebi-a-GCST90092931",Polyradio="ebi-a-GCST90092941",
            Monoradio="met-d-MUFA_pct",OMEGAradio="ebi-a-GCST90092934")
manifest <- read.delim(man)
manifest <- manifest[match(endo, manifest$phenocode),]
stopifnot(identical(manifest$phenocode,unname(endo)),all(manifest$num_cases>0),all(manifest$num_controls>0))
loadtsv <- function(n) read.delim(file.path(ext,paste0(n,".tsv")), colClasses=c(chrom="character"))
readobj <- function(n) readRDS(file.path(arc,paste0(n,".rds")))
save <- function(df,n) write.csv(df,file.path(out,paste0(n,".csv")),row.names=FALSE,na="")
harm <- readobj("all_harmonised_dat_list")
look <- setNames(lapply(names(traits),function(t) loadtsv(paste0(t,"_lookup"))),names(traits))
save(manifest,"FinnGen_R11_endpoints")

aligned <- function(s,ea,oa,d) {
  q <- d[d$SNP==s,,drop=FALSE]
  idx <- (q$ea==ea & q$oa==oa) | (q$ea==oa & q$oa==ea)
  q <- q[which(idx),,drop=FALSE]
  if (nrow(q)!=1L) return(c(beta=NA,se=NA))
  c(beta=if(q$ea==ea)q$beta else -q$beta,se=q$se)
}
validate <- function(d,q,side) {
  raw <- t(vapply(seq_len(nrow(d)),function(i)
    aligned(d$SNP[i],d[[paste0("effect_allele.",side)]][i],d[[paste0("other_allele.",side)]][i],q),
    c(beta=0,se=0)))
  ok <- rowSums(is.na(raw))==0
  mb <- abs(raw[,"beta"] - d[[paste0("beta.",side)]])
  ms <- abs(raw[,"se"] - d[[paste0("se.",side)]])
  c(matched=sum(ok),agree=sum(ok & mb<2e-6 & ms<2e-6,na.rm=TRUE))
}
checks <- do.call(rbind,lapply(names(harm),function(k){
  a <- strsplit(k,"__",fixed=TRUE)[[1]]; t<-a[1]; e<-a[2]
  d <- harm[[k]][which(harm[[k]]$mr_keep),,drop=FALSE]
  x <- validate(d,look[[t]],"exposure"); y <- validate(d,loadtsv(paste0(e,"_forward")),"outcome")
  data.frame(trait=t,endpoint=endo[e],n=nrow(d),exposure_match=x["matched"],exposure_agreement=x["agree"],
             outcome_match=y["matched"],outcome_agreement=y["agree"])
}))
save(checks,"Source_concordance")
print(checks)
if (any(checks$exposure_agreement != checks$n | checks$outcome_agreement != checks$n))
  stop("The original forward-variant lookups did not fully match the archive. Review Source_concordance.csv before interpreting new analyses.")

mufa_n <- 114999L
prevalences <- c(.001,.005,.01,.05,.10)
steig <- do.call(rbind,lapply(names(harm),function(k){
  pair <- strsplit(k,"__",fixed=TRUE)[[1]]; t <- pair[1]; e <- pair[2]; j <- match(e,names(endo))
  d <- harm[[k]][which(harm[[k]]$mr_keep),,drop=FALSE]
  n_exp <- if(t=="Monoradio") mufa_n else unique(d$samplesize.exposure)
  n_total <- manifest$num_cases[j]+manifest$num_controls[j]
  n_eff <- TwoSampleMR::effective_n(manifest$num_cases[j],manifest$num_controls[j])
  r_exp <- if(t=="Monoradio"){
    z <- d$beta.exposure/d$se.exposure
    if(is.finite(n_exp)) z/sqrt(z^2+n_exp-2) else rep(NA_real_,nrow(d))
  } else d$beta.exposure*sqrt(2*d$eaf.exposure*(1-d$eaf.exposure))
  do.call(rbind,lapply(prevalences,function(kpop){
    r_out <- TwoSampleMR::get_r_from_lor(d$beta.outcome,d$eaf.outcome,
      manifest$num_cases[j],manifest$num_controls[j],kpop)
    usable <- all(is.finite(r_exp)) && all(is.finite(r_out))
    rx2 <- if(usable) sum(r_exp^2) else NA_real_; ry2 <- sum(r_out^2)
    if (usable && (rx2>=1 || ry2>=1)) stop("Invalid aggregate R-squared")
    z <- if(usable) (atanh(sqrt(rx2))-atanh(sqrt(ry2)))/sqrt(1/(n_exp-3)+1/(n_eff-3)) else NA_real_
    logp <- log(2)+pnorm(-abs(z),log.p=TRUE)
    data.frame(trait=t,exposure_id=unname(traits[t]),endpoint=unname(endo[e]),n_snp=nrow(d),
      n_exposure=n_exp,n_cases=manifest$num_cases[j],n_controls=manifest$num_controls[j],n_outcome_total=n_total,
      n_outcome_effective=n_eff,assumed_population_prevalence=kpop,r2_exposure=rx2,r2_outcome=ry2,
      direction=if(!usable) "Not estimable: met-d sample size unverified" else if(rx2>ry2) "Fatty acid to appendicitis" else "Appendicitis to fatty acid",
      p_steiger=if(is.finite(logp) && logp>=log(1e-300)) exp(logp) else NA_real_,
      log10_p_steiger=logp/log(10),snps_supporting_forward=if(usable)sum(r_exp^2>r_out^2)else NA_integer_)
  }))
}))
save(steig,"Steiger_directionality")

rev <- list(); qc <- list(); rev_snps <- list()
for(e in names(endo)){
  all_c <- loadtsv(paste0(e,"_gws"))
  stopifnot(all(all_c$p<5e-8),all(all_c$se>0))
  cands <- all_c[!duplicated(all_c$SNP) & !duplicated(all_c$SNP,fromLast=TRUE),]
  pref <- file.path(out,paste0(e,"_reverse_ld"))
  if(nrow(cands)>0) {
    write.table(data.frame(SNP=cands$SNP,P=cands$p),paste0(pref,".tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    plargs <- c("--bfile",shQuote(ref),"--clump",shQuote(paste0(pref,".tsv")),"--clump-p1","5e-8",
      "--clump-p2","5e-8","--clump-kb","10000","--clump-r2","0.001","--threads","2","--out",shQuote(pref))
    ret <- system2(plink,plargs,stdout=TRUE,stderr=TRUE)
    if (!is.null(attr(ret,"status"))) stop(paste(ret,collapse="\n"))
    lead <- if(file.exists(paste0(pref,".clumped"))) read.table(paste0(pref,".clumped"),header=TRUE)$SNP else character()
  } else lead <- character()
  cand <- cands[match(lead,cands$SNP),,drop=FALSE]
  cand$instrument_f <- (cand$beta/cand$se)^2
  qc[[e]] <- data.frame(endpoint=endo[e],gws_biallelic_records=nrow(all_c),gws_unique_rsid=nrow(cands),clumped=nrow(cand))
  for(t in names(traits)){
    d<-look[[t]]
    d<-d[!duplicated(d$SNP)&!duplicated(d$SNP,fromLast=TRUE),,drop=FALSE]
    b <- t(vapply(seq_len(nrow(cand)),function(i)
      aligned(cand$SNP[i],cand$ea[i],cand$oa[i],d),c(beta=0,se=0)))
    ambig <- paste0(cand$ea,cand$oa) %in% c("AT","TA","CG","GC")
    keep <- !ambig & complete.cases(b)
    h <- cand[keep,,drop=FALSE]
    n <- nrow(h)
    status <- if(n>0L) "Estimated" else "No harmonized instruments"
    method <- if(n>=2L) "Inverse variance weighted" else if(n==1L) "Wald ratio" else "Not estimable"
    f <- if(n>=2L) TwoSampleMR::mr_ivw else TwoSampleMR::mr_wald_ratio
    fit <- if(n) f(h$beta,b[keep,"beta"],h$se,b[keep,"se"],TwoSampleMR::default_parameters()) else list(b=NA_real_,se=NA_real_,pval=NA_real_)
    rr <- data.frame(endpoint=unname(endo[e]),trait=t,fatty_acid_id=unname(traits[t]),method=method,status=status,clumped_snp=nrow(cand),
      palindromic_removed=sum(ambig),missing_or_incompatible=sum(!ambig & !complete.cases(b)),harmonized_snp=n,
      min_f=if(n)min(h$instrument_f)else NA_real_,beta_per_log_odds=fit$b, se=fit$se,
      ci_lower=fit$b-1.96*fit$se,ci_upper=fit$b+1.96*fit$se,p=fit$pval,Q_p=if(n>=2L)fit$Q_pval else NA_real_,
      snps=paste(h$SNP,collapse=";"))
    rev[[length(rev)+1]] <- rr
    if(n) rev_snps[[length(rev_snps)+1]] <- data.frame(endpoint=unname(endo[e]),trait=t,SNP=h$SNP,ea=h$ea,oa=h$oa,
       beta_appendicitis=h$beta,se_appendicitis=h$se,p_appendicitis=h$p,
       beta_fatty_acid=b[keep,"beta"],se_fatty_acid=b[keep,"se"],instrument_f=h$instrument_f)
  }
}
reverse <- do.call(rbind,rev)
reverse$p_bonferroni <- p.adjust(reverse$p,method="bonferroni",n=24)
save(reverse,"Reverse_MR")
save(do.call(rbind,qc),"Reverse_instrument_QC")
if(length(rev_snps)) save(do.call(rbind,rev_snps),"Reverse_instruments")
writeLines(c(capture.output(sessionInfo()),"",system2(plink,"--version",stdout=TRUE)),file.path(out,"Directionality_sessionInfo.txt"))
