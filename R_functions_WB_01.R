library(ape)
library(seqinr)
library(tidyverse)



#takes as input list returned by bin_gen_MHCBD
#calculates diversity using hillR and taking into account
#various special cases
#IMPORTANT: Xes are handled naturally by DNAbin
#IMPORTANT: in  a situation when S = 1 and > 2alleles (3 variants or 2 variants a N)
#distances between all variants are set to 1 in such cases - I think it's acceptable
#returns list with two dataframes:
#"alpha_gamma" alpha and gamma diversities for species
#"individual_alpha" alpha diversities for each individual
#both df contain xtensive metainfo

#takes as input list returned by bin_gen
#i.e. each segment or collapsed segment
#calculates diversity using hillR and taking into account
#various special cases
#IMPORTANT: Xes are handled naturally by DNAbin
#IMPORTANT: requires fixing a situation when S = 1 and > 2alleles (3 variants or 2 variants a N)
#currently distances between all variants are set to 1 in such cases
#get genotypes and strip individual names


calc_div <- function(bg, sp, typ = c("DNA", "codon", "AA", "AAGhm"), Ghm_table){
  #takes vector returned by hill_phylo (used to get individual diversities)
  #and converts it into dataframe
  div_vect_to_df <- function(l, ind) {
    ldf <- vector("list", 3)
    for (i in seq_along(l)) {
      q = i-1
      ldf[[i]] <- data.frame(id = ind, q = q, alpha = l[[i]])
    }
    return(bind_rows(ldf))
  }
  #get genotypes and strip individual names
  g <- bg[["genotypes"]] %>% select(-c(1))
  #get individual ids
  ids <- bg[["genotypes"]] %>% pull(id)
  #get allele names
  a_nam <- colnames(g)
  #summary
  s <- bg[["summary"]]
  #allele sequences as DNA/AAbin
  al_seq <- bg[["seq"]]
  #sets all diversity values in output to zero
  #for species estimates of alpha and gamma diversity
  #called later when needed
  zero_div <- function(){
    q <- c(0,1,2)
    PD_gamma <- rep(0,3)
    PD_alpha  <- rep(0,3)
    PD_beta <- rep(0,3)
    local_similarity <- rep(0,3)
    region_similarity <- rep(0,3)
    res <- data.frame(q, PD_gamma, PD_alpha, PD_beta, local_similarity, region_similarity)
    return(res)
  }
  #as above but for individual alpha diversities
  #below things with suffix _a are for individual alpha diversities
  zero_div_a <- function(ind) {
    ldf <- vector("list", 3)
    for (i in c(1:3)) {
      ldf[[i]] <- data.frame(id = ind, q = (i-1), alpha = 0)
    }
    return(bind_rows(ldf))
  }
  #if only one variant, sets diversity values into 0
  #otherwise converts df with allele sequences into DNAbin or AAbin
  if(length(a_nam) == 1){
    dsum <- zero_div()
    dsum_a <- zero_div_a(ids)
    S <- 0
    S_unamb <- 0
    dmax <- 0
  } else {
    #calculates pairwise distances between alleles using an appropriate distnce:
    #DNA: TN93 
    #codon: dS (Li 93)
    #AA: raw amino-acid p-distance
    #AAGhm: Grantham distance /100
    S <- length(seg.sites(al_seq, strict = TRUE, trailingGapsAsN = FALSE))
    S_unamb <- length(seg.sites(al_seq, strict = FALSE, trailingGapsAsN = FALSE))
    if(typ == "DNA"){
      #we may want to consider for MHC something correcting for multiple hits, like TN93
      seq_div <- dist.dna(al_seq, model = "raw", pairwise.deletion = TRUE)
    } else if(typ == "codon"){
      #kaks returns all values 10 if considers the alignment saturated
      #it'be best to get another package to do the calculations
      #but couldn't find, so the fix is described below
      seq_div <- kaks(ape::as.alignment(as.matrix(al_seq)), rmgap = FALSE)$ks
    } else if(typ == "AA" || typ == "AAGhm"){
      S <- NA
      S_unamb <- length(AAsubst(as.matrix(al_seq)))
      if (typ == "AA"){
        seq_div <- dist.aa(as.matrix(al_seq), pairwise.deletion = TRUE, scaled = TRUE)  
      } else {
        seq_div = calculate_grantham(al_seq, Ghm_table)
      }
    }
    dmax <- max(seq_div, na.rm = TRUE)
    #deals with situations when DNAdiv can't be estimated (too high divergence)
    #e.g., segment of length 1 with 1 variable site
    #currently it's quite arbitrary, but the problem appears only in dS estimation
    #for very short segments
    seq_div[is.nan(seq_div) | seq_div > 3] <- 1
    seq_div[seq_div < 0] <- 0
    if(sum(seq_div) == 0){
      #alleles differ only by missing data
      dsum_a <- zero_div_a(ids)
      dsum <- zero_div()
    } else {
      if(length(a_nam) == 2){
        #if only two alleles, phylo object is created manually 
        #this is just for consistency, so that we have tree for two alleles as well
        bl <- seq_div[1]/2
        tr <- list(edge = matrix(c(3, 3, 1, 2), 2, 2), 
                   edge.length = rep(bl, 2), 
                   tip.label = a_nam, 
                   Nnode = 1L)
        class(tr) <- "phylo"
      }
      else {
        #constructs tree from distance matrix
        tr <- bionj(seq_div)
        #sets negative branch lengths to 0
        tr$edge.length[tr$edge.length<0] <- 0
      }
      #calculates diversities for various q values
      qval <- c(0:2)
      #species alpha and gamma diversities
      #rel_ten_pool should be TRUE, because then we give equal wiegths to assemblages=individuals
      div <- map(qval, hill_phylo_parti, comm = g, tree = tr, rel_then_pool = TRUE)
      dsum <- bind_rows(div)
      #alpha diversities for each individual
      div_a <- map(qval, hill_phylo, comm = g, tree = tr, rel_then_pool = TRUE)
      dsum_a <- div_vect_to_df(div_a, ids)
    }
  }
  Ss <- data.frame(S = S, S_unamb = S_unamb, dmax = dmax)
  #combines summary with diversity values
  res_alpha_gamma <- bind_cols(s, Ss, dsum)
  print(res_alpha_gamma)
  res_ind_alpha <- bind_cols(s, Ss, dsum_a)
  print(res_ind_alpha)
  res <- list("alpha_gamma" = res_alpha_gamma, "individual_alpha" = res_ind_alpha)
  return(res)
}

#converts DNAbin or AAbin into dataframe
bin2df <- function(bin){
  temp <- as.character(bin)
  res <- as.data.frame(sapply(temp, paste, collapse = ''))
  colnames(res) <- "seq"
  res$label <- rownames(res)
  rownames(res) <- c()
  res <- res %>% select(label, seq)
  return(res)
}

#calculates_Grantham distance for AAbin object (aligned)
#Written by Tomek Gaczorek, modified by WB
#now the user has to provide the Grantham table - it's much faster
calculate_grantham <- function(protein_allignment, gram_table){
  # input - object of class AAbin
  # grantham table
  #gram_table <- read.csv("https://dl.dropbox.com/s/0awueqeyodsvveu/Gratham_Distance_table.csv?dl=1",
  #                       header = T, row.names = 1,encoding = "UTF-8")
  #outcome matrix
  out_mat <- matrix(0,nrow = length(protein_allignment),ncol = length(protein_allignment))
  pa_names <- names(protein_allignment)
  colnames(out_mat) <- pa_names
  row.names(out_mat) <- pa_names
  
  # distance between 2 sequences
  gram_dist <- function(seq1,seq2,gr_tab){
    out_dist <- c()
    for(k in 1:length(seq1)){
      one <- seq1[k]
      two <- seq2[k]
      if((one == "X" || two == "X")){
        out_dist <- c(out_dist,NA)
      }
      else if(one == "*" || two == "*"){
        next
      }
      else if(identical(one,two) == T){
        out_dist <- c(out_dist,0)
      }
      else {
        if(two == "S" || is.na(gr_tab[one,two])==T){
          out_dist <- c(out_dist,gr_tab[two,one]) 
        } else {out_dist <- c(out_dist,gr_tab[one,two])}
      }
    }
    out_dist <- na.omit(out_dist)
    sum(out_dist)/length(out_dist)
  }
  
  # calculation of distance for each pair
  for(i in c(1:length(protein_allignment))){
    for(j in c(1:length(protein_allignment))){
      if(pa_names[i] != pa_names[j]){
        out_mat[i,j] <- gram_dist(toupper(as.character(protein_allignment[[i]])),
                                  toupper(as.character(protein_allignment[[j]])),
                                  gram_table)
      }
    }
  }
  #divide values by 100 to avoid very large distances in distance matrix - 
  #WB 11.12.20
  out_mat <- out_mat/100
  as.dist(out_mat)
}


#df to DNAbin object
#expects  names in 1st and sequences in 2nd, ignores remaining columns
df2DNA <- function(df){
  DNA <- t(lapply(strsplit(df[ , 2], ""), tolower))
  if(dim(DNA)[1] == 1){
    dim(DNA) <- c(length(DNA), 1)
  }
  names(DNA) <- df[, 1]
  DNA <- as.DNAbin(DNA)
  return(DNA)
}


#df to DNAbin object
#expects  names in 1st and sequences in 2nd, ignores remaining columns
df2prot <- function(df){
  prot <- t(lapply(strsplit(df[ , 2], ""), tolower))
  if(dim(prot)[1] == 1){
    dim(prot) <- c(length(prot), 1)
  }
  names(prot) <- df[, 1]
  prot <- as.AAbin(prot)
  return(prot)
}


#countalleles from a category, which are present in at least one ind in inp
count_alleles <- function(inp, cat){
  count <- inp %>% select(one_of(cat)) %>%
    select_if(~ !is.integer(.) || sum(., na.rm = T ) > 0) %>% length()
  return(count)
}

#
#gets vector of allele names that are present in at least min_count individuals in inp
get_all_names <- function(inp, names, min_count = 1){
  all_to_keep <- inp %>% select(one_of(names)) %>% select_if(~ sum(., na.rm = T) >= min_count)
  all_names <- colnames(all_to_keep)
  return(all_names)
}

get_alleles <- function(inp, names, out_name, min_count = 1, to_trim_start = 0, to_trim_end = 0){
  all_names <- get_all_names(inp, names, min_count)
  all <- alleles %>% filter(allele %in% all_names) %>%
    transmute(seq.name = paste(allele),
              #seq.name = paste(allele, class, category, sep = "_"),
              seq.text = str_sub(sequence, start = to_trim_start + 1, end = -(to_trim_end+1)))
  dat2fasta(all, outfile = out_name)
  #return(all)
}

#takes output from gen_gen and individual info object - dataframe created from tab-delimited file like this:
#individual_id	family	genus	species	subspecies	Locality	Country
#14829	Ambystomatidae	Ambystoma	texanum		Waco	USA
#14830	Ambystomatidae	Ambystoma	texanum		Waco	USA
#produces dataframe with additional columns annotating individuals and counts converted into allele presence-absence
get_binary_annotated <- function(genotypes, ind_info){
  #transform read counts into 1/0 (presence/absence), keep NA (although there shouldn't be any)
  #skip variables other than alleles when transforming to binary
  d <- genotypes
  to_skip <-c(grep('cov_ex[2-3]$', colnames(d), value = T), 
              grep(paste0("cov_ex[2-3]_alleles"), colnames(d), value = T),
              "individual_id")
  #print(to_skip)
  to_keep <- setdiff(names(d), to_skip)
  d[,to_keep][d[,to_keep] >0 ] <- 1
  #join genotypes with individual info
  with_tax <- left_join(d, select(ind_info, -c("subspecies", "Country")), by = "individual_id")  
  #reorder column to get all annotations at the beginning
  refcol <- c("individual_id", "family", "genus", "species", "Locality", to_skip[1:2]) 
  with_tax <- with_tax %>% select(refcol, everything()) 
}


#reads genotypes in the tab-delimited format from multiple files containing genera in names:
#col1 - coverage
#col2 - coverage of alleles
#col3 - individual_id
#col4: - read counts for allelea (their names in col)
#cov_ex2	cov_ex2_alleles	individual_id	Amb_ex2_001	Amb_ex2_002	Amb_ex2_003	Amb_ex2_004	Amb_ex2_005	Amb_ex2_006	Amb_ex2_007	Amb_ex2_008
#3425	2429	14829	0	0	0	0	0	0	0	0
#4388	3852	14830	0	407	0	0	0	0	0	351
#6938	5059	14831	0	0	0	0	0	0	0	540
#returns a huge table with all alleles and all individuals, NAs coverted to 0 
get_gen <- function(genera, exon){
  a <- NULL
  for(g in genera){
    inp <- read.table(paste0(exon, "_genotypes_", g, ".txt"), header = T, sep= "\t")
    if(is.null(a)){
      a <- inp
    } else {
      a <- full_join(a, inp)
    }
  }
  a[is.na(a)] <- 0
  return(a)
}


#takes object produced by get_binary_annotated, genus, species,
#list of alleles of interest (for example all functional)
# and returns dataframe containing annotation columns (ctr - make sure that you control for exon)
#as well as genotypes for all alleles that are present in the taxon of interest
get_tax_genot <- function(d, gen, sp, alleles_to_include, 
                          ctr = c("individual_id", "family", "genus", "species", "Locality", 
                                  paste0("cov_", ex), paste0("cov_", ex, "_alleles"))){
  dtax <- d %>% filter(genus == gen & species == sp) %>% select_if(~ !is.numeric(.) || sum(., na.rm = T ) > 0)
  dtax_alle_set <- dtax %>% select(all_of(ctr), any_of(alleles_to_include))
  return(dtax_alle_set)
}


#calculate number of alleles or private alleles in each group (column-defined)
N_all_tot_priv <- function(inp, grp){
  g <- enquo(grp)
  #for each allele presence/absence (T/F) in each species
  grp_all <- inp %>% group_by(!!g) %>%
    summarize_at(vars(matches(alleles$allele)), list(~ as.logical(sum(.>0, na.rm=TRUE))))
  N_priv_all <- grp_all %>%
    select_if(~ !is.logical(.) || sum(.) == 1) %>%
    transmute(grp = !!g,
              Npriv_I_HEX = pmap_int(select(., one_of(IHEX)), sum),
              Npriv_I_LEX = pmap_int(select(., one_of(ILEX)), sum),
              Npriv_II_HEX = pmap_int(select(., one_of(IIHEX)), sum),
              Npriv_II_LEX = pmap_int(select(., one_of(IILEX)), sum))
  N_all_all <- grp_all %>%
    transmute(grp = !!g,
              Nall_I_HEX = pmap_int(select(., one_of(IHEX)), sum),
              Nall_I_LEX = pmap_int(select(., one_of(ILEX)), sum),
              Nall_II_HEX = pmap_int(select(., one_of(IIHEX)), sum),
              Nall_II_LEX = pmap_int(select(., one_of(IILEX)), sum))
  Nall_priv <- left_join(N_all_all, N_priv_all, by = "grp")
  return(Nall_priv)
}


read_data <- function(inp){
  d <- read.table(paste0("data/", inp), header = T, sep = "\t", encoding="UTF-8")
  return(d)
}

#generate vector of all allele names from a given category
sel_all_cat <- function(cl, cat){
  cl_cat <- alleles %>% filter(class == cl, category == cat) %>% select(allele) %>% as_vector()
  return(cl_cat)
}

write_lnx_head <- function(inp, outp){
  of <- file(outp, "wb")
  write.table(inp, file = of, quote = F, row.names = F, col.names= T, eol="\n", sep = "\t")
  close(of)
}

write_struct <- function(inp, outp){
  of <- file(outp, "wb")
  t <- inp %>% slice(rep(1:n(), each=2))
  r <- c(rep("", times = 3), rep(0, times = (length(inp)-3)))
  w <- rbind(r, t)
  write.table(w, file = of, quote = F, row.names = F, col.names= F, eol="\n", sep = "\t")
  close(of)
}
