# Setup

if (FALSE) {
  # Clean session
  rm(list = ls())
  
  # Install
  setRepositories() # choose CRAN and BioC
  install.packages(pkg = c("gRain", "devtools"), dep = TRUE)
  # PR that fixed the repeatPattern() data issue
  # https://github.com/hojsgaard/gRain/pull/9
  devtools::install_github(repo = "https://github.com/gregorgorjanc/gRain")
  
  # Load package
  library(package = "gRain")
}

# User defined functions

#' @rdname CheckPedigree
#' @title Check pedigree.
#'
#' @description Provides three basic checks for pedigree:
#'   any duplicated individuals,any parents that are not listed as
#'   individuals, and any parents that are mixed (fathers appear as
#'   mothers, and vice versa).
#'
#' @param x data.frame, with columns IID (individual ID), FID (father ID), and
#'   MID (mother ID) - missing parents are encoded as \code{NA}
#' @param duplicated logical, check for duplicated individuals.
#' @param listedParents logical, check if parents are listed as individuals.
#' @param mixedParents logical, check if parents are mixed (fathers appear as
#'   mothers or vice versa).
#'
#' @return Nothing, the function throws an error when pedigree has issues.
#'
#' @examples
#' ped <- data.frame(IID = c(1, 2, 3, 4, 5),
#'                   FID = c(0, 0, 1, 1, 3),
#'                   MID = c(0, 0, 2, 2, 4))
#' CheckPedigree(x= ped, duplicated = TRUE, listedParents = TRUE, mixedParents = TRUE)

CheckPedigree <- function(x,
                          duplicated = TRUE,
                          listedParents = TRUE,
                          mixedParents = TRUE) {
  # Check the data frame has columns IID, FID, and MID
  test <- !c("IID", "FID", "MID") %in% names(x)
  if (any(test)) {
    stop("The pedigreeData should have columns IID, FID, and MID!")
  }
  # CheckPedigree
  if (duplicated) {
    test <- anyDuplicated(x[[1]])
    if (test > 0) {
      stop("Duplicated records (IID) in pedigree not allowed!")
    }
  }
  if (listedParents) {
    test <- !is.na(x[[2]]) & !x[[2]] %in% x[[1]]
    if (any(test)) {
      stop("All fathers must be listed as individuals with their own row in pedigree!")
    }
    test <- !is.na(x[[3]]) & !x[[3]] %in% x[[1]]
    if (any(test)) {
      stop("All mothers must be listed as individuals with their own row in pedigree!")
    }
  }
  if (mixedParents) {
    test <- !is.na(x[[3]]) & x[[3]] %in% x[[2]]
    if (any(test)) {
      stop("Mothers must not appear as fathers in pedigree!")
    }
    test <- !is.na(x[[2]]) & x[[2]] %in% x[[3]]
    if (any(test)) {
      stop("Fathers must not appear as mothers in pedigree!")
    }
  }
}

#' @rdname PedigreeGraph
#' @title Convert pedigree data.frame to a DAG.
#'
#' @description Pedigree has a natural representation as a Directed Acyclic
#'   Graph (DGA) where parent nodes are connected to progeny nodes.
#'
#' @param x data.frame, with columns IID (individual ID), FID (father ID), and
#'   MID (mother ID).
#' @param unknown NA/numeric/character, symbol(s) that denote unknown parent(s).
#'
#' @details It is assumed that IID, FID, and MID columns are characters and this
#'   is enforced! Will throw an error if any loops in pedigree (individuals being
#'   parents of themselves).
#'
#' @return igraph graph.
#'
#' @examples
#' ped <- data.frame(IID = c(1, 2, 3, 4, 5),
#'                   FID = c(0, 0, 1, 1, 3),
#'                   MID = c(0, 0, 2, 2, 4))
#' ret <- PedigreeGraph(x = ped, unknown = 0)
#' plot(ret)

PedigreeGraph <- function(x, unknown = NA) {
  x[] <- lapply(X = x, FUN = as.character)
  # x data.frame(IID, FID, MID)
  # Father --> Individual edges c(A, C)
  Sel <- !x[[2]] %in% unknown
  FI <- cbind(x[[2]][Sel], x[[1]][Sel])
  # Mother --> Individual edges c(B, C)
  Sel <- !x[[3]] %in% unknown
  MI <- cbind(x[[3]][Sel], x[[1]][Sel])
  # Graph
  igraph::graph_from_edgelist(el = rbind(FI, MI),
                              directed = TRUE)
}

#' @rdname SortPedigree
#' @title Sort pedigree such that parents precede progeny.
#'
#' @description Many pedigree operations require that parents precede progeny.
#'   We can achieve this with topological sort of the pedigree by recognizing
#'   that pedigree can be represented as a DAG.
#'
#' @param x data.frame, with columns IID (individual ID), FID (father ID), and
#'   MID (mother ID).
#' @param unknown NA/numeric/character, symbol(s) that denote unknown parent(s)
#' @param mode character, should the pedigree be sorted from founders to
#'   non-founders ("out") or the other way around ("in").
#'
#' @return sorted \code{x}
#'
#' @examples
#' ped <- data.frame(IID = c(1, 2, 3, 4, 5),
#'                   FID = c(0, 0, 1, 1, 3),
#'                   MID = c(0, 0, 2, 2, 4))
#' ped <- ped[c(5, 1, 3, 4, 2), ]
#' ped
#' SortPedigree(x = ped, unknown = 0)
#' ped[] <- lapply(ped, FUN = as.character)
#' SortPedigree(x = ped, unknown = 0)

SortPedigree <- function(x, unknown = NA, mode = "out") {
  Order <- PedigreeGraph(x = x, unknown = unknown) |>
    igraph::topo_sort(mode = mode) |>
    names() |>
    match(table = x[[1]])
  x[Order,]
}

#' @rdname RecodePedigree
#' @title Recode pedigree ID to 1:n.
#'
#' @description Many pedigree operations require integer ID in 1:n.
#'
#' @param x data.frame, with columns IID (individual ID), FID (father ID), and
#'   MID (mother ID).
#'
#' @return data.frame with three columns.
#'
#' @examples
#' ped <- data.frame(IID = c("A", "B", "C", "D", "E"),
#'                   FID = c( NA,  NA, "A", "A", "C"),
#'                   MID = c( NA,  NA, "B", "B", "D"))
#' RecodePedigree(ped)
#' ped <- ped[c(5, 1, 3, 4, 2), ]
#' ped
#' RecodePedigree(ped)
#' ped <- SortPedigree(x = ped)
#' RecodePedigree(ped)

RecodePedigree <- function(x) {
  data.frame(
    IID = 1:nrow(x),
    FID = match(
      x = x[[2]],
      table = x[[1]],
      nomatch = NA
    ),
    MID = match(
      x = x[[3]],
      table = x[[1]],
      nomatch = NA
    )
  )
}

#' @rdname AlleleFetch
#' @title Genotype and phenotype probability calculations on a pedigree.
#'
#' @description Genotype and phenotype probability calculations on a pedigree
#'   with a considerable flexibility in terms of number of loci
#'   and gene action involved in the trait(s) of interest.
#'
#' @param pedigreeData data.frame, with columns IID (individual ID), FID (father
#'   ID, NA if unknown), MID (mother ID, NA if unknown), and YOB (year of birth,
#'   NA if unknown).
#' @param genotypeState character, list of allowed/expected genotypes.
#' @param genotypeData data.frame, with columns IID (individual ID; all
#'   individuals in this data.frame must appear in pedigree too) and Genotype
#'   (character).
#' @param phenotypeState character, list of allowed/expected phenotypes.
#' @param phenotypeData data.frame, with columns IID (individual ID; all
#'   individuals in this data.frame must appear in pedigree too) and Phenotype
#'   (character).
#' @param genoPT numeric, vector of genotype probabilities (marginal) in founders.
#' @param mendelianCPT numeric, vector of Mendelian transmission probabilities
#'   (conditional) between two parents and their progeny.
#' @param mendelianCPTOnlyFather numeric, vector of Mendelian transmission
#'   probabilities (conditional) between one parent (father) and its progeny.
#' @param mendelianCPTOnlyMother numeric, vector of Mendelian transmission
#'   probabilities (conditional) between one parent (mother) and its progeny.
#' @param genoPenetranceCPT numeric, vector of genotype penetrance
#'   probabilities (conditional) between true and observed genotype.
#' @param phenoPenetranceCPT numeric, vector of phenotype penetrance
#'   probabilities (conditional) between true genotype and phenotype.
#'
#' @return grain object (compiled and with propagated evidence).
#'
#' @examples
#' testPedigreeData <- structure(list(IID = c(5788956L, 5166631L, 9880745L, 6942882L, 5538891L, 9073689L, 6632232L, 5577159L, 5113284L, 6804872L),
#'                                    FID = c(NA, NA, NA, 5166631L, 5166631L, 6942882L, 5166631L, 5166631L, 6942882L, 6942882L),
#'                                    MID = c(NA, NA, 5788956L, NA, 9880745L, 9880745L, 5788956L, 9073689L, 5577159L, 9073689L)),
#'                               class = "data.frame", row.names = c(NA, -10L))
#'
#' testGenotypeData <- structure(list(IID = c(5538891L, 6632232L),
#'                                    Genotype = c("AA", "AB")),
#'                               class = "data.frame", row.names = c(NA, -2L))
#'
#' testPhenotypeData <- structure(list(IID = c(9880745L, 6804872L),
#'                                     Phenotype = c("OK", "NOK")),
#'                                class = "data.frame", row.names = c(NA, -2L))
#'
#' # Set the genotype values
#' testGenotypeState <- c("AA", "AB", "BB")
#'
#' # Prior allele frequency in (founder) population
#' PrA <- 2/3
#' PrB <- 1 - PrA
#'
#' # Prior genotype frequencies in (founder) population
#' # according to the Hardy-Weinberg's expectation
#' testGenoPT <- c(PrA*PrA, 2*PrA*PrB, PrB*PrB)
#'
#' # Genotype frequencies in individual given the genotype of parents according to Mendel's law
#' # and accounting for mutation (via e)
#' e <- 0.0001
#' gM_AA_AA <- c(1-2*e+e^2, 2*e-2*e^2, e^2)
#' gM_AA_AB <- c(1/2-1/2*e, 1/2, 1/2*e)
#' gM_AA_BB <- c(e-e^2, 1-2*e+2*e^2, e-e^2)
#' gM_AB_AA <- c(1/2-1/2*e, 1/2, 1/2*e)
#' gM_AB_AB <- c(1/4, 1/2, 1/4)
#' gM_AB_BB <- c(1/2*e, 1/2, 1/2-1/2*e)
#' gM_BB_AA <- c(e-e^2, 1-2*e+2*e^2, e-e^2)
#' gM_BB_AB <- c(1/2*e, 1/2, 1/2-1/2*e)
#' gM_BB_BB <- c(e^2, 2*e-2*e^2, 1-2*e+e^2)
#'
#' testMendelianCPT <- matrix(nrow=3, ncol=9)
#' testMendelianCPT[, 1] <- gM_AA_AA
#' testMendelianCPT[, 2] <- gM_AA_AB
#' testMendelianCPT[, 3] <- gM_AA_BB
#' testMendelianCPT[, 4] <- gM_AB_AA
#' testMendelianCPT[, 5] <- gM_AB_AB
#' testMendelianCPT[, 6] <- gM_AB_BB
#' testMendelianCPT[, 7] <- gM_BB_AA
#' testMendelianCPT[, 8] <- gM_BB_AB
#' testMendelianCPT[, 9] <- gM_BB_BB
#'
#' # Genotype frequencies in individual given the genotype of one parent (Mendel's law with mutation)
#' # and base population (Hardy-Weinberg's law)
#' gM_AA_unk <- c(PrA-PrA*e, PrB-PrB*e+PrA*e, PrB*e)
#' gM_AB_unk <- c(1/2*PrA, 1/2*(PrA+PrB), 1/2*PrB)
#' gM_BB_unk <- c(PrA*e, PrA-PrA*e+PrB*e, PrB-PrB*e)
#'
#' testMendelianCPTOnlyFather <- matrix(nrow=3, ncol= 3)
#' testMendelianCPTOnlyFather[,1] <- gM_AA_unk
#' testMendelianCPTOnlyFather[,2] <- gM_AB_unk
#' testMendelianCPTOnlyFather[,3] <- gM_BB_unk
#'
#' testMendelianCPTOnlyMother <- testMendelianCPTOnlyFather
#'
#' # Set of phenotype values
#' testPhenotypeState <- c("OK", "NOK")
#'
#' # Phenotype frequences in individual given the genotype of the individual - phenotype penetrance
#' e <- 0.01
#' testPhenoPenetranceCPT <- matrix(nrow=2, ncol=3)
#' #                               OK   NOK Pheno
#' testPhenoPenetranceCPT[,1] <- c(1-e, 0+e) # AA True geno
#' testPhenoPenetranceCPT[,2] <- c(1-e, 0+e) # AB
#' testPhenoPenetranceCPT[,3] <- c(0+e, 1-e) # BB
#'
#' # Genotype frequences in individual given the genotype of the individual - genotype penetrance
#' e <- 0.001
#' testGenoPenetranceCPT <- matrix(nrow=3, ncol=3)
#' #                             AA       AB     BB Observed Geno
#' testGenoPenetranceCPT[,1] <-c(1-e-e^2, e,     e^2)     # AA True Geno
#' testGenoPenetranceCPT[,2] <-c(e,       1-2*e, e)       # AB
#' testGenoPenetranceCPT[,3] <-c(e^2,     e,     1-e-e^2) # BB
#'
#' out <- AlleleFetch(pedigreeData = testPedigreeData,
#'                    genotypeState = testGenotypeState,
#'                    genotypeData = testGenotypeData,
#'                    phenotypeState = testPhenotypeState,
#'                    phenotypeData = testPhenotypeData,
#'                    genoPT = testGenoPT,
#'                    mendelianCPT = testMendelianCPT,
#'                    mendelianCPTOnlyFather = testMendelianCPTOnlyFather,
#'                    mendelianCPTOnlyMother = testMendelianCPTOnlyMother,
#'                    genoPenetranceCPT = testGenoPenetranceCPT,
#'                    phenoPenetranceCPT = testPhenoPenetranceCPT)
#' out$result
#' plot(out$model)

if (FALSE) {
  pedigreeData = testPedigreeData
  genotypeState = testGenotypeState
  genotypeData = testGenotypeData
  phenotypeState = testPhenotypeState
  phenotypeData = testPhenotypeData
  genoPT = testGenoPT
  mendelianCPT = testMendelianCPT
  mendelianCPTOnlyFather = testMendelianCPTOnlyFather
  mendelianCPTOnlyMother = testMendelianCPTOnlyMother
  genoPenetranceCPT = testGenoPenetranceCPT
  phenoPenetranceCPT = testPhenoPenetranceCPT
}
AlleleFetch <- function(pedigreeData = NULL,
                        genotypeState = NULL,
                        genotypeData = NULL,
                        phenotypeState = NULL,
                        phenotypeData = NULL,
                        genoPT,
                        mendelianCPT,
                        mendelianCPTOnlyFather,
                        mendelianCPTOnlyMother,
                        genoPenetranceCPT = NULL,
                        phenoPenetranceCPT = NULL,
                        agePhenotype = NULL,
                        phenoPenetranceCPT_Before=NULL,
                        ...) {
  # ---- Data quality checks ----
  
  # Pedigree data quality check
  if (is.null(pedigreeData)) {
    stop("pedigreeData must be provided!")
  }
  CheckPedigree(x = pedigreeData)
  pedigreeData <- SortPedigree(x = pedigreeData)
  pedigreeDataInt <- RecodePedigree(x = pedigreeData)
  colnames(pedigreeDataInt) <-
    paste0(colnames(pedigreeDataInt), "Int")
  
  # Genotype data quality check
  if (is.null(genotypeState)) {
    stop("genotypeState must be provided!")
  }
  if (!is.null(genotypeData)) {
    # Data frame check in GenotypeData
    test <- !"IID" %in% names(genotypeData)
    if (test) {
      stop("There is no IID column in the genotypeData!")
    }
    test <- !"Genotype" %in% names(genotypeData)
    if (test) {
      stop("There is no Genotype column in the genotypeData!")
    }
    test <- !genotypeData$IID %in% pedigreeData$IID
    if (any(test)) {
      stop("All individuals in genotypeData must appear in pedigreeData!")
    }
    test <- !genotypeData$Genotype %in% genotypeState
    if (any(test)) {
      stop("Invalid genotypes in genotypeData!")
    }
    test <- anyDuplicated(genotypeData$IID)
    if (test > 0) {
      stop("Duplicated records (IID) in genotypeData not allowed!")
    }
    genotypeData$IIDInt <- match(x = genotypeData$IID,
                                 table = pedigreeData$IID)
  }
  
  # Phenotype data quality check
  if (is.null(phenotypeState)) {
    stop("phenotypeState must be provided!")
  }
  if (!is.null(phenotypeData)) {
    # Data frame check in phenotypeData
    test <- !"IID" %in% names(phenotypeData)
    if (test) {
      stop("There is no IID column in the phenotypeData!")
    }
    test <- !"Phenotype" %in% names(phenotypeData)
    if (test) {
      stop("There is no Phenotype column in the phenotypeData!")
    }
    test <- !phenotypeData$IID %in% pedigreeData$IID
    if (any(test)) {
      stop("All individuals in phenotypeData must appear in pedigreeData!")
    }
    test <- !phenotypeData$Phenotype %in% phenotypeState
    if (any(test)) {
      stop("Invalid phenotypes in the phenotypeData!")
    }
    test <- anyDuplicated(phenotypeData$IID)
    if (test > 0) {
      stop("Duplicated records (IID) in phenotypeData not allowed!")
    }
    phenotypeData$IIDInt <- match(x = phenotypeData$IID,
                                  table = pedigreeData$IID)
  }
  # agePhenotype sorting
  if(!is.null(agePhenotype)){
    # Match the reordered pedigree with the age input
    agePhenotype$IIDInt <- match(x = agePhenotype$IID,
                                 table = pedigreeData$IID)
    # Merge age with phenotype data
    phenoAge <- merge(phenotypeData, agePhenotype, by="IID")
    
    # Separate into young and old at phenotype test
    # Younger than 4, older than 3 (so boundary is 4 years age)
    agePhenotypeYg <- phenoAge[phenoAge$Age < 4,]
    agePhenotypeOld <- phenoAge[phenoAge$Age > 3,]
    
    # Separate into the phenotypes, grouped by age.
    # Observed phenotypes aged below 4 years old
    ygPheno <- agePhenotypeYg[,c(1,2),]
    ygPheno$IIDInt <- match(x = ygPheno$IID,
                            table = pedigreeData$IID)
    # Observed phenotypes aged above 3 years old
    oldPheno <-agePhenotypeOld[,c(1,2),]
    oldPheno$IIDInt <- match(x = oldPheno$IID,
                             table = pedigreeData$IID)
    # When the observed phenotype was taken at an unknown age
    test <- rbind(ygPheno, oldPheno)
    NoAgePheno <- phenotypeData[!phenotypeData$IIDInt %in% test$IIDInt,]
  }
  # if no age was insertted into the function
  if(is.null(agePhenotype)){
    NoAgePheno <- phenotypeData
  }
  
  
  # ---- Model definition ----
  
  naFID <- is.na(pedigreeDataInt$FIDInt)
  naMID <- is.na(pedigreeDataInt$MIDInt)
  
  # True/inferred genotypes in founders (both parents unknown)
  genoFoundersCPT <- cptable( ~ geno[i],
                              values = genoPT,
                              levels = genotypeState)
  genoFoundersCPT <- repeatPattern(list(genoFoundersCPT),
                                   instances = pedigreeDataInt[naFID &
                                                                 naMID, "IIDInt"])
  allCPT <- genoFoundersCPT
  
  # True/inferred genotypes in founders (only father known)
  test <- !naFID & naMID
  if (any(test)) {
    genoOnlyFatherKnownCPT <-
      cptable( ~ geno[i] | geno[data[i, "FIDInt"]],
               values = mendelianCPTOnlyFather,
               levels = genotypeState)
    # browser()
    genoOnlyFatherKnownCPT <-
      repeatPattern(list(genoOnlyFatherKnownCPT),
                    instances = pedigreeDataInt[test, "IIDInt"],
                    data = pedigreeDataInt)
    allCPT <- c(allCPT, genoOnlyFatherKnownCPT)
  }
  
  # True/inferred genotypes in founders (only mother known)
  test <- naFID & !naMID
  if (any(test)) {
    genoOnlyMotherKnownCPT <-
      cptable( ~ geno[i] | geno[data[i, "MIDInt"]],
               values = mendelianCPTOnlyMother,
               levels = genotypeState)
    genoOnlyMotherKnownCPT <-
      repeatPattern(list(genoOnlyMotherKnownCPT),
                    instances = pedigreeDataInt[test, "IIDInt"],
                    data = pedigreeDataInt)
    allCPT <- c(allCPT, genoOnlyMotherKnownCPT)
  }
  
  # True/inferred genotypes in non-founders
  test <- !naFID & !naMID
  if (any(test)) {
    genoNonFoundersCPT <- cptable( ~ geno[i] | geno[data[i, "FIDInt"]] +
                                     geno[data[i, "MIDInt"]],
                                   values = mendelianCPT,
                                   levels = genotypeState)
    genoNonFoundersCPT <- repeatPattern(list(genoNonFoundersCPT),
                                        instances = pedigreeDataInt[test, "IIDInt"],
                                        data = pedigreeDataInt)
    allCPT <- c(allCPT, genoNonFoundersCPT)
  }
  
  # True/inferred phenotypes
  if (!is.null(phenotypeState)) {
    phenoCPT <- cptable( ~ pheno[i] | geno[i],
                         values = phenoPenetranceCPT,     # Using the standard phenotype penetrance only
                         levels = phenotypeState)
    phenoCPT <- repeatPattern(list(phenoCPT),
                              instances = pedigreeDataInt$IIDInt) # for all pedigreeData individuals
    allCPT <- c(allCPT, phenoCPT)
  }
  
  # Observed genotypes (=data)
  if (!is.null(genotypeData)) {
    obsGenoCPT <- cptable( ~ obsGeno[i] | geno[i],
                           values = genoPenetranceCPT,
                           levels = genotypeState)
    obsGenoCPT <- repeatPattern(list(obsGenoCPT),
                                instances = genotypeData$IIDInt)
    allCPT <- c(allCPT, obsGenoCPT)
  }
  
  # Observed phenotypes (=data)
  if (!is.null(phenotypeData)) {
    # Age-adjusted penetrance
    if (!is.null(agePhenotype)){
      # Test if any observations in OldPheno
      test <- as.numeric(nrow(oldPheno))
      # If oldPheno has 1 or more observations
      if(!test == 0){
        # Use the standard penetrance to integrate the observed phenotypes
        obsPhenoCPT <- cptable( ~ obsPheno[i] | geno[i],
                                values = phenoPenetranceCPT,
                                levels = phenotypeState)
        obsPhenoCPT <- repeatPattern(list(obsPhenoCPT),
                                     instances = oldPheno$IIDInt)
        allCPT <- c(allCPT, obsPhenoCPT)
      }
      test <- as.numeric(nrow(ygPheno))
      # If ygPheno has 1 or more observations
      if(!test == 0){
        # Use the age-adjusted penetrance to integrate the observed phenotypes
        obsPhenoCPT <- cptable( ~ obsPheno[i] | geno[i],
                                values = phenoPenetranceCPT_Before,
                                levels = phenotypeState)
        obsPhenoCPT <- repeatPattern(list(obsPhenoCPT),
                                     instances = ygPheno$IIDInt)
        allCPT <- c(allCPT, obsPhenoCPT)
      }
      test <- as.numeric(nrow(NoAgePheno))
      # If the age is unknown on observed phenotype
      if(!test == 0){
        obsPhenoCPT <- cptable( ~ obsPheno[i] | geno[i],
                                values = phenoPenetranceCPT,
                                levels = phenotypeState)
        obsPhenoCPT <- repeatPattern(list(obsPhenoCPT),
                                     instances = NoAgePheno$IIDInt)
        allCPT <- c(allCPT, obsPhenoCPT)
      }
    }
    # If the age at observed phenotype is not included in the function
    if(is.null(agePhenotype)){
      obsPhenoCPT <- cptable( ~ obsPheno[i] | geno[i],
                              values = phenoPenetranceCPT,
                              levels = phenotypeState)
      obsPhenoCPT <- repeatPattern(list(obsPhenoCPT),
                                   instances = phenotypeData$IIDInt)
      allCPT <- c(allCPT, obsPhenoCPT)
    }
  }

  
  allCPT <- compileCPT(x = allCPT)
  model <- grain(x = allCPT)
  
  # ---- Set data ----
  
  # Observed genotypes
  if (!is.null(genotypeData)) {
    model <- setEvidence(object = model,
                         nodes = paste("obsGeno", genotypeData$IIDInt, sep = ""),
                         states = genotypeData$Genotype)
  }
  
  # Observed phenotypes
  if (!is.null(phenotypeData)) {
    model <- setEvidence(object = model,
                         nodes = paste("obsPheno", phenotypeData$IIDInt, sep = ""),
                         states = phenotypeData$Phenotype)
  }
  
  # ---- Propagate data through the network and save results ----
  
  posterior <- qgrain(object = model, result = "data.frame")
  
  posteriorNames <- names(posterior)
  selGeno <- grepl(pattern = "geno", x = posteriorNames)
  posteriorGeno <- posterior[selGeno]
  posteriorGenoNames <- names(posteriorGeno)
  tmp <- names(posteriorGeno[[1]])
  tmp[1] <- "State"
  posteriorGeno <- lapply(X = posteriorGeno, FUN = function(x) setNames(x, tmp))
  posteriorGeno <- do.call(what = rbind, args = posteriorGeno)
  posteriorGeno$IIDInt <- rep(x = gsub(pattern = "geno", replacement = "",
                                       x = posteriorGenoNames),
                              each = length(genotypeState)) |> as.numeric() 
  posteriorGeno$Var <- "geno"
  
  selPheno <- grepl(pattern = "pheno", x = posteriorNames)
  posteriorPheno <- posterior[selPheno]
  posteriorPhenoNames <- names(posteriorPheno)
  tmp <- names(posteriorPheno[[1]])
  tmp[1] <- "State"
  posteriorPheno <- lapply(X = posteriorPheno, FUN = function(x) setNames(x, tmp))
  posteriorPheno <- do.call(what = rbind, args = posteriorPheno)
  posteriorPheno$IIDInt <- rep(x = gsub(pattern = "pheno", replacement = "",
                                        x = posteriorPhenoNames),
                               each = length(phenotypeState)) |> as.numeric()
  posteriorPheno$Var <- "pheno"
  
  ret <- rbind(posteriorGeno, posteriorPheno)
  ret <- merge(x = cbind(pedigreeData, pedigreeDataInt), y = ret)
  
  # Reorder rows based on IIDInt and Var
  ret <- ret[order(ret$IIDInt, ret$Var), ]
  
  # ---- Return ----
  return(list(model = model,
              result = ret))
}


#' @rdname DisplayResults
#' @title Converting the Allelefetch()$Result into a wide data frame.
#'
#' @description The conversion of the result of AlleleFetch from a tall to a
#'    wide data frame with the ability to
#'    potentially include the observed inputs.
#'
#' @param x data.frame, the full output from the AlleleFetch function.
#' @param withObserv logical, asks whether the observed values are to be
#' included in the results.
#'
#' @details Requires AlleleFetch()$Result, and the Tidyverse Package.
#'
#' @return Wide dataframe of AlleleFetch()$Results.
#'
#' @examples
#'  testPedigreeData <- structure(list(IID = c(5788956L, 5166631L, 9880745L, 6942882L, 5538891L, 9073689L, 6632232L, 5577159L, 5113284L, 6804872L),
#'                                    FID = c(NA, NA, NA, 5166631L, 5166631L, 6942882L, 5166631L, 5166631L, 6942882L, 6942882L),
#'                                    MID = c(NA, NA, 5788956L, NA, 9880745L, 9880745L, 5788956L, 9073689L, 5577159L, 9073689L)),
#'                               class = "data.frame", row.names = c(NA, -10L))
#'
#' testGenotypeData <- structure(list(IID = c(5538891L, 6632232L),
#'                                    Genotype = c("AA", "AB")),
#'                               class = "data.frame", row.names = c(NA, -2L))
#'
#' testPhenotypeData <- structure(list(IID = c(9880745L, 6804872L),
#'                                     Phenotype = c("OK", "NOK")),
#'                                class = "data.frame", row.names = c(NA, -2L))
#'
#' testGenotypeState <- c("AA", "AB", "BB")
#' testPhenotypeState <- c("OK", "NOK")
#'
#' # Prior allele frequency in (founder) population
#' PrA <- 2/3
#' PrB <- 1 - PrA
#'
#' # Prior genotype frequencies in (founder) population
#' # according to the Hardy-Weinberg's expectation
#' testGenoPT <- c(PrA*PrA, 2*PrA*PrB, PrB*PrB)
#'
#' # Genotype frequencies in individual given the genotype of parents according to Mendel's law
#' # and accounting for mutation (via e)
#' e <- 0.0001
#' gM_AA_AA <- c(1-2*e+e^2, 2*e-2*e^2, e^2)
#' gM_AA_AB <- c(1/2-1/2*e, 1/2, 1/2*e)
#' gM_AA_BB <- c(e-e^2, 1-2*e+2*e^2, e-e^2)
#' gM_AB_AA <- c(1/2-1/2*e, 1/2, 1/2*e)
#' gM_AB_AB <- c(1/4, 1/2, 1/4)
#' gM_AB_BB <- c(1/2*e, 1/2, 1/2-1/2*e)
#' gM_BB_AA <- c(e-e^2, 1-2*e+2*e^2, e-e^2)
#' gM_BB_AB <- c(1/2*e, 1/2, 1/2-1/2*e)
#' gM_BB_BB <- c(e^2, 2*e-2*e^2, 1-2*e+e^2)
#'
#' testMendelianCPT <- matrix(nrow=3, ncol=9)
#' testMendelianCPT[, 1] <- gM_AA_AA
#' testMendelianCPT[, 2] <- gM_AA_AB
#' testMendelianCPT[, 3] <- gM_AA_BB
#' testMendelianCPT[, 4] <- gM_AB_AA
#' testMendelianCPT[, 5] <- gM_AB_AB
#' testMendelianCPT[, 6] <- gM_AB_BB
#' testMendelianCPT[, 7] <- gM_BB_AA
#' testMendelianCPT[, 8] <- gM_BB_AB
#' testMendelianCPT[, 9] <- gM_BB_BB
#'
#' # Genotype frequencies in individual given the genotype of one parent (Mendel's law with mutation)
#' # and base population (Hardy-Weinberg's law)
#' gM_AA_unk <- c(PrA-PrA*e, PrB-PrB*e+PrA*e, PrB*e)
#' gM_AB_unk <- c(1/2*PrA, 1/2*(PrA+PrB), 1/2*PrB)
#' gM_BB_unk <- c(PrA*e, PrA-PrA*e+PrB*e, PrB-PrB*e)
#'
#' testMendelianCPTOnlyFather <- matrix(nrow=3, ncol= 3)
#' testMendelianCPTOnlyFather[,1] <- gM_AA_unk
#' testMendelianCPTOnlyFather[,2] <- gM_AB_unk
#' testMendelianCPTOnlyFather[,3] <- gM_BB_unk
#'
#' testMendelianCPTOnlyMother <- testMendelianCPTOnlyFather
#'
#' # Set of phenotype values
#' testPhenotypeState <- c("OK", "NOK")
#'
#' # Phenotype frequences in individual given the genotype of the individual - phenotype penetrance
#' e <- 0.01
#' testPhenoPenetranceCPT <- matrix(nrow=2, ncol=3)
#' #                               OK   NOK Pheno
#' testPhenoPenetranceCPT[,1] <- c(1-e, 0+e) # AA True geno
#' testPhenoPenetranceCPT[,2] <- c(1-e, 0+e) # AB
#' testPhenoPenetranceCPT[,3] <- c(0+e, 1-e) # BB
#'
#' # Genotype frequences in individual given the genotype of the individual - genotype penetrance
#' e <- 0.001
#' testGenoPenetranceCPT <- matrix(nrow=3, ncol=3)
#' #                             AA       AB     BB Observed Geno
#' testGenoPenetranceCPT[,1] <-c(1-e-e^2, e,     e^2)     # AA True Geno
#' testGenoPenetranceCPT[,2] <-c(e,       1-2*e, e)       # AB
#' testGenoPenetranceCPT[,3] <-c(e^2,     e,     1-e-e^2) # BB
#'
#' out <- AlleleFetch(pedigreeData = testPedigreeData,
#'                    genotypeState = testGenotypeState,
#'                    genotypeData = testGenotypeData,
#'                    phenotypeState = testPhenotypeState,
#'                    phenotypeData = testPhenotypeData,
#'                    genoPT = testGenoPT,
#'                    mendelianCPT = testMendelianCPT,
#'                    mendelianCPTOnlyFather = testMendelianCPTOnlyFather,
#'                    mendelianCPTOnlyMother = testMendelianCPTOnlyMother,
#'                    genoPenetranceCPT = testGenoPenetranceCPT,
#'                    phenoPenetranceCPT = testPhenoPenetranceCPT)
#' DisplayResults(x = out, withObserv = TRUE)

library(tidyverse)

DisplayResults <- function(x, withObserv = TRUE) {
  # Tall to Wide dataframe (using tidyverse)
  data_wide <- x$result %>%
    pivot_wider(names_from = c(State, Var),
                values_from = Freq)
  #For inclusion of observed genotypes and phenotypes.
  if (withObserv == TRUE) {
    # Create data frame of observed genotypes and phenotypes.
    test <- as.data.frame(x$model$evidence)
    # Data frame for AlleleFetch results (all individuals).
    tmp <- x$result
    # Separate the genotype states.
    tmpGeno <- tmp[tmp$Var == "geno",]
    tmpGeno <- as.vector(unique(tmpGeno$State))
    # Separate phenotype states.
    tmpPheno <- tmp[tmp$Var == "pheno",]
    tmpPheno <- as.vector(unique(tmpPheno$State))
    # Split the hard state by genotype states.
    test2 <- test %>%
      filter(hard.state %in% tmpPheno)
    test2 <- test2 %>% rename(ObsPheno = "hard.state")
    test2$IIDInt <- substr(test2$nodes, 9, nchar(test2$nodes))
    # Split the hard state by phenotype states.
    test3 <- test %>%
      filter(hard.state %in% tmpGeno)
    test3 <- test3 %>% rename(ObsGeno = "hard.state")
    test3$IIDInt <- substr(test3$nodes, 8, nchar(test3$nodes))
    # Combine the observed genotypes and phenotypes into one data frame.
    data <- merge(test2, test3, by = "IIDInt", all = TRUE)
    # Select only IIDInt, ObsPheno, and ObsGeno columns.
    data <- data %>%
      select(IIDInt, ObsPheno, ObsGeno)
    # Merge with the wide data table and return results.
    result_obs <- merge(data_wide, data, by = "IIDInt", all = TRUE)
    return(result_obs)
  } else {
    return(data_wide)
  }
  