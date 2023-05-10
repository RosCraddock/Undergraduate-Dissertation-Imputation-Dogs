# Undergraduate-Dissertation-Imputation-Dogs
 The full code used for the dissertation in fulfilment of an undergraduate degree in Agricultural Science (Animal Science) at the University of Edinburgh

Title: Evaluating the Potential of Three Models for Probabilistic Pedigree-based Imputation in the Management of Monogenic Recessive Diseases in KC-Registered Dogs


List & Description of Files:

PBGVPedigreeCleaning: This outlines the method for cleaning and preparing the data (see Section 3.X.X in the dissertation), manually completed in R then later checked using a Fortran programme named ...

FullPedigreeAnalysis: This outlines the code used for the analysis of the full pedigree imputation, including additional graphs and information not included in the dissertation.

SubPedigreeAnalysis: This outlines the analysis of the sub pedigree after running through AlphaPeel and each model in AlleleFetch (see Section 3.3.1 of the dissertation) before running through the Leave One Out Cross Validation.

SubPedigreeLOOCV: This outlines the Leave One Out Cross Validation used and the resulting analysis for comparing the three models (and verification of AlleleFetch with AlphaPeel).


AlleleFetch: Contains the Function AlleleFetch as described in Section 3.1.3 of the dissertation developed with guidance from Dr Gregor Gorjanc (my supervisor), Dr Mateja Janes (My other supervisor) and Dr Ivan Porcrnc.


