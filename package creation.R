library(roxygen2)
library(devtools)
roxygenise()
devtools::load_all(".")
getwd()
devtools::create("DNAfusion")
devtools::document()
usethis::use_package("GenomicRanges")
usethis::use_package("Rsamtools")
usethis::use_package("dplyr")
usethis::use_package("bamsignals")
usethis::use_package("IRanges")
rlang::last_error()
devtools::load_all()
ALK_sequence
?EML4_ALK_detection()
?break_position()
?break_position_depth()
?EML4_sequence
?ALK_sequence
?EML4_ALK_analysis()

?devtools::create()

#https://kbroman.org/github_tutorial/
#https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html
#https://ourcodingclub.github.io/tutorials/writing-r-package/
#https://docs.github.com/en/authentication/connecting-to-github-with-ssh
library(devtools)

getwd()
library(DNAfusion)
H3122_bam <- system.file("extdata", "H3122_EML4.bam", package = "ALKfusiondiscover")
HCC827_bam <-  system.file("extdata", "HCC827_EML4.bam", package = "ALKfusiondiscover")
head(EML4_ALK_detection(H3122_bam))
EML4_ALK_detection(HCC827_bam)
break_position(EML4_ALK_detection(H3122_bam))
break_position_depth(H3122_bam, EML4_ALK_detection(H3122_bam))
EML4_sequence(EML4_ALK_detection(H3122_bam))
ALK_sequence(EML4_ALK_detection(H3122_bam))
.rs.restartR()
getwd()
setwd("C:/Users/Christoffer/OneDrive/1PhD/ALK patienter/EML4_ALK-detection/DNAfusion")

library(GenomicRanges)
library(bamsignals)
?bamCoverage
?EML4_ALK_detection
?EML4_sequence
class("hg38")
tinytex:::install_prebuilt()
C
Sys.getlocale()
?Sys.time
?format

format(Sys.time(), "%a %b %d %X %Y")
library(devtools)
install_github("CTrierMaansson/DNAfusion", build_vignettes = TRUE)
browseVignettes("DNAfusion")
Sperson
BiocManager::install("DNAfusion")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::valid() 
BiocManager::install("GenomicAlignments", update = TRUE, ask = FALSE)
library(DNAfusion)
?EML4_ALK_detection()
devtools::check("DNAfusion")
getwd()
?inherits

colnames(H3122_res)
H3122_res[,1]

H3122_bam <- system.file("extdata", "H3122_EML4.bam", package = "DNAfusion")
HCC827_bam <- system.file("extdata", "HCC827_EML4.bam", package = "DNAfusion")
H3122_res <- EML4_ALK_detection(H3122_bam)
HCC827_res <- EML4_ALK_detection(HCC827_bam)
EML4_ALK_detection(file = H3122_bam, genome = "hg38", mates = 2)
EML4_sequence(H3122_res)
break_position(H3122_res)
break_position_depth(H3122_bam, H3122_res)
?break_position_depth

BiocCheck::BiocCheck('DNAfusion'=TRUE)
?BiocCheck
BiocManager::install("BiocCheck")
devtools::document()
library(BiocCheck)

BiocCheck::BiocCheckGitClone("DNAfusion")
install.packages("goodpractice")
library(goodpractice)
?gp
g <- gp("DNAfusion")
devtools::build("DNAfusion")


?build
?devtools::check()
g$path
g$package
g$extra_preps
g$extra_checks
g$covr
g$cyclocomp
g$description
g$lintr
g$namespace
g$rcmdcheck
results(g)
g$checks$description_url
library(testthat)
library(RUnit)
?RUnit

BiocManager::install("BiocStyle")
library(BiocStyle)
library(DNAfusion)
install.packages("tools")
texi2dvi

?build_manual
library(devtools)
devtools::install_github("hw538/cfDNAPro", build_vignettes = TRUE)
browseVignettes("cfDNAPro")


#https://contributions.bioconductor.org/readme.html

BiocManager::install("DNAfusion", force = T)
?BiocManager::install
?install.packages
utils::news(package="DNAfusion")
?DNAfusion
library(DNAfusion)
?cfDNAPro
library(cfDNAPro)
?EML4_ALK_analysis()
?warning


devtools::document()

install.packages("covr")
library(covr)
library(DNAfusion)
EML4_ALK_detection(H3122_bam)
EML4_sequence(EML4_ALK_detection(H3122_bam))
package_coverage()
covr::report()
devtools::install_github("CTrierMaansson/DNAfusion")
DNAfusion::ALK_sequence
DNAfusion:::index_helper
#https://contributions.bioconductor.org/git-version-control.html#new-package-workflow
#Brug terminal til at lave de nødvendige commands. 
#Hver gang der laves ændringer, så skal de addes gennem GitHub desktop applikation
#Plus version nummeret skal bumpes én op. Ellers registreres ændringen ikke
#Når de så skal tilføjes til selve github bioconductor siden skrives der:
#git push upstream main:master
#Koden er MyGitHub
#Og for at tilføje det til min egen github:
#git push origin main
#Husk at opdatere NEWS.md

setwd("C:/Users/Christoffer/OneDrive/1PhD/ALK_patienter/EML4_ALK-detection/")
devtools::check("DNAfusion") #Outer
BiocCheck::BiocCheckGitClone("DNAfusion") #Outer
devtools::build_manual("DNAfusion") #Outer
covr::report() #Outer
devtools::build("DNAfusion") #Outer
setwd("C:/Users/Christoffer/OneDrive/1PhD/ALK_patienter/EML4_ALK-detection/DNAfusion")
BiocCheck::BiocCheck('DNAfusion'=TRUE) #Inner
getwd()
library(devtools)
.rs.restartR()
devtools::document("DNAfusion")

