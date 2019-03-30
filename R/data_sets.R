#' Arthropod collection site information
#'
#' A dataset with geologic and geographic information on the colleciton sites
#'
#' @format A data frame with 5 rows and 5 columns:
#' \describe{
#'   \item{code}{collection site ID as it appears in \code{grun}}
#'   \item{name}{full name of the collection site}
#'   \item{age}{substrate age of the collection site in millions of years}
#'   \item{island}{island where the collection site is located}
#'   \item{island_age}{geologic age of the island in millions of years}
#' }
# @source \url{}
"grunSite"


#' Arthropod community abundance and body size data from across the Hawaiian archipelago
#'
#' A dataset containing the community abundace data for individuals, as well 
#' as their body mass.
#'
#' @format A data frame with 15997 rows and 19 columns:
#' \describe{
#'   \item{SpecimenCode}{specimen ID}
#'   \item{XID}{foo}
#'   \item{CollectionCode}{collection ID}
#'   \item{Site}{collection site}
#'   \item{Tree}{collection tree ID}
#'   \item{Tray}{collection tray ID}
#'   \item{DateCollected}{collection date}
#'   \item{Class}{taxonomic class}
#'   \item{Order}{taxonomic class}
#'   \item{Family}{taxonomic class}
#'   \item{SpeciesCode}{unique species ID}
#'   \item{Genus}{taxonomic genus}
#'   \item{species}{species epithet}
#'   \item{trophic}{trophic level: `D` = detritivore, `H` = herbivore, `P` = predator, `T` = transient, `U` = unknown}
#'   \item{Origin}{geographic origin}
#'   \item{Stage}{life stage}
#'   \item{Lengthmm}{length in mm}
#'   \item{Abundance}{abundance in the sample}
#'   \item{IND_BIOM}{estimated individual biomass in mg}
#' }
#' @source Gruner, D. S. 2007. Geological age, ecosystem development, and local resource constraints on arthropod community structure in the Hawaiian Islands. Biological Journal of the Linnean Society, 90: 551--570.
"grun"
