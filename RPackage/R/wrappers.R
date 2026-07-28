#' Generate dataframe and optional CSV of colours for requested taxa (global mode or dataset-specific mode)
#'
#' @param database_name (string) the name of the database you would like all phyla for (ncbi or gtdb)
#' @param database_version (string or list) string: the release version, or empty for the latest version, list: files for your database (see https://github.com/pirovc/multitax for more information)
#' @param requested_taxa (string or DataFrame): path to a CSV containing column 'Taxon' that is the scientific name of all taxa you need a colour for OR
#'           a dataframe of the same format
#' @param name_save_colours (string): name to save requested colours (see OUTPUT) as CSV. If do not want to save CSV, use ""
#' @param name_provided_colours (string): OPTIONAL: A path to a CSV of parent colours (at least phyla needed) in format of 2 columns:
#'           - 'Taxon' contains scientific name of phyla (and additional parent colours at lower taxonomic levels if requested)
#'           - 'Colour' contains hex code (with leading "#") of the corresponding taxon
#' @param mode (string): OPTIONAL: can be global (always generate the same colour for the same taxon) or "dataset-specific" (generate colours such that all
#'           taxa in the same phyla are maximally different). Default is "global"
#' @return A dataframe of 2 columns: 1) 'Taxon' with values from name_requested_colours and 2)'Colour' the corresponding hex colour code
#' @export
getRequestedColours <- function(database_name, database_version, requested_taxa, name_save_colours, name_provided_colours = "", mode = "global") {
  PhylaChromeMain_mod$getRequestedColours(database_name = database_name, database_version = database_version, requested_taxa = requested_taxa, name_save_colours = name_save_colours, name_provided_colours = name_provided_colours, mode = mode)
}

#' Extract CSV of default phyla colours for ncbi or gtdb
#'
#' @param database_name (string) the name of the database you would like default phyla colours for (ncbi or gtdb)
#' @param save_name (string) the full path to where the CSV should be saved
#' @return Nothing returned, but a CSV is saved
#' @export
getDefaultColours <- function(database_name, save_name) {
  PhylaChromeMain_mod$getDefaultColours(database_name = database_name, save_name = save_name)
}

#' Extract CSV of all phyla names and an empty column for Colour (ncbi or gtdb)
#'
#' @param database_name (string) the name of the database you would like phyla from (ncbi or gtdb)
#' @param database_version (string or list) string: the release version, or empty for the latest version, list: files for your database (see https://github.com/pirovc/multitax for more information)
#' @param save_name (string) the full path to where the CSV should be saved
#' @return Nothing returned, but a CSV is saved
#' @export
getPhylaCSV <- function(database_name, database_version, save_name) {
  PhylaChromeMain_mod$getPhylaCSV(database_name = database_name, database_version = database_version, save_name = save_name)
}

#' Extract CSV of all phyla names and an empty column for Colour (ncbi or gtdb)
#'
#' @param database_name (string) the name of the database you would like bacterial phyla from (ncbi or gtdb)
#' @param database_version (string or list) string: the release version, or empty for the latest version, list: files for your database (see https://github.com/pirovc/multitax for more information)
#' @param save_name (string) the full path to where the CSV should be saved
#' @return Nothing returned, but a CSV is saved
#' @export
getPhylaCSVBacteria <- function(database_name, database_version, save_name) {
  PhylaChromeMain_mod$getPhylaCSVBacteria(database_name = database_name, database_version = database_version, save_name = save_name)
}

