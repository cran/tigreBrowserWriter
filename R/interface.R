# tigreBrowserWriter user-facing API code
# Copyright (C) 2016 Antti Honkela
#
# Portions copied from tigre Bioconductor package
# Copyright (C) 2010 Miika-Petteri Matikainen, Antti Honkela and
# Aalto University
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Create and initialize a database
#'
#' @param dbPath Path to the database file to be created. Empty string
#    creates a temporary database that will be deleted at the end.
#' @param datasetName Name of the data set
#' @param datasetSpecies Optional data set metadata: species
#' @param datasetSource Optional data set metadata: source
#' @param datasetPlatform Optional data set metadata: platform
#' @param datasetDescription Optional data set metadata: description
#' @param datasetSaveLocation Optional data set metadata: save location
#' @param datasetFigureFilename Optional data set metadata: figure file name
#' @return A database object db needed by other tigreBrowserWriter functions
#' @examples
#' \dontrun{
#'   # Create a real database to a file
#'   db <- initializeDb("/path/to/the/database.sqlite", "My Dataset")
#'   closeDb(db)
#' }
#'
#'   # Create a temporary database to be deleted at the end
#'   db <- initializeDb("", "My Dataset")
#'   closeDb(db)
#' @export
initializeDb <- function(dbPath, datasetName, datasetSpecies='', datasetSource='', datasetPlatform='', datasetDescription='', datasetSaveLocation='', datasetFigureFilename='') {
    db <- .openConnection(dbPath)
    .createTables(db)
    datasetId <- .addAndGetDatasetId(db, datasetName, datasetSpecies, datasetSource, datasetPlatform, datasetDescription, datasetSaveLocation, datasetFigureFilename)
    return (list(db=db, datasetId=datasetId, datasetName=datasetName,
                 experimentIds=list(), regulatorIds=list()))
}

#' Insert aliases
#'
#' @param db Database object created by \code{\link{initializeDb}}
#' @param aliasType Name of the alias
#' @param aliases A vector of aliases with names giving the primary identifiers
#' @param aliasSource Optional alias metadata: source
#' @param aliasDescription Optional alias metadata: description
#' @return An updated database object db
#' @examples
#'   db <- initializeDb("", "My Dataset")
#'   aliases <- c("aliasA", "aliasB", "aliasC")
#'   names(aliases) <- c("A", "B", "C")
#'   db <- insertAliases(db, "testalias", aliases)
#'   closeDb(db)
#' @export
insertAliases <- function(db, aliasType, aliases, aliasSource='', aliasDescription='') {
    alias_id <- .addAndGetAliasId(db$db, db$datasetId, aliasType, aliasSource, aliasDescription)

    gene_ids <- .addAndGetProbeGeneIds(db$db, names(aliases))

    dbBegin(db$db)
    for (i in seq_along(gene_ids[[1]])) {
        probe <- gene_ids[i, 1]
        gene_id <- gene_ids[i, 2]
        .addGeneAliases(db$db, alias_id, gene_id, list(probe=aliases[probe]))
    }
    dbCommit(db$db)
    return (db)
}

#' Insert results
#'
#' @param db Database object created by \code{\link{initializeDb}}
#' @param experimentName Experiment name
#' @param regulatorName Regulator name (more detailed experiment identifier)
#' @param figurePath URL path to figures
#' @param loglikelihoods A vector of log-likelihoods of elements identified by names
#' @param baselineloglikelihoods A vector of baseline log-likelihoods of elements identified by names (optional)
#' @param experimentDesc Optional experiment description
#' @param loopVariable Optional: Loop variable (1=regulator, 2=target (default))
#' @param modelTranslation Optional: For gene regulation models, is translation modelled
#' @param numberOfParameters Optional: Number of parameters
#' @param parameterNames Optional: Parameter names
#' @param experimentProducer Optional: Experiment producer (string)
#' @param experimentTimestamp Optional: Experiment timestamp (string)
#' @param parameters Optional: A vector of parameter values for the model
#' @return An updated database object db
#' @examples
#'   db <- initializeDb("", "My Dataset")
#'   logl <- c(-4.0, -2.0, 0.0)
#'   names(logl) <- c("A", "B", "C")
#'   baselogl <- c(1.0, -1.0, 4.0)
#'   names(baselogl) <- names(logl)
#'   db <- insertResults(db, "testexperiment", "testregulator", "",
#'                      logl, baselineloglikelihoods=baselogl)
#'   closeDb(db)
#' @export
insertResults <- function(db, experimentName, regulatorName, figurePath, loglikelihoods, baselineloglikelihoods=NA, experimentDesc='', loopVariable=2, modelTranslation=FALSE, numberOfParameters=NA, parameterNames=NA, experimentProducer='', experimentTimestamp='', parameters=NA) {
    if (experimentDesc == '')
        experimentDesc <- experimentName
    regId <- .addAndGetRegulatorId(db$db, regulatorName, db$datasetId)
    experimentId <- .addAndGetExperimentId(db$db, experimentName, experimentDesc, db$datasetId, regulatorId=regId, loopVariable=loopVariable, modelTranslation=modelTranslation, numberOfParameters=numberOfParameters, parameterNames=parameterNames, producer=experimentProducer, timestamp=experimentTimestamp)
    .addResults(db$db, experimentId, loglikelihoods, baselineloglikelihoods, parameters)
    db$experimentIds[[experimentName]] <- experimentId
    return (db)
}

#' Insert figures
#'
#' @param db Database object created by \code{\link{initializeDb}}
#' @param experimentName Experiment name
#' @param regulatorName Regulator name (more detailed experiment identifier)
#' @param filename URL path to figures. The path can contain the
#'   special form \code{${probe_name}} which will be substituted
#'   by the name of the entity (gene, ...) by the browser.
#' @param name Optional figure name
#' @param description Optional figure description
#' @param priority Integer priority used for sorting figures (default: 0)
#' @return An updated database object db
#' @examples
#'   db <- initializeDb("", "My Dataset")
#'   logl <- c(-4.0, -2.0, 0.0)
#'   names(logl) <- c("A", "B", "C")
#'   baselogl <- c(1.0, -1.0, 4.0)
#'   names(baselogl) <- names(logl)
#'   db <- insertResults(db, "testexperiment", "testregulator", "",
#'                      logl, baselineloglikelihoods=baselogl)
#'   db <- insertFigures(db, "testexperiment", "testregulator",
#'                       "http://foo.invalid/path/${probe_name}_fit.png")
#'   closeDb(db)
#' @export
insertFigures <- function(db, experimentName, regulatorName, filename, name='', description='', priority=0) {
    tryCatch(experimentId <- db$experimentIds[[experimentName]],
             error = function(e) stop("Insert results for experiment before inserting figures."))
    .addFigures(db$db, experimentId, filename=filename, name=name, description=description, priority=priority, figureData=NULL)
    return (db)
}

#' Insert supplementary data
#'
#' @param db Database object created by \code{\link{initializeDb}}
#' @param name Name of the supplementary data
#' @param suppData A vector of supplementary data of elements identified by names
#' @param regulatorName Regulator name the data links to (optional)
#' @param source Optional annotation: source
#' @param platform Optional annotation: platform
#' @param description Optional annotation: description
#' @return An updated database object db
#' @examples
#'   db <- initializeDb("", "My Dataset")
#'   suppdata <- c(1, 2, 3)
#'   names(suppdata) <- c("A", "B", "C")
#'   db <- insertSupplementaryData(db, "supptest", suppdata)
#'   boolsupp <- c(TRUE, TRUE, FALSE)
#'   names(boolsupp) <- names(suppdata)
#'   db <- insertSupplementaryData(db, "supptest_bool", boolsupp)
#'   closeDb(db)
#' @export
insertSupplementaryData <- function(db, name, suppData, regulatorName=NA, source='', platform='', description='') {
    if (!is.na(regulatorName)) {
        regulatorId <- .addAndGetRegulatorId(db$db, regulatorName, db$datasetId)
    } else {
        regulatorId <- NA
    }
    if (description == '')
        description <- name
    if (is.logical(suppData)) {
        type = 0
    } else if (is.factor(suppData)) {
        type = 1
    } else {
        type = 2
    }
    suppDatasetId <- .addAndGetSupplementaryDataId(db$db, name, regulatorId, type, source, platform, description)
    .addSupplementaryData(db$db, suppDatasetId, suppData)    
    return (db)
}

#' Insert data z-scores used for filtering
#'
#' @param db Database object created by \code{\link{initializeDb}}
#' @param zscores A vector of z-scores of elements identified by names
#' @return An updated database object db
#' @examples
#'   db <- initializeDb("", "My Dataset")
#'   zscores <- c(1, 2, 3)
#'   names(zscores) <- c("A", "B", "C")
#'   db <- insertZScores(db, zscores)
#'   closeDb(db)
#' @export
insertZScores <- function(db, zscores) {
    .addZscores(db$db, db$datasetId, zscores)
    return (db)
}

#' Finalise and close the database
#'
#' @param db Database object created by \code{\link{initializeDb}}
#' @param experimentSet Name of the experiment set for all the experiments (optional)
#' @examples
#'   db <- initializeDb("", "My Dataset")
#'   # ...
#'   closeDb(db)
#' @export
closeDb <- function(db, experimentSet='') {
    rootId <- .addAndGetExperimentSetId(db$db, 'All experiments', NA)
    if (experimentSet != '')
        setId <- .addAndGetExperimentSetId(db$db, paste(experimentSet, ' (', db$datasetName, ')', sep=''), rootId)
    else
        setId <- rootId
    for (i in seq_along(db$experimentIds)) {
        .addExperimentSetExperiments(db$db, setId, db$experimentIds[[i]])
    }
    dbDisconnect(db$db)
    return (list())
}
