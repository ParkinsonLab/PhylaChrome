# PhylaChrome: A systematic and deterministic tool for assigning colours that capture taxonomic relationships in microbiome datasets 

Please cite: https://www.biorxiv.org/content/10.64898/2026.07.29.741550v1

## Highlights

- From a set of phylum-level colours, assign colours to lower-level bacterial taxa within minutes.
- Use default phylum-level colours, or supply your own!
- multitax was used to support custom and existing NCBI and GTDB versions (https://github.com/pirovc/multitax#local)

## Install
### Required Python packages on your device:
- python >=3.8 (tested with version 3.14.4)
- numpy <code>pip install numpy</code>
- attrs <code>pip install attrs</code>
- chromato <code>pip install chromato</code>
- pandas <code>pip install pandas</code>
- multitax <code>pip install multitax</code>
- script <code>pip install scipy</code>

### Use R package:
#### Pull from GitHub
```
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("ParkinsonLab/PhylaChrome", subdir = "RPackage")
```

#### Using Zip file:
1. Download Zip file PhylaChrome_1.0.0.tar.gz
2. Install the package: <code>remotes::install_local("your_path/PhylaChrome_1.0.0.tar", build = FALSE)</code>

### Use as Python Script:
- Download 3 files from PythonScripts and place each file **in your working directory** (or edit the paths for default CSVs in PhylaChromeMain.py - not recommended)
- To use functions in own script, call <code>from PhylaChromeMain import * </code>

## How-To Use 
### <code>getRequestedColours</code>
The main function is <code>getRequestedColours</code>. This function generates taxonomically-informed colours for bacterial taxa from phylum-level to species-level. All taxa of a phylum are assigned the same hue, lower-level taxa are differenitated by changes to saturation and lightness values. Default bacterial phylum-level colours are provided for NCBI and GTDB databases as of May 2026. 

#### Input Arguments
1. **database name** (string): the database your taxa were generated from - <code>ncbi</code> or <code>gtdb</code>
2. **database version** (string or list):

   string: the release version of the database your taxa were generated from. For the most recent supported version, use <code>""</code>. For NCBI, only the most recent version is available. See https://github.com/pirovc/multitax for a complete list of supported versions.

   list: For any release version or a custom NCBI or GTDB database. A list of paths (as strings) to files for the database 
         (ex. ["nodes.dmp", "names.dmp"] or ["tarDump.tar.gz"]). See https://github.com/pirovc/multitax for more information.
3. **path to requested taxa** (string): path to a CSV containing column 'Taxon' that contains the scientific names of taxa you want colours assigned to. One row per taxon.
4. **path to save generated requested taxa colours** (string): name to save requested colours (see OUTPUT) as CSV. If do not want to save a CSV, use <code>""</code>
5. **OPTIONAL: path to parent colours file** (string):  A path to a CSV of parent colours (at least phyla needed) in format of 2 columns:
    - 'Taxon' contains scientific name of phyla (and additional parent colours at lower taxonomic levels if requested which will be applied to all children of the included taxon)
    - 'Colour' contains hex code (with leading "#") of the corresponding taxon
6. **OPTIONAL mode** (string): The mode you would like used to generate colours for requested taxa.
    - global: the entire bacterial taxonomic tree is considered. For a given database and version, the same colour will always be assigned to the same taxon.
    - dataset-specific: within each phylum hue is constant but taxa are maximally separated by saturation and lightness values. 

 
#### Output
1. **generated requested taxa colours**  (dataframe): a dataframe of 2 columns 
      - 'Taxon' same as Taxon column from input requested taxa
      - 'Colour' assigned colour (value is "N/A" if the taxon could not be located in the requested database)
2. If a path to save the generated colours is provided, the data frame will also be saved as a CSV

### Example Function Calls (check these)
<code>getRequestedColours(database name, database version, path to requested taxa, path to saved colours, path to parent colours, mode)</code>

| Description  | Function Call |
| ------------- | ------------- |
| Call GTDB release 226, use default phyla colours, save CSV, default global mode  |  <code>colours = getRequestedColours("gtdb", "226", "requestedGTDB.csv", "savedColours.csv")</code> |
| Call GTDB release 226, use personal parent colours, save CSV, default global mode |  <code>colours = getRequestedColours("gtdb", "226", "requestedGTDB.csv", "savedColours.csv", "presonalParentColours.csv")</code> |
| Use latest NCBI version, use default phyla colours, do not save CSV, default global mode  |  <code>colours = getRequestedColours("ncbi", "", "requestedNCBI.csv", "")</code>  |
| Use latest NCBI version, use default phyla colours, do not save CSV, use default parent colours, dataset-specific mode  |  <code>colours = getRequestedColours("ncbi", "", "requestedNCBI.csv", "", "", "dataset-specific")</code>  |
| Use custom GTDB version, save CSV, use custom phyla colours, dataset-specific mode  |  <code>colours = getRequestedColours("gtdb", ["nodes.dmp", "names.dmp"], "requestedGTDB.csv", "outputGTDB.csv", "updatedDefault.csv", "dataset-specific")
</code>  |

## Additional Functions
### <code>getDefaultColours</code>
Extract NCBI or GTDB default phylum-level colours as a CSV. Files can also be downloaded from this Github page (phylaNCBIMay2026, and phylaGTDB226.csv)

#### Arguments
1. **database name** (string): the database your taxa were generated from ("ncbi" or "gtdb")
2. **path to save CSV** (string): full path and name for the extracted CSV
   
#### Output
- Saved CSV with columns 'Taxon' and 'Colour' for the default phyla-level colours of the requested database.

### <code>getPhylaCSV</code>
For a given database, create and save a CSV of all phyla and an empty column of colours. This can be used when defining your own phylum-level colours for an NCBI or GTDB database. Any supported NCBI or GTDB version can be requested.

#### Arguments
1. **database name** (string): the database you your taxa were generated from ("ncbi" or "gtdb")
2. **database version** (string): the release version of the database your taxa were generated from. For the latest database version, an empty string can be supplied (<code>""</code>). See https://github.com/pirovc/multitax#local for a list of all supported versions
3. **path to save CSV** (string): full path and name for the created CSV
   
#### Output
- Saved CSV with columns 'Phylum' and 'Colour', where Colour is empty

### <code>getPhylaCSVBacteria</code>
Same as getPhylaCSV but only phyla in domain Bacteria are included.

## Tips and Tricks
- For requested colours the taxa names can be used with or without taxonomy prefixes (ex. "p__Actinomycetota" or "Actinomycetota")
- If creating your own default colours, exact scientific names as they appear in the datatbase are expected (no prefix).
- In R generate a colour palette called <code>colourMap</code> that can be used with ggplot:
  ```
  # Read CSV (or use dataframe output from getRequestedColours)
  savedColoursCSV <- read.csv("savedColours.csv")
  
  #extract the colours (Hex values)
  colourMap <- savedColoursCSV$Colour
  
  #Assign taxa names to corresponding colours for plotting
  names(colourMap) <- savedColoursCSV$Taxon
  ```
- If there is a lower-level taxon you want to assign a specific colour to, add the taxon and its colour to the parent colour CSV. Example R code for NCBI taxonomy is:
  
  ```
  # Save a copy of default phyla CSV (or create your own) - see section Additional Functions for more information
  getDefaultColours("ncbi", "parentTaxaColours.csv"):

  # Read the CSV
  parentTaxaColours = read.csv("parentTaxaColours.csv")

  # Create a new row defining class Blastocatellia as bright pink 
  newRow <- data.frame(Taxon = "Blastocatellia" , Colour = "#FF007F")
  parentTaxaColours <- rbind(parentTaxaColours, newRow)

  # Save the updated CSV
  write.csv(parentTaxaColours, "parentTaxaColoursNew.csv")
  ```
  You can now provide "parentTaxaColoursNew.csv" as the phylum-level colours file in <code>getRequestedColours</code>. All children of the added taxon will be generated based on the new colour instead of the parent phylum colour.

  Note: The program converts hex codes to hue saturation and lightness (HSL) for calculations. If the added colour has the same hue as another parent colour, the assigned colours of all taxa with the same hue may change.

- Default parent phyla colours were tested using the NCBI database as of May 2026 and GTDB release 226. If your taxa are generated from a database version with different phyla, you may need to create your own phylum-level colours. A template CSV for any supported database and version can be produced for you using <code>getPhylaCSVBacteria</code> or <code>getPhylaCSV</code> (see Addition Functions section for more details). Where applicable, phylum-level colours can be transferred from the default phylum-level colours file that can be extracted using <code>getDefaultColours</code>.
- Use of TaxIDs in place of scientific names is currently not supported. We recommend using taxize (R) to convert taxIDS to scientific names (see https://docs.ropensci.org/taxize/articles/taxize.html for examples).

