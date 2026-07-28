# Purpose: This script contains the helper functions for the microbiome colours package PhylaChrome
from chromato.spaces import * # This contains the functions to convert between colour spaces
import hashlib
import pandas as pd
from multitax import GtdbTx, NcbiTx
import numpy as np
import warnings
from typing import Union
from scipy.optimize import minimize, Bounds, linear_sum_assignment
from scipy.spatial.distance import pdist, squareform, cdist
import math
import matplotlib.pyplot as plt


min_S = 0.2 # min saturation allowed
max_S = 1 # max saturationa allowed
min_L = 0.3 # min lightness allowed
max_L = 0.8 # may lightness allowed
min_range = 0.05 # given parent P, S and L must have range P+/- min_range (change back if bad --> test in R)

def pairwise_distance_objective(points_1D: np.array):
    # PURPOSE: calculate the summed pairwise distance between points
    # INPUT:
    #   points_1D (np.array): a 1D array ex. [x1, y1, x2, y2, ..., xn, yn]
    # calculate pairwise distance between all points
    # OUTPUT:
    #   final_val (float): the sum of squared Euclidean distances

    # convert flattened array back to 2D array (N x 2 columns)
    penalty = 0

    # add large penalty if points are the same

    points_2D = points_1D.reshape(-1, 2)

    # find if 2 points are the same, if so add a penalty term

    # calculate pair-wise distances between points
    dist_matrix = pdist(points_2D, metric='euclidean')

    # if any 2 points are within 0.01 distance of each other, then add a penalty 
    if sum(dist_matrix < 0.01) > 1:
        penalty = 10

    # if any points are the same, make this value inf
    if min(dist_matrix) == 0:
        final_val = math.inf
    else:
        final_val = 1/min(dist_matrix)

    # return the sum of squared distances
    return final_val**2 + penalty

def optimizeDistance(num_children: int):
    # PURPOSE: for a given number of  datapoints, maximize the distance in a 2D space (here maxmimize distance of children in S/L space)
    # INPUT:
    #   num_children(int):number of children to maximize distance between points
    # OUTPUT:
    #   optimized_points (np.array): 2D array of size num_children x 2 which are L and S coordinates maximally spread apart in S and L space
    global min_S
    global max_S
    global min_L
    global max_L
    global min_range 

    N = num_children
    # repeat boundaries so appliesd to every data point
    lower_bounds = np.tile([min_S, min_L], N)
    upper_bounds = np.tile([max_S, max_L], N)
    all_bounds = Bounds(lower_bounds, upper_bounds)

    # Randomly place N points within the boundaries
    np.random.seed(42)
    initial_points = np.random.uniform([min_S, min_L], [max_S, max_L], size=(N, 2)).flatten() # [[x1, y1], [x2, y2]] to [x1, y1, x2, y2]

    # Run the optimization algorithm
    optimized_points = minimize(
        fun=pairwise_distance_objective,
        x0=initial_points,
        method='SLSQP', 
        bounds=all_bounds
    )

    # matlab plot of the pointsL
   
    # reshape the points back to 2D
    return optimized_points.x.reshape(-1, 2)

def generateJaccardDistanceMatrix(tax, nodes: list[str], lineages: list[set]):
    # PURPOSE: for a set of nodes, calculate the Jaccard similarity matrix of lineages
    # INPUT:
    #   nodes(list[str]): a list of nodes (as strings)
    #   lineages(list[set]): lineage for each corresponding node in nodes that has been converted to a set
    # OUTPUT:
    #   node_distances(np.array): an array of num_nodes x num_nodes where entry i,j is the Jaccard similarity of lineage for nodes[i] and nodes[j]

    num_nodes = len(nodes)

    # create a 2D array of 0s of size num_nodes x num_nodes
    node_distances = np.zeros((num_nodes, num_nodes))

    # for every node pairing, calcuklate Jaccard distance of lineage and add to node_distances
    for node_A in range(0, num_nodes):
        for node_B in range(node_A + 1, num_nodes):
            # try only including lineages if of rank phylum, class, order, family, genus, or species
            
            lineage_A = set([t for t in lineages[node_A] if tax.rank(t) == 'phylum' or tax.rank(t) == 'class' or tax.rank(t) == 'order' or tax.rank(t) == 'family' or tax.rank(t) == 'genus' or tax.rank(t) == 'species'])
            lineage_B = set([t for t in lineages[node_B] if tax.rank(t) == 'phylum' or tax.rank(t) == 'class' or tax.rank(t) == 'order' or tax.rank(t) == 'family' or tax.rank(t) == 'genus' or tax.rank(t) == 'species'])
            intersection_lineages = lineage_A.intersection(lineage_B)
            union_lineages = lineage_A.union(lineage_B)
            similarity = len(intersection_lineages) / len(union_lineages)
            node_distances[node_A][node_B] = 1/similarity # for difference
            node_distances[node_B][node_A] = 1/similarity # for difference

    # for every pair of nodes, calculate the jaccard similarity score and place it in matrix

    return node_distances

def optimizeMatrixPairing(lineages_distance: np.array, points_distance: np.array):
    # PURPOSE: for a set of maximally separated points in the L/S space, assign taxa to S/L values such that lineage differences are respected
    # INPUT: 
    #   lineages_distance(np.array): a numpy array of size NxN where entry i,j is the Jaccard similarity between i and j
    #   points_distance(np.array): a numpy array of size NxN, N being the number of taxa, where entry i,j is the distance between SL point i and j
    # OUTPUT:
    #   point_assignment(np.array): a numpy array of size N where lineage j was assigned to point number at entry j

    num_nodes = len(points_distance) # returns number of rows

    # calculate cost matrix
    cost_matrix = np.zeros((num_nodes, num_nodes))


    cost_matrix = cdist(lineages_distance, points_distance, metric='seuclidean')

    # perform distacne profile matching using linear sum assignment
    lineage_assignment, point_assignment = linear_sum_assignment(1/cost_matrix)

    return point_assignment


def getOpenSLChannels(parent_HSL: list[float], hue_range: float, dict_HSL: dict):
    # PURPOSE: For a given parent HSL values, check saturation and lightness values available for the child.
    # INPUT: 
    #   parent_HSL (list[float): [hue, saturation, lightness] of parent node (assuming using scales from chromato.spaces)
    #       hue_range (float): the parent hue +/- the hue_range will be checked for current S and L values in use
    #       dict_HSL (dict): containing all HSL values (nodes are hue, values are 2D array where [0] are S values and [1] are L values )
    # OUTPUT: [min_curr_S, max_curr_S, min_curr_L, max_curr_L] where these are minimum and maximum saturation and lightness 
    #          values that haven't been used yet and can be used for the children of the parent node

    # define needed global
    curr_S = []
    curr_L = []

    global min_S
    global max_S
    global min_L
    global max_L
    global min_range 

    min_curr_S = min_S
    max_curr_S = max_S
    min_curr_L = min_L
    max_curr_L = max_L

    # extract saturation and lightness values from neighboring hues (do not)
    hue_check_min = parent_HSL[0] - hue_range
    hue_check_max = parent_HSL[0] + hue_range
    
    hues_to_check = [k for k in dict_HSL if k > hue_check_min and k < hue_check_max] 

    for h in hues_to_check:
            # if the hue has been used, extract saturation and lightness values
            curr_S = curr_S + dict_HSL.get(h)[0] # S values
            curr_L = curr_L + dict_HSL.get(h)[1] # L values
        
    # sort the lists in ascending order
    curr_S.sort()
    curr_L.sort()

    # Find the value closest to parent S that is smaller and the value closest that is larger
    # split the lists in 2 (< hue S and > hue S)
    less_than_S = [val for val in curr_S if val < parent_HSL[1]]
    greater_than_S = [val for val in curr_S if val > parent_HSL[1]]

    # if there is a S less than parent S, take point halfway between
    if less_than_S:
        min_curr_S = (parent_HSL[1] + less_than_S[-1])/2 # get value halfway
        min_curr_S = min(parent_HSL[1] - min_range, min_curr_S)
    
    # if there is a S greater than parent S, take point halfway between
    if greater_than_S:
        max_curr_S = (greater_than_S[0] + parent_HSL[1])/2 # get value halfway
        max_curr_S = max(parent_HSL[1] + min_range, max_curr_S)


    # Find the value closest to parent L that is smaller and the value closest that is larger
    # split the lists in 2 (< hue L and > hue L)
    less_than_L = [val for val in curr_L if val < parent_HSL[2]]
    greater_than_L = [val for val in curr_L if val > parent_HSL[2]]

    # if there is a L less than parent L, take point halfway between
    if less_than_L:
        min_curr_L = (parent_HSL[2] + less_than_L[-1])/2 # get value halfway
        min_curr_L = min(parent_HSL[2] - min_range, min_curr_L)

    # if there is a L greater than parent L, take point halfway between
    if greater_than_L:
        max_curr_L = (greater_than_L[0] + parent_HSL[2])/2 # get value halfway
        max_curr_L = max(parent_HSL[2] + min_range, max_curr_L)

    # check the min boundaries are enforced. 
    if min_curr_S < min_S:
        diff_min_S = min_S - min_curr_S
        # increase min_curr_S to min_S and add the difference to the max side (if not already at max)
        min_curr_S = min_S
        max_curr_S = min(max_curr_S + diff_min_S, max_S)

    if min_curr_L < min_L:
        diff_min_L = min_L - min_curr_L
        # increase min_curr_L to min_L and add the difference to the max side (if not already at max)
        min_curr_L = min_L
        max_curr_L = min(max_curr_L + diff_min_L, max_L)    


    # check that max boundaries are enforced
    if max_curr_S > max_S:
        diff_max_S = max_curr_S - max_S
        # decrease max_curr_S to max_S and subtract the difference from min side (if not already min)
        max_curr_S = max_S
        min_curr_S = max(min_curr_S - diff_max_S, min_S)

    if max_curr_L > max_L:
        diff_max_L = max_curr_L - max_L
        # decrease max_curr_L to max_L and subtract the difference from min side (if not already min)
        max_curr_L = max_L
        min_curr_L = max(min_curr_L - diff_max_L, min_L)

    # return the ranges
    return min_curr_S, max_curr_S, min_curr_L, max_curr_L

def getChildNodecolourHSL(child_name: str, range_S: list[float], range_L: list[float]):
    # Purpose: Given the allowable saturation and lightness range, calculate S and L values for the child
    # INPUT:
    #       child_name (string): scientific name of node
    #       range_S (list[floats]): minimum and maximum S given parent S and currently in-use S values in format [min_S, max_S]
    #       range_L (list[floats]): minimum and maximum L given parent L and currently in-use L values in format [min_L, max_L]
    # OUTPUT:
    #       
    # Convert name to 4 digit integer used sha256 hash function 
    hash_4_dig = int(hashlib.sha256(child_name.encode('utf-8')).hexdigest(), 16) % 10**4

    hash_4_dig_str = f"{hash_4_dig:04d}"

    # first 2 digits of integer = location on scale between min and max saturation
    digs_S = (int(hash_4_dig_str[0:2]))
    # second 2 digits of integer = location on scale between min and max saturation
    digs_L = (int(hash_4_dig_str[2:4]))

    # Convert digits to a percentage and calculate that percentage of the allowable saturation, repeat for lightness.
    # Ex. min is 25, max is 45, digit is 50, the returned value will be 35 (half way between the min and max)
    val_S = ((range_S[1] - range_S[0]) * digs_S/100) + range_S[0]
    val_L = ((range_L[1] - range_L[0]) * digs_L/100) + range_L[0]

    return val_S, val_L

def getPhylumcolour(name: str, provided_colours: pd.Dataframe):
    # PURPOSE: from the provided phyla level colours, extract the hex value for the phylum name
    # INPUT
    #   name (str): the exact name of the phylum as it appears in the CSV
    #   provided_colours (pandas Dataframe): dataframe with columns
    #       - 'Phylum' which contains the the strring phyla names
    #       - 'Colour' which are the hex codes for the corresponding phyla
    # OUTPUT:
    #       hex_val (str): provided hex value for requested phylum or colour for 'Other' if

    # find if the colour is in the list of parent colours
    hex_val = provided_colours.loc[provided_colours['Taxon'].str.contains(name), 'Colour']

    # if not replace it with the colour for Other
    if hex_val.empty:
        hex_val = provided_colours.loc[provided_colours['Taxon'].str.contains('Other'), 'Colour'].values[0]
    
        # if no colour for Other provided, use default light gray
        if hex_val.empty:
            hex_val = "030303"

        # if Other colour was provided, remove the leading "#"
        else:
            hex_val = hex_val[1:]

    # if a colour was provided, extract the colour
    else:
        hex_val = hex_val.values[0]
        # remove the leading "#"
        hex_val = hex_val[1:]

    return hex_val

def getPhylaCSV(database_name: str, database_version: Union[str, list], save_name: str):
    # PURPOSE: for a given database, create and save a CSV of all phyla and an empty column of Colours. This can be used
    #        when defining your own phyla-level colours for an NCBI or GTDB database.
    # INPUT:
    #       database_name (string): the name of the database you would like all phyla for (ncbi or gtdb)
    #       database_version (string or list): 
    #           string: the release version, or empty for the latest version 
    #           list: for a custom database.  A list of paths (as strings) to files for database 
    #           (ex. ["nodes.dmp", "names.dmp"] or ["tarDump.tar.gz"]) - See https://github.com/pirovc/multitax for more information.
    #       save_name (string): full path for saving CSV
    # OUTPUT: none
    # SAVED RESULTS: CSV with columns 'Phylum' and 'Colour', where Colour is empty

    # load the requested database (default is most recent)
    if database_name == "ncbi":
        # check if the latest version was requested
        if isinstance(database_version, str):
            # if empty string, load latest version
            if database_version == "":
                print("--- Loading NCBI database, latest version available (see https://github.com/pirovc/multitax#local) ---")
                tax = NcbiTx(version=database_version)
            else:
                print("--- Loading NCBI version " + database_version + " ---")
                tax = NcbiTx(version=database_version)

        # check if a custom database was provided
        elif isinstance(database_version, list):
            if not database_version:
                print("ERROR: Please provide pathways to database files")
                return 1
            else:
                print("--- Loading custom NCBI database ---")
                tax = NcbiTx(files = database_version)
        else:
            print("ERROR: database version must be a string, empty string or list")
            return 1

        
    elif database_name == "gtdb":
        # check if the latest version was requested
        if isinstance(database_version, str):
            # if empty string, load latest version
            if database_version == "":
                print("--- Loading GTDB database, latest version available (see https://github.com/pirovc/multitax#local) ---")
                tax = GtdbTx(version=database_version)
            else:
                print("--- Loading GTDB r" + database_version + " ---")
                tax = GtdbTx(version=database_version)

        # check if a custom database was provided
        elif isinstance(database_version, list):
            if not database_version:
                print("ERROR: Please provide pathways to database files")
                return 1
            else:
                print("--- Loading custom GTDB database ---")
                tax = GtdbTx(files = database_version)
        else:
            print("ERROR: database version must be a string, empty string or list")
            return 1
        
    else:
        print("Error: Requested database " + database_name + " version " + database_version + " is not available")
        return 1


    # get list of all phyla nodes
    phyla_arr = [t for t in list(tax.__dict__.get('_nodes')) if tax.rank(t) == 'phylum']
    # get scientific names of all phyla
    phyla_arr = [tax.name(n) for n in phyla_arr]

    # save the CSV
    np.savetxt(save_name, phyla_arr, delimiter=",", fmt='%s', header="Name, Colour")

def getPhylaCSVBacteria(database_name: str, database_version: Union[str, list], save_name: str):
    # PURPOSE: for a given database, create and save a CSV of all bacterial phyla and an empty column of Colours. This can be used
    #        when defining your own phyla-level colours for an NCBI or GTDB database.
    # INPUT:
    #       database_name (string): the name of the database you would like all phyla for (ncbi or gtdb)
    #       database_version (string or list): 
    #           string: the release version, or empty for the latest version 
    #           list: for a custom database.  A list of paths (as strings) to files for database 
    #           (ex. ["nodes.dmp", "names.dmp"] or ["tarDump.tar.gz"]) - See https://github.com/pirovc/multitax for more information.    
    #       save_name (string): the full path to where the CSV should be saved 
    # OUTPUT: none
    # OUTPUT: none
    # SAVED RESULTS: CSV with columns 'Taxon' and 'Colour', where Colour is empty

    # load the requested database (default version is most recent)
        # load the requested database (default is most recent)
    if database_name == "ncbi":
        # check if the latest version was requested
        if isinstance(database_version, str):
            # if empty string, load latest version
            if database_version == "":
                print("--- Loading NCBI database, latest version available (see https://github.com/pirovc/multitax#local) ---")
                tax = NcbiTx(version=database_version)
            else:
                print("--- Loading NCBI version " + database_version + " ---")
                tax = NcbiTx(version=database_version)

        # check if a custom database was provided
        elif isinstance(database_version, list):
            if not database_version:
                print("ERROR: Please provide pathways to database files")
                return 1
            else:
                print("--- Loading custom NCBI database ---")
                tax = NcbiTx(files = database_version)
        else:
            print("ERROR: database version must be a string, empty string or list")
            return 1

        
    elif database_name == "gtdb":
        # check if the latest version was requested
        if isinstance(database_version, str):
            # if empty string, load latest version
            if database_version == "":
                print("--- Loading GTDB database, latest version available (see https://github.com/pirovc/multitax#local) ---")
                tax = GtdbTx(version=database_version)
            else:
                print("--- Loading GTDB r" + database_version + " ---")
                tax = GtdbTx(version=database_version)

        # check if a custom database was provided
        elif isinstance(database_version, list):
            if not database_version:
                print("ERROR: Please provide pathways to database files")
                return 1
            else:
                print("--- Loading custom GTDB database ---")
                tax = GtdbTx(files = database_version)
        else:
            print("ERROR: database version must be a string, empty string or list")
            return 1
        
    else:
        print("ERROR: Requested database " + database_name + " version " + database_version + " is not available")
        return 1

    
    # get list of all phyla nodes that are in the domain "Bacteria"
    phyla_arr = [t for t in list(tax.__dict__.get('_nodes')) if (tax.rank(t) == 'phylum' and tax.name(tax.closest_parent(t, ranks = ["domain"])) == 'Bacteria') ]
    # get scientific names of all phyla
    phyla_arr = [tax.name(n) for n in phyla_arr]

    # CHECK HERE ********************************************************
    # if database_name == "gtdb":
    #     phyla_arr = ["p__" + n for n in phyla_arr]

    np.savetxt(save_name, phyla_arr, delimiter=",", fmt='%s', header="Name, Colour")

def getDefaultColours(database_name: str, save_name: str):
    # PURPOSE: for a given database, Extract NCBI or GTDB default phylum-level colours as a CSV. Files can also 
    #       be downloaded from the Github repository (phylaNCBIMay2026, and phylaGTDB226.csv)
    # INPUT:
    #       database_name (string): the name of the database you would like all phyla for (ncbi or gtdb)
    #       save_name (string): the full path to where the CSV should be saved 
    # OUTPUT: none
    # SAVED RESULTS: CSV with columns 'Taxon' and 'Colour' for the default phyla-level colours for requested database

    if database_name == "ncbi":
        provided_colours = pd.read_csv("phylaNCBIMay2026.csv", header='infer')
        provided_colours.to_csv(save_name, index=False)

    elif database_name == "gtdb":
        provided_colours = pd.read_csv("phylaGTDB226.csv", header='infer')
        provided_colours.to_csv(save_name, index=False)

    else:
        print("Error: Requested database " + database_name + " is not available")
        return 1

def generateColoursGlobal(provided_colours:pd.DataFrame, requested_colours:pd.DataFrame, database_name: str, database_version: Union[str, list]):
    # PURPOSE: For a given taxonomic database, generate colours for requested taxa based on phyla colours. For a given taxon, the same colour will always be assigned.
    # INPUT:
    #       provided_colours (pd.DataFrame): A data frame of parent colours (at least phyla needed) in format of 2 columns:
    #           - 'Taxon' contains scientific name of phyla (and additional parent colours at lower taxonomic levels if requested)
    #           - 'Colour' contains hex code (with leading "#") of the corresponding taxon
    #       requested_colours (pd.DataFrame): A dataframe containing column 'Taxon' that is the scientific name of all taxa you need a colour for
    #       database_name (string): the name of the database you would like all phyla for (ncbi or gtdb)
    #       database_version (string or list): 
    #           string: the release version, or empty for the latest version 
    #           list: for a custom database.  A list of paths (as strings) to files for database 
    #           (ex. ["nodes.dmp", "names.dmp"] or ["tarDump.tar.gz"]) - See https://github.com/pirovc/multitax for more information.
    # OUTPUT:
    #       requested_colours (pd.DataFrame): a dataframe of 2 columns 
    #           - 'Taxon' same as Taxon column from requested_colours
    #           - 'Colour' assigned colour


    hex_colour_map = {} # each key will be a node and the value will the hex code string
    used_HSL = {} # each key will be a hue, and the value will be an array of S,L values
    parent_HSL_range = {} # key is a node (parent) and the value is the returned allowable SL ranges --> unless phylum level, uses parent decided range 

    # load the requested database (default is most recent)
    if database_name == "ncbi":
        # check if the latest version was requested
        if isinstance(database_version, str):
            # if empty string, load latest version
            if database_version == "":
                print("--- Loading NCBI database, latest version available (see https://github.com/pirovc/multitax#local) ---")
                tax = NcbiTx(version=database_version)
            else:
                print("--- Loading NCBI version " + database_version + " ---")
                tax = NcbiTx(version=database_version)

        # check if a custom database was provided
        elif isinstance(database_version, list):
            if not database_version:
                print("ERROR: Please provide pathways to database files")
                return 1
            else:
                print("--- Loading custom NCBI database ---")
                tax = NcbiTx(files = database_version)
        else:
            print("ERROR: database version must be a string, empty string or list")
            return 1

        
    elif database_name == "gtdb":
        # check if the latest version was requested
        if isinstance(database_version, str):
            # if empty string, load latest version
            if database_version == "":
                print("--- Loading GTDB database, latest version available (see https://github.com/pirovc/multitax#local) ---")
                tax = GtdbTx(version=database_version)
            else:
                print("--- Loading GTDB r" + database_version + " ---")
                tax = GtdbTx(version=database_version)

        # check if a custom database was provided
        elif isinstance(database_version, list):
            if not database_version:
                print("ERROR: Please provide pathways to database files")
                return 1
            else:
                print("--- Loading custom GTDB database ---")
                tax = GtdbTx(files = database_version)
        else:
            print("ERROR: database version must be a string, empty string or list")
            return 1
        
    else:
        print("ERROR: Requested database " + database_name + " version " + database_version + " is not available")
        return 1

    # get all phyla nodes
    phyla_nodes = [t for t in list(tax.__dict__.get('_nodes')) if tax.rank(t) == 'phylum']

    # only inlcude phyla that are in the provided list of phyla colours
    # these will be our current list of nodes to assign (have to include all in case of overlap)
    nodes_to_assign = [item for item in phyla_nodes if provided_colours['Taxon'].str.contains(tax.name(item)).any()]

    # extract list of requested taxa colours
    colours_needed = requested_colours['Taxon'].tolist()
    
    # get the lowest taxonomic level of requested colours
    colours_needed_nodes = []
    for item in colours_needed:
        
        # get the possible list of scientific names for the requested colour (needed for "p__" type pre-fixes)
        node = [k for k in list(tax.__dict__.get('_names')) if (tax.name(k) in item)]

        # if multiple, then choose the closest one 
        if len(node) > 1:
            min_diff = 1000
            best = node[0]
            # for every node in the list, count the difference in length, and keep track of most similar node
            for sub_node in node:
                diff = abs(len(tax.name(sub_node)) - len(item))
                if diff < min_diff:
                    min_diff = diff
                    best = sub_node
            node = best
        elif node:
            node = node[0]
        if node:
            colours_needed_nodes.append(node)

    colours_needed_ranks = [tax.rank(r) for r in colours_needed_nodes]

    unique_colours_needed_ranks = np.unique(colours_needed_ranks)
    level_stop = ""
    if "species" in unique_colours_needed_ranks:
        level_stop = "species"
    elif "genus" in unique_colours_needed_ranks:
        level_stop = "genus"
    elif "family" in unique_colours_needed_ranks:
        level_stop = "family"
    elif "order" in unique_colours_needed_ranks:
        level_stop = "order"
    elif "class" in unique_colours_needed_ranks:
        level_stop = "class"
    else:
        level_stop = "phylum"

    print("--- Generating colours ---")
    print("Generating colours to level " + level_stop +  ".")
    print("Note: Lower taxonomic levels will take more time (up to a few minutes for 16GB RAM, 10-core CPU)")
    nodes_range_to_assign = [] # in a given iteration will contain list of children from the same parent that need to have ranges assigned
    curr_parent = ""
    rank_reached = []

    # while colours still need to be assigned:
    warnings.filterwarnings("ignore", category=UserWarning, message="This pattern is interpreted as a regular expression")
    while nodes_to_assign:

        # remove the first item from the list
        curr_node = nodes_to_assign.pop(0)

        # If want to be updated on taxonomic level reached for generation
        # if tax.rank(curr_node) not in rank_reached:
        #     print(tax.rank(curr_node))
        #     rank_reached.append(tax.rank(curr_node))
        #     print("Nodes to assign at current level:")
        #     print(len(nodes_to_assign))

        # If the node is not yet at the lowest taxonomic level, add it's children to the list nodes_to_assign
        if level_stop != tax.rank(curr_node):

            # Get all immediate children and append to nodes_to_assign
            nodes_to_assign.extend(tax.children(curr_node))

        # Identify the parent node of the current node, this will be used to 
        parent_node = tax.parent(curr_node)

        # --------------------------------------------------------------------------------------
        # HANDLE ASSIGNING ALLOWABLE RANGES
        # --------------------------------------------------------------------------------------

        # if this is the first parent seen, update curr_parent and nodes that need colours assigned
        if curr_parent == "":
            curr_parent = parent_node
            nodes_range_to_assign.extend([curr_node])

        # check if the parent node is different from the previous one and if so, update the acceptable ranges for all children of the previous parent (all children have been seen)
        elif parent_node != curr_parent and nodes_range_to_assign:
            
            # get the acceptable ranges for all children of the previous parent (preparation for when they are the parent)
            for new_parent in nodes_range_to_assign:
                # get the HSL colour of the new parent node (child of previous parent)  
                new_parent_hex = hex_colour_map.get(new_parent)
                new_parent_HSL = convert.hex_to_hsl(new_parent_hex)
                new_parent_HSL = [float(getattr(new_parent_HSL, 'h')), float(getattr(new_parent_HSL, 's')), float(getattr(new_parent_HSL, 'l'))]
                allowable_SL = getOpenSLChannels(new_parent_HSL, 0.025, used_HSL) # check hues ~+/-10 current hue (visibly the same)
                parent_HSL_range[new_parent] = allowable_SL

            # add the current node as a parent whose range needs to be assigned (it's siblings will be added)
            nodes_range_to_assign = [curr_node]
            curr_parent = parent_node


        # the children of the parent are still be assigned colours, add the current node as a range to be assigned
        else:
            nodes_range_to_assign.extend([curr_node])

        # --------------------------------------------------------------------------------------
        # ASSIGNING COLOR TO CURRENT NODE 
        # --------------------------------------------------------------------------------------
        
        # check if at the lowest requested level and if so, only calculate the colour if it is a requested colour
        calculate_colour = True
        if level_stop == tax.rank(curr_node):

            # if at the lowest level, only calculate colour if it is requested
            if any(tax.name(curr_node) in s for s in colours_needed):
                calculate_colour = True
                
            else:
                calculate_colour = False

        # if the colour needs to be calculated:
        if calculate_colour:

            # if the name is in the provided list of colours, use the provided colour (even if there is an overlap with an existing)
            if provided_colours['Taxon'].str.contains(tax.name(curr_node)).any():
                assigned_hex_colour = getPhylumcolour(tax.name(curr_node), provided_colours)
                hsl_colour = convert.hex_to_hsl(assigned_hex_colour)
                assigned_hsl_colour = [float(getattr(hsl_colour, 'h')), float(getattr(hsl_colour, 's')), float(getattr(hsl_colour, 'l'))]
            
            # otherwise calculate the colour of the current node based on the allowable range from the parent node
            else:
            
                # find the parents assigned colour in HSL
                parent_hex = hex_colour_map.get(tax.parent(curr_node))
                parent_HSL = convert.hex_to_hsl(parent_hex)
                parent_HSL = [float(getattr(parent_HSL, 'h')), float(getattr(parent_HSL, 's')), float(getattr(parent_HSL, 'l'))]

                # Find the allowed range for the current node based on it's parent
                allowable_SL = parent_HSL_range.get(tax.parent(curr_node))

                # get the child nodes colour
                SL_val = getChildNodecolourHSL(tax.name(curr_node), allowable_SL[0:2], allowable_SL[2:4])
                assigned_hsl_colour = [parent_HSL[0], SL_val[0], SL_val[1]]
                assigned_hex_colour = convert.hsl_to_hex(assigned_hsl_colour)
                
            # convert to HSL and store in colour map
            hex_colour_map[curr_node] = str(assigned_hex_colour)
            
            
            # update used HSL ----------------------------------------------------------------------
            if assigned_hsl_colour[0] in used_HSL:
                arrays_SL = used_HSL.get(assigned_hsl_colour[0])
                S_val = arrays_SL[0]
                L_val = arrays_SL[1]
                S_val.append(assigned_hsl_colour[1])
                L_val.append(assigned_hsl_colour[2])
                used_HSL[assigned_hsl_colour[0]] = [S_val, L_val]
            
            else:
                used_HSL[assigned_hsl_colour[0]] = [[assigned_hsl_colour[1]], [assigned_hsl_colour[2]]]

    # --------------------------------------------------------------------------------------
    # EXTRACT REQUESTED COLORS AND SAVE IN DATAFRAME
    # --------------------------------------------------------------------------------------
    requested_colours["Colour"] = ""
    requested_colours = requested_colours.set_index('Taxon')

    print("--- Saving requested colours ---")
    for item in colours_needed:

        if "Other" in item:
            requested_colours.at[item, "Colour"] = "#" + getPhylumcolour("Other", provided_colours)
        
        else:
            # get the possible list of scientific names for the requested colour (needed for "p__" type pre-fixes)
            node = [k for k in list(tax.__dict__.get('_names')) if (tax.name(k) in item)]

            # if multiple, then choose the closest one 
            if len(node) > 1:
                min_diff = 1000
                best = node[0]

                # for every node in the list, count the difference in length, and keep track of most similar node
                for sub_node in node:
                    diff = abs(len(tax.name(sub_node)) - len(item))
                    if diff < min_diff:
                        min_diff = diff
                        best = sub_node

                node = best
    
            elif node:
                node = node[0]
    
            # if a node was found with the name matching the name in requested_colours, then include the colour
            if node:
                if hex_colour_map.get(node):
                    requested_colours.at[item, "Colour"] = "#" + hex_colour_map.get(node)
                # if a colour was not assigned for this node return NA)
    
                else:
                    requested_colours.at[item, "Colour"] = "NA"

            # if a matching name was not found, put NA
            else:
                requested_colours.at[item, "Colour"] = "NA"
        # end if-else for "Other

    requested_colours = requested_colours.reset_index()
    return requested_colours

def generateColoursSpecific(provided_colours:pd.DataFrame, requested_colours:pd.DataFrame, database_name: str, database_version: Union[str, list]):
    # PURPOSE: For a given taxonomic database, generate colours for requested taxa based on phyla colours with maximally differen colours within a phyla for the requested taxa.
    #   All requested taxa will be separated, rank will not be considered.
    #   NOT REPRODUCIBLE WITH DIFFERENT SETS OF TAXA.
    # INPUT:
    #       provided_colours (pd.DataFrame): A data frame of parent colours (at least phyla needed) in format of 2 columns:
    #           - 'Taxon' contains scientific name of phyla (and additional parent colours at lower taxonomic levels if requested)
    #           - 'Colour' contains hex code (with leading "#") of the corresponding taxon
    #       requested_colours (pd.DataFrame): A dataframe containing column 'Taxon' that is the scientific name of all taxa you need a colour for
    #       database_name (string): the name of the database you would like all phyla for (ncbi or gtdb)
    #       database_version (string or list): 
    #           string: the release version, or empty for the latest version 
    #           list: for a custom database.  A list of paths (as strings) to files for database 
    #           (ex. ["nodes.dmp", "names.dmp"] or ["tarDump.tar.gz"]) - See https://github.com/pirovc/multitax for more information.
    # OUTPUT:
    #       requested_colours (pd.DataFrame): a dataframe of 2 columns 
    #           - 'Taxon' same as Taxon column from requested_colours
    #           - 'Colour' assigned colour


    # load the requested database (default is most recent)
    if database_name == "ncbi":
        # check if the latest version was requested
        if isinstance(database_version, str):
            # if empty string, load latest version
            if database_version == "":
                print("--- Loading NCBI database, latest version available (see https://github.com/pirovc/multitax#local) ---")
                tax = NcbiTx(version=database_version)
            else:
                print("--- Loading NCBI version " + database_version + " ---")
                tax = NcbiTx(version=database_version)

        # check if a custom database was provided
        elif isinstance(database_version, list):
            if not database_version:
                print("ERROR: Please provide pathways to database files")
                return 1
            else:
                print("--- Loading custom NCBI database ---")
                tax = NcbiTx(files = database_version)
        else:
            print("ERROR: database version must be a string, empty string or list")
            return 1

        
    elif database_name == "gtdb":
        # check if the latest version was requested
        if isinstance(database_version, str):
            # if empty string, load latest version
            if database_version == "":
                print("--- Loading GTDB database, latest version available (see https://github.com/pirovc/multitax#local) ---")
                tax = GtdbTx(version=database_version)
            else:
                print("--- Loading GTDB r" + database_version + " ---")
                tax = GtdbTx(version=database_version)

        # check if a custom database was provided
        elif isinstance(database_version, list):
            if not database_version:
                print("ERROR: Please provide pathways to database files")
                return 1
            else:
                print("--- Loading custom GTDB database ---")
                tax = GtdbTx(files = database_version)
        else:
            print("ERROR: database version must be a string, empty string or list")
            return 1
        
    else:
        print("ERROR: Requested database " + database_name + " version " + database_version + " is not available")
        return 1
    
    # get requested colours as a list of strings
    colours_needed_str = requested_colours['Taxon'].tolist()
    colours_needed_parent_phyla = []
    colours_needed_node = []
    colours_needed_lineages = []  

    # For every requested colour, find the parent phylum
    for item in colours_needed_str:
        
        # get the possible list of scientific names for the requested colour (needed for "p__" type pre-fixes)
        node = [k for k in list(tax.__dict__.get('_names')) if (tax.name(k) in item)]

        # if multiple, then choose the closest one 
        if len(node) > 1:
            min_diff = 1000
            best = node[0]

            # for every node in the list, count the difference in length, and keep track of most similar node
            for sub_node in node:
                diff = abs(len(tax.name(sub_node)) - len(item))
                if diff < min_diff:
                    min_diff = diff
                    best = sub_node
            node = best
        elif node:
            node = node[0]
        if node:
            # get lineage information for the node
            lineage = tax.lineage(node)
            lineage_names = [tax.name(val) for val in lineage]
            lineage_ranks = [tax.rank(val) for val in lineage]
            colours_needed_node.append(node)

            # append name of parent phylum
            colours_needed_parent_phyla.append(lineage_names[lineage_ranks.index('phylum')])
            colours_needed_lineages.append(set(lineage))

        else:
            colours_needed_parent_phyla.append("None")
            colours_needed_node.append('')
            colours_needed_lineages.append(set(['']))

    # set-up colours dataframe to return
    requested_colours_save = pd.DataFrame()
    requested_colours_save["Taxon"] = colours_needed_str
    requested_colours_save["Colour"] = "NA" # default is NA, replaced if available
    
    # Find any entries containing "Other" and replace
    requested_colours_save.loc[requested_colours_save['Taxon'].str.contains('Other', na=False), 'Colour'] = '#' + getPhylumcolour("Other", provided_colours)


    print("--- Generating and Saving Requested Colours ---")

    # For each parent phylum, select all that are in that phylum and assign colours, save to new CSV
    for parent_name in np.unique(colours_needed_parent_phyla):

        # only calculate phyla for scientific names
        if parent_name != "None":
            # get phylum colour
            phylum_hex_colour = getPhylumcolour(parent_name, provided_colours)
            phylum_hsl_colour = convert.hex_to_hsl(phylum_hex_colour)
            phylum_hue = float(getattr(phylum_hsl_colour, 'h'))
        
            # get all children with the parent phylum
            children_idx = [i for i, x in enumerate(colours_needed_parent_phyla) if x == parent_name]
            children_nodes = [colours_needed_node[idx] for idx in children_idx]
            children_lineages = [colours_needed_lineages[idx] for idx in children_idx]
            
            # get the data points spread out in 2D space (S and L space)
            if len(children_idx) > 1:
                points_mat = optimizeDistance(len(children_idx))

                # get distance between points
                points_distance_condensed = pdist(points_mat, metric='euclidean')
                points_distance= squareform(points_distance_condensed)   

                # get the lineage distances between children
                children_tax_similarity = generateJaccardDistanceMatrix(tax, children_nodes, children_lineages)

                # Get the assignment of data points in SL space to children based on taxonomy similarity
                point_assignment = optimizeMatrixPairing(children_tax_similarity, points_distance)

                # go through all of the datapoints, and if one of them maps to  white, make it colour ['red']

                # PLOT FOR PAPER *************************************************************************
                # colours_ex = []
                # for point_idx in range(0,len(points_mat)):
                #     curr_hex = convert.hsl_to_hex([phylum_hue, points_mat[point_idx][0], points_mat[point_idx][1]])
                #     colours_ex.append('#' + str(curr_hex)) # probable change this to be the hex
                
                # fig, ax = plt.subplots()
                # ax.set_xlim(0, 120)
                # ax.set_ylim(0, 100)     
                # ax.set_xlabel("Saturation")
                # ax.set_ylabel("Lightness")
                # # plot it so we see points
                # x = [20, 100, 100, 20]
                # y = [80, 80, 30, 30]
                # ax.fill(x, y, colour="lightgray", edgecolour="gray", alpha=0.6)
                # ax.scatter(points_mat[:, 0]*100, points_mat[:,1]*100, c = colours_ex, s = 800, marker = 's')                
                # **************************************************************************************

                # for each child, assign the SL colours from the corresponding node and parent hue. Update colours
                for child_idx in range(0, len(children_idx)):
                    # for this child get the point to assign
                    child_SL = points_mat[point_assignment[child_idx]]
                    child_hex = convert.hsl_to_hex([phylum_hue, child_SL[0], child_SL[1]])

                    # save to the dataframe (at the overall index across all requested, updated the "Colour" to the child's hex
                    requested_colours_save.at[children_idx[child_idx], "Colour"] =  "#" + str(child_hex)

                # PLOT FOR PAPER *************************************************************************
                    # get the lineage of the child
                #     label = [tax.name(t) for t in children_lineages[child_idx] if tax.rank(t) == "order"]


                #     plt.annotate(label, # this is the text
                #                 (child_SL[0], child_SL[1]), # these are the coordinates to position the label
                #                 textcoords="offset points", # how to position the text
                #                 xytext=(0,10), # distance from text to points (x,y)
                #                 ha='center') # horizontal alignment can be left, right or center
                # plt.show()
                # plt.savefig()
                # ****************************************************************************************
        # if there is only 1 child, give phylum colour exactly
            else:
                requested_colours_save.at[children_idx[0], "Colour"] =  "#" + str(phylum_hex_colour)


    # return the dataframe of requested colours
    return requested_colours_save


def getRequestedColours(database_name: str, database_version: Union[str, list], requested_taxa: Union[str, pd.DataFrame], name_save_colours: str, name_provided_colours = "", mode = "global"):
    # PURPOSE: call colour generation code
    # INPUT:
    #       database_name (string): the name of the database you would like all phyla for (ncbi or gtdb)
    #       database_version (string or list): 
    #           string: the release version, or empty for the latest version .
    #           list: for a custom database.  A list of paths (as strings) to files for database.
    #           (ex. ["nodes.dmp", "names.dmp"] or ["tarDump.tar.gz"]) - See https://github.com/pirovc/multitax for more information.
    #       requested_taxa (string or dataframe): path to a CSV containing column 'Taxon' that is the scientific name of all taxa you need a colour for OR 
    #           a dataframe of the same format.
    #       name_save_colours (string): name to save requested colours (see OUTPUT) as CSV. If do not want to save CSV, use "".
    #       name_provided_colours (string): OPTIONAL: A path to a CSV of parent colours (at least phyla needed) in format of 2 columns:
    #           - 'Taxon' contains scientific name of phyla (and additional parent colours at lower taxonomic levels if requested).
    #           - 'Colour' contains hex code (with leading "#") of the corresponding taxon.
    #       mode (string): OPTIONAL: can be global (always generate the same colour for the same taxon) or "dataset-specific" (generate colours such that all 
    #           taxa in the same phyla are maximally different). Default is "global".
    # OUTPUT: 
    #   requested (pd.DataFrame): a dataframe of 2 columns 
    #           - 'Taxon' same as Taxon column from requested_colours.
    #           - 'Colour' assigned colour.
    # ensure database name is lowercase
    database_name = database_name.lower()


    # if a correct database name was supplied
    if database_name == "ncbi" or database_name == "gtdb":

        # if no parent colours were provided, use default colours and warn user
        if (not name_provided_colours) and database_name == "ncbi":
            provided_colours = pd.read_csv("phylaNCBIMay2026.csv", header='infer')
            print("Using default colours for NCBI as of May 2026")
        
        elif (not name_provided_colours)  and database_name == "gtdb":
            provided_colours = pd.read_csv("phylaGTDB226.csv", header='infer')
            print("Using default colours for GTDB release 226")

        # if parent colours were provided, use the colours requested by user
        else:
            provided_colours = pd.read_csv(name_provided_colours, header='infer')


    # if a database name was requested that is not supported stop the program and return an error
    else:
        print("Error: Requested database" +  database_name + " is not available")
        return 1
    
    if isinstance(requested_taxa, pd.DataFrame):
        requested_colours = requested_taxa
    elif isinstance(requested_taxa, str):
        requested_colours = pd.read_csv(requested_taxa, header='infer')
    else:
        print("ERROR: requested taxa must be either a string that is the path to a CSV or a dataframe")
    

    if mode == "global":
        requested = generateColoursGlobal(provided_colours, requested_colours, database_name, database_version)
    elif mode == "dataset-specific":
         requested = generateColoursSpecific(provided_colours, requested_colours, database_name, database_version)
    else:
        print("Requested mode not available. Please use global or dataset-specific")


    # if a save name was provided, save the dataframe
    if name_save_colours:
        requested.to_csv(name_save_colours, index=False)

    # return the dataframe
    return requested