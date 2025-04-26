import os

import numpy as np
import xml.etree.ElementTree as ET

from matplotlib import pyplot as plt

# -------------------------------Settings----------------------------------------------
BASE_PATH = r"files/"
SINGLE_FILE_NAME = "0.xml"
SINGLE_MITO = False
SINGLE_MITO_ID = 15
# https://matplotlib.org/stable/users/explain/colors/colormaps.html
COLOR_SCHEME = "binary"

bin_sizes = {
    1: {
        "1": (58, 1271*2),
        "2": (58, 1271*2),
        "3": (58, 1271*2),
        "4": (58, 1271*2),
        "5": (58, 58),
        "6": (58, 58),
    },
    0: {
        "1": (58, 58),
        "2": (58, 58),
        "3": (58, 58),
        "4": (58, 58),
        "5": (58, 58),
        "6": (58, 58)
    }
}

mito_sizes = {
    "1": {
        "1": [0.580003, 25.4281],
        "2": [0.580003, 25.4281],
        "3": [0.580003, 25.4281],
        "4": [0.580003, 25.4281],
        "5": [0.580003, 0.580003],
        "6": [0.580003, 0.580003]
    },
    "0": {
        "1": [0.580003, 0.580003],
        "2": [0.580003, 0.580003],
        "3": [0.580003, 0.580003],
        "4": [0.580003, 0.580003],
        "5": [0.580003, 0.580003],
        "6": [0.580003, 0.580003]
    }

}


def get_all_mitos():
    """
    Fetch all Mitos with all base parameters.
    :return:
    """
    mito_infos = dict()
    tree = ET.parse(BASE_PATH + SINGLE_FILE_NAME)
    root = tree.getroot()

    for xagent in root.findall('xagent'):
        if xagent.find('name').text == "Mito":
            mito_id = int(xagent.find('id').text)
            mito_infos[mito_id] = dict()
            mito_infos[mito_id]["mito_size"] = int(xagent.find("current_mitosize").text)
            mito_infos[mito_id]["x-length"] = float(xagent.find('mitolength').text)
            mito_infos[mito_id]["y-length"] = float(xagent.find('mitoheight').text)
            mito_infos[mito_id]["z-length"] = float(xagent.find('mitowidth').text)
            for i in range(1, 9):
                mito_infos[mito_id][f"x{i}"] = float(xagent.find(f"x{i}").text)
                mito_infos[mito_id][f"y{i}"] = float(xagent.find(f"y{i}").text)
                mito_infos[mito_id][f"z{i}"] = float(xagent.find(f"z{i}").text)

    return mito_infos


index_to_ignore = {
    "1": 1,
    "2": 0,
    "3": 1,
    "4": 0,
    "5": 2,
    "6": 2,
}


def get_all_points():
    tree = ET.parse(BASE_PATH + SINGLE_FILE_NAME)
    root = tree.getroot()

    points_receptors = dict()

    for xagent in root.findall('xagent'):
        if xagent.find('name').text == "Receptor":
            mito_id = xagent.find('mitoid').text
            if SINGLE_MITO and mito_id != str(SINGLE_MITO_ID):
                continue

            x = float(xagent.find('x').text)
            y = float(xagent.find('y').text)
            z = float(xagent.find('z').text)

            edge = xagent.find('foundedge').text

            # if the edge does not yet exist, add a dictionary underneath in which Mito_ID
            # is used as the key and a list with coordinates is saved underneath it
            if points_receptors.get(edge) is None:
                points_receptors[edge] = dict()
            # add empty list
            if points_receptors[edge].get(mito_id) is None:
                points_receptors[edge][mito_id] = list()

            # removes the coordinate, which is the same for every edge
            # in order to be able to make 2-D histograms later on
            coords = [x, y, z]
            del coords[index_to_ignore[edge]]
            points_receptors[edge][mito_id].append(coords)

    return points_receptors


dict_receptors = get_all_points()
dict_mitos = get_all_mitos()

# assign each edge as key with a list of min point and max point id (aka xi,yi,zi)
min_max_points_edge = {
    "1": ["1", "6"],
    "2": ["2", "8"],
    "3": ["3", "8"],
    "4": ["1", "7"],
    "5": ["1", "4"],
    "6": ["5", "8"]
}


def count_sublists(data):
    """
    This function counts the occurrences of each sublist in a list.

    Args:
        data: A list of sublists.

    Returns:
        A dictionary where the keys are the unique sublists in the list and the values are their corresponding counts.
    """

    sublist_counts = {}
    for sublist in data:
        # Convert the sublist to a hashable tuple for dictionary lookup
        hashable_sublist = tuple(sublist)
        if hashable_sublist in sublist_counts:
            sublist_counts[hashable_sublist] += 1
        else:
            sublist_counts[hashable_sublist] = 1
    return sublist_counts


point_list = list()
dict_histograms = {
    "1": {"0": [], "1": []},
    "2": {"0": [], "1": []},
    "3": {"0": [], "1": []},
    "4": {"0": [], "1": []},
    "5": {"0": [], "1": []},
    "6": {"0": [], "1": []}
}

# Put all points in buckets and write them to dict_histograms
for edge in sorted(dict_receptors.keys()):
    for mito in sorted(list(dict_receptors[edge].keys()), key=lambda x: int(x)):
        mito_size = dict_mitos[int(mito)]["mito_size"]

        point_list = np.array(dict_receptors[edge][mito])

        # delete useless coord-index
        coords = ["x", "y", "z"]
        del coords[index_to_ignore[edge]]
        # fetch lower bound
        lower_bound = list()
        for axis in coords:
            lower_bound.append(dict_mitos[int(mito)][axis + min_max_points_edge[edge][0]])

        # scale points by subtracting the lower bound from all points
        point_list_scaled = np.subtract(point_list, lower_bound)

        # asd
        bucket_size = bin_sizes[mito_size][edge]
        factor_1 = dict_mitos[int(mito)][f"{coords[0]}-length"]
        factor_2 = dict_mitos[int(mito)][f"{coords[1]}-length"]

        bucket_factors = np.divide([factor_1, factor_2], bucket_size)

        point_list_buckets = np.round(np.divide(point_list_scaled, bucket_factors)).astype(int)

        # count for each bucket how often items are in them
        count_sublist = count_sublists(list(point_list_buckets))

        # create new array with bucket size as dimensions filled with 0
        grid = np.zeros(bucket_size)

        for sublist in count_sublist.keys():
            grid[sublist[0] - 1][sublist[1] - 1] = count_sublist[sublist]

        if not SINGLE_MITO:
            dict_histograms[edge][str(mito_size)].append(grid)
        else:
            dict_histograms[edge][str(mito_size)] = grid

# iterate over all edge / size combinations and get the mean
for edge in sorted(dict_histograms.keys()):
    for mito_size in sorted(list(dict_histograms[edge].keys()), key=lambda x: int(x)):
        # to scale bins variably to particles per surface area
        y_size, x_size = mito_sizes[mito_size][edge]
        y_bin_amount, x_bin_amount = bin_sizes[int(mito_size)][edge]
        divisor = (y_size / y_bin_amount) * (x_size / x_bin_amount)
        if not SINGLE_MITO:
            mean_per_bucket = np.mean(dict_histograms[edge][mito_size], axis=0)
            # apply scaling factor
            mean_per_bucket /= divisor
        else:
            mean_per_bucket = np.array(dict_histograms[edge][mito_size])
            mean_per_bucket /= divisor

        dict_histograms[edge][mito_size] = mean_per_bucket


for edge in sorted(dict_histograms.keys(), key=lambda x: int(x)):
    for mito_size in sorted(list(dict_histograms[edge].keys()), key=lambda y: int(y)):
        H = dict_histograms[edge][mito_size]
        shape = H.shape
        try:
            X, Y = np.meshgrid(list(range(shape[1])), list(range(shape[0])))
        except IndexError:
            continue

        figure_size = (12.8, 4.8) if edge not in ["5", "6"] and mito_size == '1' else (6.4, 4.8)
        fig = plt.figure(figsize=figure_size, layout='tight')
        plt.pcolormesh(X, Y, H, cmap=COLOR_SCHEME, vmax=500)
        plt.colorbar()
        plt.title(f"{'fragmented' if mito_size == '0' else 'non-fragmented'} mitochondria | Edge {edge}")
        plt.savefig(
            f"./output/{os.path.basename(os.path.normpath(BASE_PATH))}_{SINGLE_FILE_NAME.replace('.xml', '')}_edge{edge}_size{mito_size}.png")
        plt.clf()
