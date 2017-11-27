import math
import random

import time


def dist(node1, node2):
    """
    L2 euclidean distance between coord1 and coord2
    :param node1: tuple (x1, y1)
    :param node2: tuple (x2, y2)
    :return: Euclidean distance between coord1 and coord2
    """

    x1, y1 = node1
    x2, y2 = node2
    dx = x2 - x1
    dy = y2 - y1
    return math.sqrt(dx ** 2 + dy ** 2)


def create_distances_dict(coord):
    """
    Function for creating distance dict
    :param coord: coordinates, tuple [(x1, y1), (x2, y2), ..]
    :return: dict of distances between all coordinates
    """
    distances = dict()  # dictionary to hold n times n distances
    for i in range(len(coord) - 1):
        for j in range(i + 1, len(coord)):
            distances[i, j] = dist(coord[i], coord[j])
            distances[j, i] = distances[i, j]
    return distances


def closest_dist(distances, n):
    """
    Create sorted list of closest distances
    :param distances: dictionary of distances between nodes
    :param n: number of nodes
    :return: Sorted list of closest distances
    """

    closest = []
    for i in range(n):
        dists = [(distances[i, j], j) for j in range(n) if j != i]
        dists.sort()
        closest.append(dists)

    return closest


def length(tour, distances):
    """
    Calculate total length of the tour
    :param tour: list of coords in visiting order
    :param distances: distances dict
    :return: total distance
    """

    total_dist = distances[tour[-1], tour[0]]
    for i in range(1, len(tour)):
        total_dist += distances[tour[i], tour[i - 1]]
    return total_dist


def rand_order(nodes):
    """
    Return list in random order of n elements
    :param nodes: nodes list
    :return: shuffeled list of n elements
    """
    rand = list(range(nodes))
    random.shuffle(rand)  # place it in a random order
    return rand


def exchange(tour, pos_city, i, j):
    """
    Calculate new tour where vertices (i, i+1) and (j, j+1) are changed to
    (i, j), (i+1, j+1)
    :param tour: tour with list of nodes
    :param pos_city: position of each city in tour
    :param i: position of the first vertice
    :param j: position of the second vertice
    """

    if i > j:
        i, j = j, i
    assert 0 <= i < j - 1 and j < len(tour)
    path = tour[i + 1:j + 1]
    path.reverse()
    tour[i + 1:j + 1] = path
    for k in range(i + 1, j + 1):
        pos_city[tour[k]] = k


# Partly from https://stackoverflow.com/questions/30552656/python-traveling-salesman-greedy-algorithm
def improve(tour, dist, distances, closest):
    """
    Try to improve our current tour by making local improvements to improve the global tour
    :param tour: tour with list of nodes
    :param dist: length of the initial tour
    :param distances: distances dictionary
    :param closest: Sorted list of distances to each node
    :return: length of the tour
    """

    n = len(tour)
    pos_city = [0 for i in tour]
    for k in range(n):
        pos_city[tour[k]] = k  # position of each city in tour
    for i in range(n):
        a, b = tour[i], tour[(i + 1) % n]
        dist_ab = distances[a, b]
        improved = False
        for dist_ac, c in closest[a]:
            if dist_ac >= dist_ab:  # There is no improvement even in local scale
                break
            j = pos_city[c]
            d = tour[(j + 1) % n]
            dist_cd = distances[c, d]
            dist_bd = distances[b, d]
            delta = (dist_ac + dist_bd) - (dist_ab + dist_cd)
            if delta < 0:  # Better solution found
                exchange(tour, pos_city, i, j)
                dist += delta
                improved = True
                break
        if improved:
            continue
        for dist_bd, d in closest[b]:
            if dist_bd >= dist_ab:
                break
            j = pos_city[d] - 1
            if j == -1:
                j = n - 1
            c = tour[j]
            dist_cd = distances[c, d]
            dist_ac = distances[a, c]
            delta = (dist_ac + dist_bd) - (dist_ab + dist_cd)
            if delta < 0:  # exchange decreases length
                exchange(tour, pos_city, i, j)
                dist += delta
                break
    return dist


def localsearch(tour, dist, distances, closes_dist=None):
    """
    Perform a local search until we reach a local minimum
    :param tour: Current tour
    :param dist: Total dist of tour
    :param distances: distances dictionary
    :param closes_dist: dict of closest distances
    :return: Best tour length 
    """
    n = len(tour)
    if closes_dist is None:
        closes_dist = closest_dist(distances, n)  # create a sorted list of distances to each node
    while True:
        new_dist = improve(tour, dist, distances, closes_dist)
        if new_dist < dist:
            dist = new_dist
        else:
            break
    return dist


def multistart_localsearch(iter, n, distances):
    """
    Do multiple iteration of local search starting from random solutions
    :param iter: number of iterations
    :param n: number of nodes
    :param distances: distances dictionary
    :return: best tour and it's total distance
    """
    closest = closest_dist(distances, n)  # create a sorted list of distances to each node
    best_tour = None
    best_dist = None

    for i in range(0, iter):
        tour = rand_order(n)
        dist = length(tour, distances)
        dist = localsearch(tour, dist, distances, closest)
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_tour = list(tour)

    return best_tour, best_dist


def nearest_node(curr_node, unvisited, distances):
    """
    Greedy implementation to find nearest neighbour to current node
    :param curr_node: Current node int
    :param unvisited: list of unvisited nodes
    :param distances: distances between the nodes
    :return: nearest node
    """
    nearest = unvisited[0]
    min_dist = distances[curr_node, nearest]
    for site in unvisited:
        if distances[curr_node, site] < min_dist:
            nearest = site
            min_dist = distances[curr_node, nearest]
    return nearest


def nearest_neighbor(nodes, start, distances):
    """
    Nearest neighbour implementation to find greedy smallest path to all nodes
    :param nodes: list of node id-s
    :param start: starting point of the tour, node id: int
    :param distances: distances dictionary
    :return: Greedy solution of TSP tour
    """

    unvisited = nodes
    unvisited.remove(start)
    curr_node = start
    tour = [start]
    while unvisited:
        next_site = nearest_node(curr_node, unvisited, distances)
        tour.append(next_site)
        unvisited.remove(next_site)
        curr_node = next_site
    return tour

def tour_to_coord(tour, coord):
    """
    Convert tour in indexes back to coordinates of points
    :param tour: list of nodes in tour
    :param coord: list of coordinates
    :return: tour of coordinates
    """
    output = []
    for elem in tour:
        output.append(coord[elem])

    return output


def read_coords(filename, header=False):
    """
    Read coordinates from file
    :param filename: file path
    :param header: if file contains header
    :return: list of coordinates tuple
    """
    coords = []
    with open(filename) as f:
        if header:
            next(f)
        for line in f:
            x, y = line.strip().split(" ")
            coords.append(tuple((float(x), float(y))))
    return coords


def write_swog_lines(tour, filename):
    """
    Create swog file of the tour lines.
    :param tour: tour of indexes list
    """
    f = open(filename, "w")
    for i in range(len(tour) - 1):
        curr_node = tour[i]
        next_node = tour[i + 1]
        f.write("line (p" + str(curr_node + 1) + ") (p" + str(next_node + 1) + ")\n")
    f.write("line (p" + str(tour[-1] + 1) + ") (p" + str(tour[0] + 1) + ")\n")
    f.close()


def multistart_tsp(coord_filename, output_file_name, niter):
    """
    Calculate multistart tsp path
    :param coord_filename: filepath for coords file
    :param output_file_name: filename for outputs of swog lines
    :param niter: number of iterations
    :return: tour
    """
    coord_list = read_coords(coord_filename, True)
    dist_dict = create_distances_dict(coord_list)
    start = time.time()
    tour, z = multistart_localsearch(niter, len(coord_list), dist_dict)
    end = time.time()
    print("Total dist, " + str(len(coord_list)) + " cities:", z)
    print("Total time:", end - start)
    write_swog_lines(tour, output_file_name)
    return tour


def nearest_neighbour_tsp(coord_filename, output_filename):
    """
    Calculate nearest neighbour tsp path
    :param coord_filename: filepath for coords file
    :param output_filename: filename for outputs of swog lines
    :return: tour
    """
    coord_list = read_coords(coord_filename, True)
    dist_dict = create_distances_dict(coord_list)
    start = time.time()
    tour = nearest_neighbor(list(range(len(coord_list))), 0, dist_dict)
    z = length(tour, dist_dict)
    end = time.time()
    print("NN Total dist, " + str(len(coord_list)) + " cities:", z)
    print("NN Total time:", end - start)
    write_swog_lines(tour, output_filename)
    return tour


if __name__ == "__main__":

    # 1000 Cities case
    multistart_tsp("data/TSP_1000.txt", "swog_path_1000.txt", 100)

    # 100 Cities case
    multistart_tsp("data/TSP_100.txt", "swog_path_100.txt", 100)

    # 20 Cities case
    multistart_tsp("data/TSP_20.txt", "swog_path_20.txt", 100)

    # ____________________________________________________________

    # 1000 Cities case
    nearest_neighbour_tsp("data/TSP_1000.txt", "swog_path_1000_NN.txt")

    # 100 Cities case
    nearest_neighbour_tsp("data/TSP_100.txt", "swog_path_100_NN.txt")

    # 20 Cities case
    nearest_neighbour_tsp("data/TSP_20.txt", "swog_path_20_NN.txt")