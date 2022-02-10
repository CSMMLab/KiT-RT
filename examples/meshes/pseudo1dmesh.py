### manually create pseudo 1d su2 format mesh
import numpy as np
from optparse import OptionParser


def main():
    print("---------- Start Creating the mesh ------------")
    print("| Parsing options")
    # --- parse options ---
    parser = OptionParser()
    parser.add_option("-o", "--output_name", dest="output_name", default="1dmesh")
    parser.add_option("-c", "--char_length", dest="char_length", default=0.01)
    parser.add_option("-s", "--start_pt", dest="start_pt", default=0)
    parser.add_option("-l", "--x_length", dest="x_length", default=1)
    (options, args) = parser.parse_args()

    options.output_name = str(options.output_name)
    options.char_length = float(options.char_length)
    options.x_length = float(options.x_length)
    char_length = options.char_length
    lengthX = options.x_length
    start_pt = float(options.start_pt)
    n_elem = int(lengthX / char_length)
    n_pts = 2 * (int(lengthX / char_length) + 1)
    list_nodes = create_nodes(start_pt, lengthX, char_length).tolist()
    list_elems = create_elements(lengthX, char_length).tolist()
    list_markerdirichlet, list_markerneurmann = create_markers(lengthX, char_length)
    with open(options.output_name + ".su2", 'w') as meshfile:
        meshfile.write("NDIME= 2\n")
        meshfile.write("NELEM= " + str(n_elem) + "\n")
        for elem in list_elems:
            meshfile.write(
                str(elem[0]) + " " + str(elem[1]) + " " + str(elem[2]) + " " + str(elem[3]) + " " + str(elem[4]) + "\n")
        meshfile.write("NPOIN= " + str(n_pts) + "\n")
        for point in list_nodes:
            meshfile.write(
                str(point[0]) + " " + str(point[1]) + " " + str(int(point[2])) + "\n")
        meshfile.write("NMARK= 2" + "\n")
        meshfile.write("MARKER_TAG= dirichlet" + "\n")
        meshfile.write("MARKER_ELEMS= " + str(len(list_markerdirichlet)) + "\n")
        for marker in list_markerdirichlet:
            meshfile.write(
                str(marker[0]) + " " + str(marker[1]) + " " + str(marker[2]) + "\n")
        meshfile.write("MARKER_TAG= wall_neumann" + "\n")
        meshfile.write("MARKER_ELEMS= " + str(len(list_markerneurmann)) + "\n")
        for marker in list_markerneurmann:
            meshfile.write(
                str(marker[0]) + " " + str(marker[1]) + " " + str(marker[2]) + "\n")
    print("| Mesh created with " + str(n_elem) + " elements and " + str(n_pts) + " nodes")
    print("---------- Successfully created the mesh ------------")


def create_nodes(start, length, char_length) -> np.ndarray:
    n_pts = int(length / char_length) + 1
    x_coord = np.linspace(start, start + length, n_pts)
    y_coords = np.concatenate([np.zeros(shape=(n_pts)), np.ones(shape=(n_pts))], axis=0)
    x_coords = np.concatenate([x_coord, x_coord], axis=0)
    indices = np.asarray(list(range(n_pts * 2)))
    list_nodes = np.stack([x_coords, y_coords, indices], axis=1)

    return list_nodes


def create_elements(length, char_length) -> np.ndarray:
    n_elems = int(length / char_length)
    n_pts = int(length / char_length) + 1

    list_elems = []
    for i in range(n_elems):
        list_elems.append([9, i, i + 1, i + n_pts, i + n_pts + 1, i])
    return np.asarray(list_elems)


def create_markers(length, char_length) -> tuple:
    n_elems = int(length / char_length)
    n_pts = int(length / char_length) + 1
    list_neumann = []
    for i in range(n_elems):
        list_neumann.append([3, i, i + 1])
    for i in range(n_elems):
        list_neumann.append([3, i + n_pts, i + n_pts + 1])

    list_dirichlet = []
    list_dirichlet.append([3, 0, n_pts])
    list_dirichlet.append([3, n_pts - 1, 2 * n_pts - 1])

    return np.asarray(list_dirichlet), np.asarray(list_neumann)


if __name__ == '__main__':
    main()
