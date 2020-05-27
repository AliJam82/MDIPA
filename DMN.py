import networkx as nx
import numpy
from operator import itemgetter, attrgetter
from datetime import datetime



def DMN_N():
    graph = nx.Graph()
    m = []
    d = []
    dmn_graph = []
    dmn_final = []
    m_name_array = []
    d_name_array = []
    element_array = []
    threshold = 3
    index = []
    index_array = []

    with open('mirna_drug.txt') as f:
        all_lines = f.readlines()  # get all the lines
        all_lines = [x.strip() for x in all_lines]
        cols = all_lines[0].split(',')
        d_name_array = cols
        y_matrix = numpy.zeros((len(all_lines) - 1, len(cols) - 1))
        y_matrix_t = numpy.transpose(y_matrix)
        all_lines = numpy.delete(all_lines, 0)
        x_count = 0
        y_row_count = len(all_lines)

        for i in range(len(all_lines)):
            rows = all_lines[i].split(',')
            m_name_array.append(rows[0])
            rows = numpy.delete(rows, 0)
            y_col_count = len(rows)
            for j in range(len(rows)):
                y_matrix[i][j] = rows[j]  # (i-1) and (j-1) indexes used for skipped first row and col
    # m_name_array = numpy.delete(m_name_array, 0)
    d_name_array = numpy.delete(d_name_array, 0)
    print(d_name_array)
    print(m_name_array)

    for i in range(len(m_name_array)):
        for j in range(len(d_name_array)):
            if y_matrix[i][j] == 1:
                element_array.append([m_name_array[i], d_name_array[j], (i*len(d_name_array))+(j+1)])
    numpy.savetxt('element_array.txt', element_array, '%s')

    with open('DM_network.txt') as f:
        all_lines = f.readlines()  # get all the lines
        all_lines = [x.strip() for x in all_lines]  # strip away newline chars
        for i in range(len(all_lines)):
            hh = all_lines[i].split(',')
            try:
                m.append(hh[0])
                d.append(hh[1])
                dmn_graph.append([hh[0], hh[1]])
            except:
                continue

    for i in range(len(dmn_graph)):
        edge_i = dmn_graph[i]
        for j in range(i + 1, len(dmn_graph)):
            edge_j = dmn_graph[j]
            if edge_i[1] == edge_j[1]:
                dmn_graph.append([edge_i[0], edge_j[0]])

    for x in range(len(dmn_graph)):
        edge_x = dmn_graph[x]
        graph.add_edge(edge_x[0], edge_x[1])
    for mm in range(len(m)):
        m_name = m[mm]
        for dd in range(len(d)):
            d_name = d[dd]
            try:
                spath = nx.shortest_path_length(graph, source=m_name, target=d_name)
            except:
                spath = 0
            if spath == 0:
                status = 'Absent'
            elif spath < threshold:
                status = 'interacting'
            else:
                status = 'non-interacting'
                for m_index in range (len (m_name_array)):
                    if m_name_array [m_index] == m_name:
                        mirna_index = m_index
                for d_index in range (len (d_name_array)):
                    if d_name_array [d_index] == d_name:
                        drug_index = d_index
                element_array.append([m_name, d_name, (mirna_index*(len(d_name_array)))+(drug_index+1)])
            dmn_final.append([m_name, d_name, spath, status])
        print(mm, "This is printed on: ", str(datetime.now()))
    numpy.savetxt("element_array_ALL.txt", element_array, "%s")
    dmn_final_sorted = sorted(dmn_final, key=itemgetter(2))
    numpy.savetxt('DMN_All.txt', dmn_final_sorted, '%s')
    for i_index in range (len(element_array)):
        element_array_col = element_array[i_index]
        index_array.append(int(element_array_col[2]))
    numpy.savetxt("index_array_ALL.txt", index_array, '%d')
    return index_array


if __name__ == '__main__':
    DMN_N()

