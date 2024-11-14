"""
Topological Sort Implementation

  File: topo_sort.py
  Description:

  Student Name: Marissa Shuchart
  Student UT EID: ms87339

  Partner Name: Trinity Thompson
  Partner UT EID: tyt242

  Course Name: CS 313E
  Unique Number: 50165
  Date Created: 11/13/2024
  Date Last Modified:

"""

import sys

class Stack():
    """
    Stack class; you can use this for your search algorithms  
    """
    def __init__(self):
        self.stack = []

    def push(self, item):
        """
        Add an item to the top of the stack
        """
        self.stack.append(item)

    def pop(self):
        """
        Remove an item from the top of the stack
        """
        return self.stack.pop()

    def peek(self):
        """
        Check the item on the top of the stack
        """
        return self.stack[-1]

    def is_empty(self):
        """
        Check if the stack if empty
        """
        return len(self.stack) == 0

    def size(self):
        """
        Return the number of elements in the stack
        """
        return len(self.stack)


class Queue():
    """
    Queue class; you can use this for your search algorithms
    """
    def __init__(self):
        self.queue = []

    def enqueue(self, item):
        """
        Add an item to the end of the queue
        """
        self.queue.append(item)

    def dequeue(self):
        """
        Remove an item from the beginning of the queue
        """
        return self.queue.pop(0)

    def peek(self):
        """
        Checks the item at the top of the Queue
        """
        return self.queue[0]

    def is_empty(self):
        """
        Check if the queue is empty
        """
        return len(self.queue) == 0

    def size(self):
        """
        Return the size of the queue
        """
        return len(self.queue)

class Vertex():
    """Vertex Class"""
    def __init__(self, label):
        self.label = label
        self.visited = False

    
    def was_visited(self):
        """Determine if a vertex was visited"""
        return self.visited

    def get_label(self):
        """Determine the label of the vertex"""
        return self.label

    def __str__(self):
        """String representation of the vertex"""
        return str(self.label)


class Graph(object):
    """A Class to present Graph."""

    def __init__(self):
        self.vertices = []  # a list of vertex objects
        self.adj_mat = []  # adjacency matrix of edges

    def has_vertex(self, label):
        """Check if a vertex is already in the graph"""
        n_vert = len(self.vertices)
        for i in range(n_vert):
            if label == (self.vertices[i]).get_label():
                return True
        return False

    def get_index(self, label):
        """Given a label get the index of a vertex"""
        n_vert = len(self.vertices)
        for i in range(n_vert):
            if label == (self.vertices[i]).get_label():
                return i
        return -1

    def add_vertex(self, label):
        """Add a Vertex with a given label to the graph"""
        if self.has_vertex(label):
            return

        # add vertex to the list of vertices
        self.vertices.append(Vertex(label))

        # add a new column in the adjacency matrix
        n_vert = len(self.vertices)
        for i in range(n_vert - 1):
            (self.adj_mat[i]).append(0)

        # add a new row for the new vertex
        new_row = []
        for i in range(n_vert):
            new_row.append(0)
        self.adj_mat.append(new_row)

    def add_directed_edge(self, start, finish, weight=1):
        """Add weighted directed edge to graph"""
        self.adj_mat[start][finish] = weight

    def add_undirected_edge(self, start, finish, weight=1):
        """Add weighted undirected edge to graph"""
        self.adj_mat[start][finish] = weight
        self.adj_mat[finish][start] = weight

    def get_adj_unvisited_vertex(self, v):
        """Return an unvisited vertex adjacent to vertex v (index)"""
        for i in range(len(self.vertices)):
            if (self.adj_mat[v][i] > 0) and (
                    not (self.vertices[i]).was_visited()):
                return i
        return -1

    def get_adj_vertexes(self, v):
        """Return an adjacent vertex  to vertex v (index)"""
        verts = []
        n_vert = len(self.vertices)
        for i in range(n_vert):
            if self.adj_mat[v][i] > 0:
                verts.append(i)
        return verts

    def get_adj_back_forth_vertex(self, v):
        """get back adj vertex."""
        verts = []
        n_vert = len(self.vertices)
        for i in range(n_vert):
            if (self.adj_mat[v][i] > 0) or self.adj_mat[i][v] > 0:
                verts.append(i)
        return verts


    def dfs(self, v):
        """Do the depth first search in a graph from vertex v (index)"""
        # create the Stack
        the_stack = Stack()

        # cycle check
        cyclic = False

        # mark the vertex v as visited and push it on the stack
        (self.vertices[v]).visited = True
        # print(self.Vertices[v])
        the_stack.push(v)

        # visit the other vertices according to depth
        while (not the_stack.is_empty()) and not cyclic:
            # get an adjacent unvisited vertex
            u = self.get_adj_unvisited_vertex(the_stack.peek())
            # print(u)
            # print(theStack.__str__())
            adjacents = self.get_adj_vertexes(u)
            # print(adjacents)
            if v in adjacents:
                # print(v)
                final_adjacents = self.get_adj_back_forth_vertex(v)
                # print(final_adjacents)
                # print(u)
                if u in final_adjacents:
                    cyclic = True
            if u == -1:
                u = the_stack.pop()
            else:
                (self.vertices[u]).visited = True
                the_stack.push(u)


        # the stack is empty, let us reset the flags
        for i in range(len(self.vertices)):
            (self.vertices[i]).visited = False

        return cyclic
        # determine if a directed graph has a cycle
        # this function should return a boolean and not print the result

    def has_cycle(self):
        """If the graph has cycle"""
        # Call existing DFS function for each vertex to find any possible cycles
        for v_index, _ in enumerate(self.vertices):
            if self.dfs(v_index):
                return True
        return False


    def bfs(self, v):
        """Do the breadth first search in a graph"""
        the_queue = Queue()

        (self.vertices[v]).visited = True
        print(self.vertices[v])
        the_queue.enqueue(v)

        # visit the other vertices according to depth
        while not the_queue.is_empty():
            # get an adjacent unvisited vertex
            u = self.get_adj_unvisited_vertex(the_queue.current())
            if u == -1:
                u = the_queue.dequeue()
            else:
                (self.vertices[u]).visited = True
                print(self.vertices[u])
                the_queue.enqueue(u)

        # the stack is empty, let us reset the flags
        for i in range( len(self.vertices)):
            (self.vertices[i]).visited = False


    def delete_edge(self, from_vertex_label, to_vertex_label):
        """
        Delete an edge from the adjacency matrix
        Delete a single edge if the graph is directed
        Delete two edges if the graph is undirected
        """
        start = self.get_index(from_vertex_label)
        finish = self.get_index(to_vertex_label)
        if self.adj_mat[start][finish] == self.adj_mat[finish][start]:
            self.adj_mat[start][finish] = 0
            self.adj_mat[finish][start] = 0
        else:
            self.adj_mat[start][finish] = 0

    def delete_vertex(self, vertex_label):
        """
        Delete a vertex from the vertex list and all edges from and
        to it in the adjacency matrix    
        """

        idx = self.get_index(vertex_label)
        for row in self.adj_mat:
            row.pop(idx)
        self.adj_mat.pop(idx)
        self.vertices.pop(idx)

    def toposort(self):
        """
        Return a list of vertices after a topological sort
        this function should not print the list
        """

        # Calculate the in-degrees of all vertices
        in_degrees = {i: 0 for i in range(len(self.vertices))}
        for u in range(len(self.vertices)):
            for v in self.get_adj_vertexes(u):
                in_degrees[v] += 1

        # Initialize queue 
        queue = deque([v for v in range(len(self.vertices)) if in_degrees[v] == 0])
        topo_list = []

        # Process each node 
        while queue:
            v = queue.popleft()
            topo_list.append(self.vertices[v].get_label())
            for neighbor in self.get_adj_vertexes(v):
                in_degrees[neighbor] -= 1
                if in_degrees[neighbor] == 0:
                    queue.append(neighbor)
        
        # Check if sort is possible
        if len(topo_list) != len(self.vertices):
            return None 
      
        return topo_list

    def get_index2(self, label, vertices_copy):
        """Given a label get the index of a vertex"""
        for i in range(len(vertices_copy)):
            if label == (vertices_copy[i]).get_label():
                return i
        return -1

    def delete_vertex2(self, vertex_label, adj_mat_copy, vertices_copy):
        """delete vertex """
        idx = self.get_index2(vertex_label, vertices_copy)
        # print(self.adjMat[0])
        # print(self.Vertices)
        for row in adj_mat_copy:
            row.pop(idx)
        adj_mat_copy.pop(idx)
        vertices_copy.pop(idx)


def main():
    """A main function to run a test."""
    # create a Graph object
    the_graph = Graph()

    # read the number of vertices
    line = sys.stdin.readline()
    line = line.strip()
    num_vertices = int(line)

    # read the vertices and insert them into the graph
    for _ in range(num_vertices):
        line = sys.stdin.readline()
        vertex = line.strip()
        the_graph.add_vertex(vertex)

    # read the number of edges
    line = sys.stdin.readline()
    line = line.strip()
    num_edges = int(line)

    # read the edges and insert them into the graph
    for _ in range(num_edges):
        line = sys.stdin.readline()
        line = line.strip()
        edge = line.split()
        # print(edge)
        start = the_graph.get_index(edge[0])
        finish = the_graph.get_index(edge[1])
        # print(start, finish)

        the_graph.add_directed_edge(start, finish, 1)

    # print(num_edges)
    # test if a directed graph has a cycle
    if the_graph.has_cycle():
        print("The Graph has a cycle.")
    else:
        print("The Graph does not have a cycle.")

    # test topological sort
    if not the_graph.has_cycle():
        vertex_list = the_graph.toposort()
        print("\nList of vertices after toposort")
        print(vertex_list)

main()
