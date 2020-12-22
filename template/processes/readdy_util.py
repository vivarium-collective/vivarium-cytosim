import numpy as np
import scipy.linalg as linalg
import random
import readdy

class ReaddyUtil:
    '''
    Utilities for ReaDDy models
    '''
    def __init__(self):
        self.bond_pairs = []
        self.angle_triples = []
        self.dihedral_quads = []
        self.repulse_pairs = []

    @staticmethod
    def normalize(v):
        '''
        normalize a vector
        '''
        return v / np.linalg.norm(v)

    @staticmethod
    def get_angle_between_vectors(v1, v2):
        '''
        get the angle between two vectors in radians
        '''
        return np.arccos(np.clip(np.dot(ReaddyUtil.normalize(v1), ReaddyUtil.normalize(v2)), -1., 1.))

    @staticmethod
    def rotate(v, axis, angle):
        '''
        rotate a vector around axis by angle (radians)
        '''
        rotation = linalg.expm(np.cross(np.eye(3), ReaddyUtil.normalize(axis)*angle))
        return np.dot(rotation,np.copy(v))

    @staticmethod
    def vertex_to_string(topology, vertex):
        '''
        get string with type and id for vertex
        '''
        return "{} ({})".format(topology.particle_type_of_vertex(vertex),
                        topology.particle_id_of_vertex(vertex))

    @staticmethod
    def get_non_periodic_boundary_position(pos1, pos2, box_size):
        '''
        if the distance between two positions is greater than box_size,
        move the second position across the box (for positioning calculations)
        '''
        for dim in range(3):
            if abs(pos2[dim] - pos1[dim]) > box_size / 2.:
                pos2[dim] -= pos2[dim]/abs(pos2[dim]) * box_size

        return pos2

    @staticmethod
    def calculate_diffusionCoefficient(r0,eta,T):
        '''
        calculates the theoretical diffusion constant of a spherical particle
            with radius r0[nm]
            in a media with viscosity eta [cP]
            at temperature T [Kelvin]

            returns nm^2/s
        '''
        return ((1.38065*10**(-23) * T)/
            (6*np.pi*eta*10**(-3)*r0*10**(-9))*10**18/10**9)

    @staticmethod
    def calculate_nParticles(C, dim):
        '''
        calculates the number of particles for a species
            at concentration C [uM]
            in cube container with dimensions dim [nm]

            returns unitless number
        '''
        return int(round(C * 1e-30 * 6.022e23 * np.power(dim, 3.)))

    @staticmethod
    def get_vertex_of_type(topology, vertex_type, exact_match):
        '''
        get the first vertex with a given type
        '''
        for vertex in topology.graph.get_vertices():
            pt = topology.particle_type_of_vertex(vertex)
            if ((not exact_match and vertex_type in pt) or
                (exact_match and pt == vertex_type)):
                return vertex
        return None

    @staticmethod
    def get_first_vertex_of_types(topology, vertex_types):
        '''
        get the first vertex with any of the given types
        '''
        for vertex in topology.graph.get_vertices():
            if topology.particle_type_of_vertex(vertex) in vertex_types:
                return vertex
        return None

    @staticmethod
    def get_vertex_with_id(topology, vertex_id):
        '''
        get the first vertex with a given id
        '''
        for vertex in topology.graph.get_vertices():
            if topology.particle_id_of_vertex(vertex) == vertex_id:
                return vertex
        return None

    @staticmethod
    def get_first_neighbor(topology, vertex, exclude_vertices):
        '''
        get the first neighboring vertex
        '''
        exclude_ids = []
        for v in exclude_vertices:
            exclude_ids.append(topology.particle_id_of_vertex(v))

        for neighbor in vertex:
            if topology.particle_id_of_vertex(neighbor) in exclude_ids:
                continue
            return neighbor.get()
        return None

    @staticmethod
    def get_neighbor_of_type(topology, vertex, vertex_type, exact_match):
        '''
        get the first neighboring vertex of type vertex_type
        '''
        for neighbor in vertex:
            v_neighbor = neighbor.get()
            pt = topology.particle_type_of_vertex(v_neighbor)
            if ((not exact_match and vertex_type in pt) or
                (exact_match and pt == vertex_type)):
                return v_neighbor
        return None

    @staticmethod
    def get_neighbor_of_types(topology, vertex, vertex_types, exclude_vertices):
        '''
        get the first neighboring vertex with any of the given types,
        excluding particles with the given ids
        '''
        exclude_ids = []
        for v in exclude_vertices:
            exclude_ids.append(topology.particle_id_of_vertex(v))

        for neighbor in vertex:
            if topology.particle_id_of_vertex(neighbor) in exclude_ids:
                continue
            v_neighbor = neighbor.get()
            pt = topology.particle_type_of_vertex(v_neighbor)
            if pt in vertex_types:
                return v_neighbor
        return None

    @staticmethod
    def get_vertices_of_type(topology, vertex_type, exact_match):
        '''
        get all vertices with a given type
        '''
        v = []
        for vertex in topology.graph.get_vertices():
            pt = topology.particle_type_of_vertex(vertex)
            if ((not exact_match and vertex_type in pt) or
                (exact_match and pt == vertex_type)):
                v.append(vertex)
        return v

    @staticmethod
    def get_neighbors_of_type(topology, vertex, vertex_type, exact_match):
        '''
        get all neighboring vertices with a given type
        '''
        v = []
        for neighbor in vertex:
            v_neighbor = neighbor.get()
            pt = topology.particle_type_of_vertex(v_neighbor)
            if ((not exact_match and vertex_type in pt) or
                (exact_match and pt == vertex_type)):
                v.append(v_neighbor)
        return v

    @staticmethod
    def get_random_vertex_of_type(topology, vertex_type, exact_match):
        '''
        get a random vertex with a given type
        '''
        vertices = ReaddyUtil.get_vertices_of_type(topology, vertex_type, exact_match)

        if len(vertices) == 0:
            return None
        return random.choice(vertices)

    @staticmethod
    def get_random_vertex_of_types(topology, vertex_types):
        '''
        get a random vertex with any of the given types
        '''
        v = []
        for vertex_type in vertex_types:
            v += ReaddyUtil.get_vertices_of_type(topology, vertex_type, True)
        if len(v) == 0:
            return None
        return random.choice(v)

    @staticmethod
    def vertex_satisfies_type(vertex_type, types_include, types_exclude):
        '''
        check if vertex satisfies the type requirements
        '''
        for t in types_include:
            if t not in vertex_type:
                return False
        for t in types_exclude:
            if t in vertex_type:
                return False
        return True

    @staticmethod
    def set_flags(topology, recipe, vertex, add_flags,
                remove_flags, reverse_sort=False):
        '''
        set flags on a vertex
        '''
        pt = topology.particle_type_of_vertex(vertex)

        if "#" not in pt:

            for f in range(len(add_flags)):
                pt = "{}{}{}".format(pt, "_" if f > 0 else "#", add_flags[f])

            recipe.change_particle_type(vertex, pt)

        else:

            flag_string = pt[pt.index("#")+1:]
            flags = flag_string.split("_")
            polymer_indices = ""

            if "tubulin" in pt and len(flags) > 1:
                polymer_indices = "_{}_{}".format(flags[-2], flags[-1])
                flags = flags[:-2]

            for flag in remove_flags:
                if flag in flags:
                    flags.remove(flag)

            for flag in add_flags:
                if flag not in flags:
                    flags.append(flag)

            if len(flags) < 1:
                recipe.change_particle_type(vertex, pt[:pt.index("#")])
                return

            flags.sort(reverse=reverse_sort)

            flag_string = ""
            for f in range(len(flags)):
                flag_string = "{}{}{}".format(
                    flag_string, "_" if f > 0 else "", flags[f])

            recipe.change_particle_type(vertex,
                "{}#{}{}".format(pt[:pt.index("#")], flag_string, polymer_indices))

    @staticmethod
    def calculate_polymer_number(number, offset):
        '''
        calculates the polymer number
            from number
            by offset in [-2, 2]

            returns number in [1,3]
        '''
        n = number + offset
        if n > 3:
            n -= 3
        if n < 1:
            n += 3
        return n

    @staticmethod
    def get_vertex_position(topology, vertex):
        '''
        get the position of a vertex
        '''
        pos = topology.position_of_vertex(vertex)
        return np.array([pos[0], pos[1], pos[2]])

    @staticmethod
    def vertices_are_equal(topology, vertex1, vertex2):
        '''
        check if references are the same vertex
        '''
        return (topology.particle_id_of_vertex(vertex1) ==
            topology.particle_id_of_vertex(vertex2))

    @staticmethod
    def vertices_are_connected(topology, vertex1, vertex2):
        '''
        check if the vertices share an edge
        '''
        for neighbor in vertex1:
            if ReaddyUtil.vertices_are_equal(topology, vertex2, neighbor.get()):
                return True
        return False

    @staticmethod
    def reaction_function_reset_state(topology):
        '''
        reaction function for removing flags from a topology
        '''
        recipe = readdy.StructuralReactionRecipe(topology)

        tt = topology.type
        recipe.change_topology_type(tt[:tt.index("#")])

        return recipe

    @staticmethod
    def rate_function_infinity(topology):
        '''
        rate function for a reaction that should trigger immediately
        whenever reactants are available
        '''
        return 1e30

    def add_bond(self, types1, types2, force_const, bond_length, system):
        '''
        adds a bond to the system (if it hasn't been added already)
            from each type in types1
            to each type in types2
            with force constant force_const
            and length bond_length [nm]
        '''
        for t1 in types1:
            for t2 in types2:
                if (t1, t2) not in self.bond_pairs and (t2, t1) not in self.bond_pairs:

                    system.topologies.configure_harmonic_bond(
                        t1, t2, force_const, bond_length)
                    self.bond_pairs.append((t1, t2))
                    self.bond_pairs.append((t2, t1))

    def add_angle(self, types1, types2, types3, force_const, angle, system):
        '''
        adds an angle to the system (if it hasn't been added already)
            from each type in types1
            through each type in types2
            to each type in types3
            with force constant force_const
            and angle [radians]
        '''
        for t1 in types1:
            for t2 in types2:
                for t3 in types3:
                    if ((t1, t2, t3) not in self.angle_triples and
                        (t3, t2, t1) not in self.angle_triples):

                        system.topologies.configure_harmonic_angle(
                            t1, t2, t3, force_const, angle)
                        self.angle_triples.append((t1, t2, t3))
                        self.angle_triples.append((t3, t2, t1))

    def add_dihedral(self, types1, types2, types3, types4, force_const, angle, system):
        '''
        adds a cosine dihedral to the system (if it hasn't been added already)
            from each type in types1
            through each type in types2
            through each type in types3
            to each type in types4
            with force constant force_const
            and angle [radians]
        '''
        for t1 in types1:
            for t2 in types2:
                for t3 in types3:
                    for t4 in types4:
                        if ((t1, t2, t3, t4) not in self.dihedral_quads and
                            (t4, t3, t2, t1) not in self.dihedral_quads):

                            system.topologies.configure_cosine_dihedral(
                                t1, t2, t3, t4, force_const, 1., angle)
                            self.dihedral_quads.append((t1, t2, t3, t4))
                            system.topologies.configure_cosine_dihedral(
                                t4, t3, t2, t1, force_const, 1., angle)
                            self.dihedral_quads.append((t4, t3, t2, t1))

    def add_repulsion(self, types1, types2, force_const, distance, system):
        '''
        adds a pairwise repulsion to the system (if it hasn't been added already)
            between each type in types1
            and each type in types2
            with force constant force_const
            with equilibrium distance [nm]
        '''
        for t1 in types1:
            for t2 in types2:
                if (t1, t2) not in self.repulse_pairs and (t2, t1) not in self.repulse_pairs:

                    system.potentials.add_harmonic_repulsion(
                        t1, t2, force_const, distance)
                    self.repulse_pairs.append((t1, t2))
                    self.repulse_pairs.append((t2, t1))

    @staticmethod
    def get_random_unit_vector():
        '''
        get a random unit vector
        '''
        return normalize(np.array(
            [random.random(), random.random(), random.random()]))

    @staticmethod
    def get_random_boundary_position(box_size):
        '''
        get a random position on one of the boundary faces
        '''
        pos = box_size * np.random.uniform(size=(3)) - box_size / 2

        face = random.randint(0,5)
        if face == 0:
            pos[0] = -box_size / 2
        elif face == 1:
            pos[0] = box_size / 2
        elif face == 2:
            pos[1] = -box_size / 2
        elif face == 3:
            pos[1] = box_size / 2
        elif face == 4:
            pos[2] = -box_size / 2
        else:
            pos[2] = box_size / 2

        return pos

    @staticmethod
    def try_remove_edge(topology, recipe, vertex1, vertex2):
        '''
        try to remove an edge
        '''
        if not ReaddyUtil.vertices_are_connected(topology, vertex1, vertex2):
            return False, "tried to remove non-existent edge! {} -- {}".format(
                ReaddyUtil.vertex_to_string(topology, vertex1),
                ReaddyUtil.vertex_to_string(topology, vertex2))

        recipe.remove_edge(vertex1, vertex2)
        return True, ""
