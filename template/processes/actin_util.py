import numpy as np
import readdy
import random
from readdy_util import ReaddyUtil
from actin_constants import *

class ActinUtil:
    '''
    Utilities for ReaDDy actin models
    '''
    def __init__():
        self.readdy_util = ReaddyUtil()

    def create_actin_simulation(self, parameters):
        '''
        Create the ReaDDy simulation for actin
        '''
        system = self.create_actin_system(parameters)
        simulation = system.simulation("CPU")
        simulation.kernel_configuration.n_threads = parameters['n_cpu']
        return simulation, system

    @staticmethod
    def add_free_monomers(parameters, simulation):
        '''
        Add randomly distributed actin monomers, Arp2/3 dimers, 
        and capping protein according to concentrations and box size
        '''
        ActinUtil.add_actin_monomers(
            ReaddyUtil.calculate_nParticles(
                parameters['actin_concentration'], box_size), 
            [box_size] * 3, 
            simulation
        )
        ActinUtil.add_arp23_dimers(
            ReaddyUtil.calculate_nParticles(
                parameters['arp23_concentration'], box_size), 
            [box_size] * 3, 
            simulation
        )
        ActinUtil.add_capping_protein(
            ReaddyUtil.calculate_nParticles(
                parameters['cap_concentration'], box_size), 
            [box_size] * 3, 
            simulation
        )

    def create_actin_system(self, parameters):
        '''
        Create the ReaDDy system for actin 
        including particle types, constraints, and reactions
        '''
        set_box_size(parameters['box_size'])
        system = readdy.ReactionDiffusionSystem([box_size] * 3)
        parameters['temperature_K'] = parameters['temperature_C'] + 273.15
        system.temperature = parameters['temperature_K']
        ActinUtil.add_actin_types(parameters, system)
        self.add_actin_constraints(parameters, system)
        ActinUtil.add_actin_reactions(parameters, system)
        return system

    @staticmethod
    def add_actin_types(parameters, system):
        '''
        Add particle and topology types for actin particles 
        to the ReaDDy system
        '''
        temperature = parameters['temperature_K']
        eta = parameters['eta']
        actin_diffCoeff = calculate_diffusionCoefficient(
            parameters['actin_radius'], 
            eta, 
            temperature
        ) #nm^2/s
        arp23_diffCoeff = calculate_diffusionCoefficient(
            parameters['arp23_radius'], 
            eta, 
            temperature
        ) #nm^2/s
        cap_diffCoeff = calculate_diffusionCoefficient(
            parameters['actin_radius'], 
            eta, 
            temperature
        ) #nm^2/s
        
        system.topologies.add_type("Arp23-Dimer")
        system.add_topology_species("arp2", arp23_diffCoeff)
        system.add_topology_species("arp2#ATP", arp23_diffCoeff)
        system.add_topology_species("arp2#new", arp23_diffCoeff)
        system.add_topology_species("arp2#new_ATP", arp23_diffCoeff)
        system.add_topology_species("arp3", arp23_diffCoeff)
        system.add_topology_species("arp3#branched", arp23_diffCoeff)

        system.topologies.add_type("Cap")
        system.add_topology_species("cap", cap_diffCoeff)
        system.add_topology_species("cap#new", cap_diffCoeff)
        system.add_topology_species("cap#bound", cap_diffCoeff)

        system.topologies.add_type("Actin-Monomer")
        system.topologies.add_type("Actin-Dimer")
        system.topologies.add_type("Actin-Dimer#Fail")
        system.topologies.add_type("Actin-Trimer")
        system.topologies.add_type("Actin-Trimer#Growing")
        system.topologies.add_type("Actin-Trimer#Shrinking")
        system.topologies.add_type("Actin-Trimer#Fail")
        system.topologies.add_type("Polymer")
        system.topologies.add_type("Polymer#GrowingPointed")
        system.topologies.add_type("Polymer#GrowingBarbed")
        system.topologies.add_type("Polymer#Shrinking")
        system.topologies.add_type("Polymer#Branching")
        system.topologies.add_type("Polymer#Branch-Nucleating")
        system.topologies.add_type("Polymer#Capping")
        system.topologies.add_type("Polymer#Fail-Pointed-Shrink-ATP")
        system.topologies.add_type("Polymer#Fail-Pointed-Shrink-ADP")
        system.topologies.add_type("Polymer#Fail-Barbed-Shrink-ATP")
        system.topologies.add_type("Polymer#Fail-Barbed-Shrink-ADP")
        system.topologies.add_type("Polymer#Fail-Hydrolysis-Actin")
        system.topologies.add_type("Polymer#Fail-Hydrolysis-Arp")
        system.topologies.add_type("Polymer#Fail-Branch-ATP")
        system.topologies.add_type("Polymer#Fail-Branch-ADP")
        system.topologies.add_type("Polymer#Fail-Arp-Bind-ATP")
        system.topologies.add_type("Polymer#Fail-Arp-Bind-ADP")
        system.topologies.add_type("Polymer#Fail-Debranch-ATP")
        system.topologies.add_type("Polymer#Fail-Debranch-ADP")
        system.topologies.add_type("Polymer#Fail-Arp-Unbind-ATP")
        system.topologies.add_type("Polymer#Fail-Arp-Unbind-ADP")
        system.topologies.add_type("Polymer#Fail-Nucleotide-Exchange-Actin")
        system.topologies.add_type("Polymer#Fail-Nucleotide-Exchange-Arp")
        system.topologies.add_type("Polymer#Fail-Cap-Unbind")
        system.add_topology_species("actin#free", actin_diffCoeff)
        system.add_topology_species("actin#free_ATP", actin_diffCoeff)
        system.add_topology_species("actin#new", actin_diffCoeff)
        system.add_topology_species("actin#new_ATP", actin_diffCoeff)
        for i in range(3):
            system.add_topology_species(f"actin#{i}", actin_diffCoeff)
            system.add_topology_species(f"actin#ATP_{i}", actin_diffCoeff)
            system.add_topology_species(f"actin#pointed_{i}", actin_diffCoeff)
            system.add_topology_species(f"actin#pointed_ATP_{i}", actin_diffCoeff)
            system.add_topology_species(f"actin#barbed_{i}", actin_diffCoeff)
            system.add_topology_species(f"actin#barbed_ATP_{i}", actin_diffCoeff)
        system.add_topology_species("actin#branch_1", actin_diffCoeff)
        system.add_topology_species("actin#branch_ATP_1", actin_diffCoeff)
        system.add_topology_species("actin#branch_barbed_1", actin_diffCoeff)
        system.add_topology_species("actin#branch_barbed_ATP_1", actin_diffCoeff)

    def add_actin_constraints(self, parameters, system):
        '''
        Add geometric constraints for connected actin particles, 
        including bonds, angles, and repulsions, to the ReaDDy system
        '''
        force_constant = parameters['force_constant']

        self.add_bonds_between_actins(force_constant, system)
        self.add_branch_bonds(force_constant, system)
        self.add_temporary_bonds(force_constant, system)

        self.add_filament_twist_angles(force_constant, system)
        self.add_filament_twist_dihedrals(force_constant, system)
        self.add_branch_angles(force_constant, system)
        self.add_branch_dihedrals(force_constant, system)

        self.add_cap_bonds(force_constant, system)
        self.add_cap_angles(force_constant, system)
        self.add_cap_dihedrals(force_constant, system)

        self.add_repulsions(force_constant, system)

    @staticmethod
    def add_actin_reactions(parameters, system):
        '''
        Add reactions to the ReaDDy system
        '''
        actin_radius = parameters['actin_radius']
        reaction_distance = parameters['reaction_distance']
        ActinUtil.add_dimerize_reaction(
            system, 
            parameters['dimerize_rate'], 
            2 * actin_radius + reaction_distance
        )
        ActinUtil.add_dimerize_reverse_reaction(system)
        ActinUtil.add_trimerize_reaction(
            system, 
            parameters['trimerize_rate'], 
            2 * actin_radius + reaction_distance
        )
        ActinUtil.add_trimerize_reverse_reaction(system)
        ActinUtil.add_nucleate_reaction(
            system, 
            parameters['nucleate_ATP_rate'], 
            parameters['nucleate_ADP_rate'],
            2 * actin_radius + reaction_distance
        )
        ActinUtil.add_pointed_growth_reaction(
            system, 
            parameters['pointed_growth_ATP_rate'], 
            parameters['pointed_growth_ADP_rate'],
            2 * actin_radius + reaction_distance
        )
        ActinUtil.add_pointed_shrink_reaction(system)
        ActinUtil.add_barbed_growth_reaction(
            system, 
            parameters['barbed_growth_ATP_rate'], 
            parameters['barbed_growth_ADP_rate'],
            2 * actin_radius + reaction_distance
        )
        ActinUtil.add_barbed_shrink_reaction(system)
        ActinUtil.add_hydrolyze_reaction(system)
        ActinUtil.add_actin_nucleotide_exchange_reaction(system)
        ActinUtil.add_arp23_nucleotide_exchange_reaction(system)
        ActinUtil.add_arp23_bind_reaction(
            system, 
            parameters['arp_bind_ATP_rate'], 
            parameters['arp_bind_ADP_rate'],
            actin_radius + parameters['arp23_radius'] + reaction_distance
        )
        ActinUtil.add_arp23_unbind_reaction(system)
        ActinUtil.add_nucleate_branch_reaction(
            system, 
            parameters['barbed_growth_branch_ATP_rate'], 
            parameters['barbed_growth_branch_ADP_rate'],
            actin_radius + parameters['arp23_radius'] + reaction_distance
        )
        ActinUtil.add_debranch_reaction(system)
        ActinUtil.add_cap_bind_reaction(
            system, 
            parameters['cap_bind_rate'], 
            actin_radius + parameters['cap_radius'] + reaction_distance
        )
        ActinUtil.add_cap_unbind_reaction(system)

    @staticmethod
    def get_new_arp23(topology):
        '''
        get an arp2 and its unbranched arp3 neighbor,
        meaning the arp2/3 dimer has just bound
        '''
        for vertex in topology.graph.get_vertices():
            pt = topology.particle_type_of_vertex(vertex)
            if "arp2#new" in pt:
                for neighbor in vertex:
                    if topology.particle_type_of_vertex(neighbor.get()) == "arp3":
                        return vertex, neighbor.get()
        return None, None

    @staticmethod
    def get_first_actin_neighbor(topology, vertex, last_actin):
        '''
        get the first neighboring vertex that is some type of actin,
        excluding last_actin
        '''
        for neighbor in vertex:
            v_neighbor = neighbor.get()
            pt = topology.particle_type_of_vertex(v_neighbor)
            if v_neighbor != last_actin and "actin" in pt:
                return v_neighbor
        return None

    @staticmethod
    def cancel_failed_branch_reaction(topology, recipe, actin_arp2, arp2):
        '''
        Undo the branching spatial reaction if the structural reaction fails
        '''
        if VERBOSE:
            print("Canceling branch reaction")
        pt = topology.particle_type_of_vertex(actin_arp2)
        recipe.remove_edge(actin_arp2, arp2)
        ReaddyUtil.set_flags(topology, recipe, arp2, [], ["new"], True)
        recipe.change_topology_type(
            "Polymer#Fail-Branch-{}".format("ATP" if ("ATP" in pt) else "ADP"))

    @staticmethod
    def get_actin_number(topology, vertex, offset):
        '''
        get the type number for an actin plus the given offset in range [-1, 1]
        (i.e. return 3 for type = "actin#ATP_1" and offset = -1)
        '''
        pt = topology.particle_type_of_vertex(vertex)
        if not "actin" in pt:
            raise Error("Failed to get actin number: {} is not actin".format(pt))

        return ReaddyUtil.calculate_polymer_number(int(pt[-1]), offset)

    @staticmethod
    def get_actin_number_difference(topology, vertex1, vertex2):
        '''
        get the difference in actin number for two adjacent actins
        '''
        n1 = ActinUtil.get_actin_number(topology, vertex1, 0)
        n2 = ActinUtil.get_actin_number(topology, vertex2, 0)

        d = n2 - n1
        if abs(d) > 1:
            d *= -1 / abs(d)

        return d

    @staticmethod
    def get_all_polymer_actin_types(vertex_type):
        '''
        get a list of all numbered versions of a type
        (e.g. for "actin#ATP" return ["actin#ATP_1", "actin#ATP_2", "actin#ATP_3"])
        '''
        spacer = "_"
        if not "#" in vertex_type:
            spacer = "#"

        return ["{}{}1".format(vertex_type, spacer),
                "{}{}2".format(vertex_type, spacer),
                "{}{}3".format(vertex_type, spacer)]

    @staticmethod
    def get_actin_orientation(positions):
        '''
        orthonormalize and cross the vectors to an actin's two neighbor actins
        to get a basis local to the actin,
        positions = [previous actin position, this actin position, next actin position]
        '''
        v1 = ReaddyUtil.normalize(positions[0] - positions[1])
        v2 = ReaddyUtil.normalize(positions[2] - positions[1])
        v2 = ReaddyUtil.normalize(v2 - (np.dot(v1, v2) / np.dot(v1, v1)) * v1)
        v3 = np.cross(v2, v1)

        return np.matrix([[v1[0], v2[0], v3[0]],
                        [v1[1], v2[1], v3[1]],
                        [v1[2], v2[2], v3[2]]])

    @staticmethod
    def get_actin_rotation(positions):
        '''
        get the difference in the actin's current orientation
        compared to the initial orientation as a rotation matrix
        positions = [previous actin position, middle actin position, next actin position]
        '''
        positions[0] = ReaddyUtil.get_non_periodic_boundary_position(
            positions[1], positions[0], box_size)
        positions[2] = ReaddyUtil.get_non_periodic_boundary_position(
            positions[1], positions[2], box_size)

        current_orientation = ActinUtil.get_actin_orientation(positions)

        return np.matmul(current_orientation, np.linalg.inv(INITIAL_ORIENTATION))

    @staticmethod
    def get_position_for_new_vertex(neighbor_positions, offset_vector):
        '''
        get the offset vector in the local space for the actin at neighbor_positions[1]
        neighbor_positions = [
            previous actin position,
            middle actin position,
            next actin position
        ]
        '''
        rotation = ActinUtil.get_actin_rotation(neighbor_positions)
        if rotation is None:
            return None

        vector_to_new_pos = np.squeeze(np.array(np.dot(rotation, offset_vector)))

        return (neighbor_positions[1] + vector_to_new_pos).tolist()

    @staticmethod
    def get_next_arp3(topology, vertex, last_vertex_id, max_edges):
        '''
        recurse along the chain until arp3 is found or max_edges is reached
        '''
        for neighbor in vertex:

            n_id = topology.particle_id_of_vertex(neighbor)
            if n_id == last_vertex_id:
                continue

            pt = topology.particle_type_of_vertex(neighbor)
            if pt == "arp3#branched" or pt == "arp3":
                return neighbor.get(), max_edges
            else:
                if max_edges <= 1:
                    return None, max_edges
                return ActinUtil.get_next_arp3(topology, neighbor.get(), n_id, max_edges-1)

        return None, max_edges

    @staticmethod
    def get_branch_orientation_vertices_and_offset(topology, vertex):
        '''
        get orientation vertices [actin, actin_arp2, actin_arp3]
        for a new actin within 3 actins of a branch,
        as well as the offset vector
        '''
        v_arp3, edges = ActinUtil.get_next_arp3(topology, vertex, None, 4)
        if v_arp3 is None:
            if VERBOSE:
                print("couldn't set position of new vertex: failed to find arp3")
            return None, None
        offset_index = 4 - edges

        v_actin_arp3 = ReaddyUtil.get_neighbor_of_types(
            topology, v_arp3, ["actin#arp3", "actin#arp3_ATP"], [])
        if v_actin_arp3 is None:
            if VERBOSE:
                print("failed to find actin_arp3 so using arp3 to position vertex instead")
            v_actin_arp3 = v_arp3

        v_arp2 = ReaddyUtil.get_neighbor_of_types(
            topology, v_arp3, ["arp2", "arp2#ATP"], [])
        if v_arp2 is None:
            if VERBOSE:
                print("couldn't set position of new vertex: failed to find arp2")
            return None, None

        v_actin_arp2 = ReaddyUtil.get_neighbor_of_types(
            topology, v_arp2, ["actin#arp2", "actin#arp2_ATP"], [])
        if v_actin_arp2 is None:
            if VERBOSE:
                print("couldn't set position of new vertex: failed to find actin_arp2")
            return None, None

        actin_types = (ActinUtil.get_all_polymer_actin_types("actin")
            + ActinUtil.get_all_polymer_actin_types("actin#ATP")
            + ActinUtil.get_all_polymer_actin_types("actin#pointed")
            + ActinUtil.get_all_polymer_actin_types("actin#pointed_ATP")
            + ActinUtil.get_all_polymer_actin_types("actin#barbed")
            + ActinUtil.get_all_polymer_actin_types("actin#barbed_ATP")
            + ActinUtil.get_all_polymer_actin_types("actin#branch")
            + ActinUtil.get_all_polymer_actin_types("actin#branch_ATP"))

        v1 = ReaddyUtil.get_neighbor_of_types(
            topology, v_actin_arp2, actin_types, [v_actin_arp3])
        if v1 is None:
            if VERBOSE:
                print("couldn't set position of new vertex: failed to find v1")
            return None, None

        return ([v1, v_actin_arp2, v_actin_arp3],
            VECTOR_TO_NEW_BRANCH_ACTIN[offset_index])

    @staticmethod
    def set_end_vertex_position(topology, recipe, v_new, barbed):
        '''
        set the position of a new pointed or barbed vertex
        '''
        actin_types = (ActinUtil.get_all_polymer_actin_types("actin")
            + ActinUtil.get_all_polymer_actin_types("actin#ATP")
            + ActinUtil.get_all_polymer_actin_types("actin#pointed")
            + ActinUtil.get_all_polymer_actin_types("actin#pointed_ATP")
            + ActinUtil.get_all_polymer_actin_types("actin#barbed")
            + ActinUtil.get_all_polymer_actin_types("actin#barbed_ATP")
            + ActinUtil.get_all_polymer_actin_types("actin#branch")
            + ActinUtil.get_all_polymer_actin_types("actin#branch_ATP")
            + ActinUtil.get_all_polymer_actin_types("actin#branch_barbed")
            + ActinUtil.get_all_polymer_actin_types("actin#branch_barbed_ATP"))

        vertices = []
        offset_vector = VECTOR_TO_NEW_BARBED if barbed else VECTOR_TO_NEW_POINTED
        at_branch = False

        vertices.append(ReaddyUtil.get_neighbor_of_types(
            topology, v_new, actin_types, []))
        if vertices[0] is None:

            vertices, offset_vector = ActinUtil.get_branch_orientation_vertices_and_offset(
                topology, v_new)
            if vertices is None:
                return
            at_branch = True

        else:
            vertices.append(ReaddyUtil.get_neighbor_of_types(
                topology, vertices[0], actin_types, [v_new]))
            if vertices[1] is None:

                vertices, offset_vector = ActinUtil.get_branch_orientation_vertices_and_offset(
                    topology, v_new)
                if vertices is None:
                    return
                at_branch = True

            else:
                vertices.append(ReaddyUtil.get_neighbor_of_types(
                    topology, vertices[1], actin_types, [vertices[0]]))
                if vertices[2] is None:

                    vertices, offset_vector = ActinUtil.get_branch_orientation_vertices_and_offset(
                        topology, v_new)
                    if vertices is None:
                        return
                    at_branch = True

        positions = []
        for v in vertices:
            positions.append(ReaddyUtil.get_vertex_position(topology, v))

        if barbed and not at_branch:
            positions = positions[::-1]

        pos = ActinUtil.get_position_for_new_vertex(positions, offset_vector)
        if pos is not None:
            recipe.change_particle_position(v_new, pos)

    @staticmethod
    def set_new_trimer_vertex_position(topology, recipe,
                                    v_new, v_pointed, v_barbed):
        '''
        set the position of an actin monomer just added to a dimer to create a trimer
        '''
        pos_new = ReaddyUtil.get_vertex_position(topology, v_new)
        pos_pointed = ReaddyUtil.get_vertex_position(topology, v_pointed)
        pos_barbed = ReaddyUtil.get_vertex_position(topology, v_barbed)

        v_barbed_to_pointed = pos_pointed - pos_barbed
        v_barbed_to_new = pos_new - pos_barbed
        current_angle = ReaddyUtil.get_angle_between_vectors(
            v_barbed_to_pointed, v_barbed_to_new)
        angle = ANGLE_BETWEEN_ACTINS - current_angle
        axis = np.cross(v_barbed_to_pointed, v_barbed_to_new)
        pos = pos_barbed + ReaddyUtil.rotate(
            4.27 * ReaddyUtil.normalize(v_barbed_to_new), axis, angle)

        recipe.change_particle_position(v_new, pos)

    @staticmethod
    def set_arp23_vertex_position(topology, recipe, v_arp2, v_arp3,
                                v_actin_arp2, v_actin_arp3):
        '''
        set the position of new arp2/3 vertices
        '''
        actin_types = (ActinUtil.get_all_polymer_actin_types("actin")
            + ActinUtil.get_all_polymer_actin_types("actin#ATP")
            + ActinUtil.get_all_polymer_actin_types("actin#pointed")
            + ActinUtil.get_all_polymer_actin_types("actin#pointed_ATP")
            + ActinUtil.get_all_polymer_actin_types("actin#branch")
            + ActinUtil.get_all_polymer_actin_types("actin#branch_ATP"))

        v1 = ReaddyUtil.get_neighbor_of_types(
            topology, v_actin_arp2, actin_types, [v_actin_arp3])
        if v1 is None:
            if VERBOSE:
                print("couldn't set position of new vertex: failed to find v1")
            return

        pos1 = ReaddyUtil.get_vertex_position(topology, v1)
        pos2 = ReaddyUtil.get_vertex_position(topology, v_actin_arp2)
        pos3 = ReaddyUtil.get_vertex_position(topology, v_actin_arp3)

        pos_arp2 = ActinUtil.get_position_for_new_vertex(
            [pos1, pos2, pos3], VECTOR_TO_NEW_ARP2)
        if pos_arp2 is not None:
            recipe.change_particle_position(v_arp2, pos_arp2)

        pos_arp3 = ActinUtil.get_position_for_new_vertex(
            [pos1, pos2, pos3], VECTOR_TO_NEW_ARP3)
        if pos_arp3 is not None:
            recipe.change_particle_position(v_arp3, pos_arp3)

    @staticmethod
    def get_random_arp3(topology, with_ATP, with_branch):
        '''
        get a random bound arp3 with the given arp2 nucleotide state
        and with or without a branch attached to the arp3
        '''
        v_arp2s = ReaddyUtil.get_vertices_of_type(
            topology, "arp2#ATP" if with_ATP else "arp2", True)
        if len(v_arp2s) < 1:
            if VERBOSE:
                print("failed to find arp2 with {}".format(
                    "ATP" if with_ATP else "ADP"))
            return None

        v_arp3s = []
        for v_arp2 in v_arp2s:
            v_arp3 = ReaddyUtil.get_neighbor_of_types(
                topology, v_arp2, ["arp3#branched" if with_branch else "arp3"], [])
            if v_arp3 is not None:
                v_arp3s.append(v_arp3)
        if len(v_arp3s) < 1:
            if VERBOSE:
                print("failed to find arp3 with{} branch".format(
                    "out" if not with_branch else ""))
            return None

        return random.choice(v_arp3s)

    @staticmethod
    def add_branched_seed(simulation):
        '''
        add a branched seed
        '''
        positions = np.array([[15, 0, 0],[16.31, -2.9, -2.85],[15.7, 0.22, -5.7],
            [15.57, -2.96, -8.55],[16.43, 0.11, -11.4],[14.89, -2.68, -14.25],
            [17.03, -0.32, -17.1],[14.4, -2.12, -19.95],[17.38, -0.98, -22.8],
            [14.22, -1.4, -25.65],[17.39, -1.72, -28.5],[14.38, -0.68, -31.35],
            [12.4, -5.71, -19.28],[15.13, -8.57, -20.61],[12.26, -11.19, -19.05],
            [15.37, -14.1, -19.49],[12.18, -16.94, -19.25]])
        types = ["actin#pointed_ATP_1", "actin#ATP_2", "actin#ATP_3",
                "actin#ATP_1", "actin#ATP_2", "actin#ATP_3",
                "actin#ATP_1", "actin#ATP_2", "actin#ATP_3",
                "actin#ATP_1", "actin#ATP_2", "actin#barbed_ATP_3",
                "arp2#ATP", "arp3#branched", "actin#branch_ATP_1",
                "actin#ATP_2", "actin#barbed_ATP_3"]

        top = simulation.add_topology("Polymer", types, positions)

        for i in range(1,len(positions)):
            if i != 12:
                top.get_graph().add_edge(i-1, i)
        top.get_graph().add_edge(7, 12)
        top.get_graph().add_edge(8, 13)

    @staticmethod
    def add_seed(simulation):
        '''
        add a 6-mer seed
        '''
        positions = np.array([[15, 0, 0],[16.31, -2.9, -2.85],[15.7, 0.22, -5.7],
            [15.57, -2.96, -8.55],[16.43, 0.11, -11.4],[14.89, -2.68, -14.25]])
        types = ["actin#pointed_ATP_1", "actin#ATP_2", "actin#ATP_3",
            "actin#ATP_1", "actin#ATP_2", "actin#barbed_ATP_3"]

        top = simulation.add_topology("Polymer", types, positions)

        for i in range(1,len(positions)):
            top.get_graph().add_edge(i-1, i)

    @staticmethod
    def add_dimer(position, simulation):
        '''
        add an actin dimer seed
        '''
        positions = np.array(
            [[0, 0, 0], 4.27 * ReaddyUtil.get_random_unit_vector()])
        types = ["actin#pointed_ATP_1", "actin#barbed_ATP_2"]

        top = simulation.add_topology("Actin-Dimer", types, position + positions)
        top.get_graph().add_edge(0, 1)

    @staticmethod
    def add_dimers(n, box_dims, simulation):
        '''
        add actin dimers
        '''
        positions = np.random.uniform(size=(n,3)) * box_dims - box_dims * 0.5
        for p in range(len(positions)):
            ActinUtil.add_dimer(positions[p], simulation)

    @staticmethod
    def add_actin_monomers(n, box_dims, simulation):
        '''
        add free actin
        '''
        positions = np.random.uniform(size=(n,3)) * box_dims - box_dims * 0.5
        for p in range(len(positions)):
            simulation.add_topology(
                "Actin-Monomer", ["actin#free_ATP"], np.array([positions[p]]))

    @staticmethod
    def add_arp23_dimers(n, box_dims, simulation):
        '''
        add arp2/3 dimers
        '''
        positions = np.random.uniform(size=(n,3)) * box_dims - box_dims * 0.5
        for p in range(len(positions)):
            top = simulation.add_topology("Arp23-Dimer", ["arp2#ATP", "arp3"],
                np.array([positions[p], positions[p]
                + 4. * ReaddyUtil.get_random_unit_vector()]))
            top.get_graph().add_edge(0, 1)

    @staticmethod
    def add_capping_protein(n, box_dims, simulation):
        '''
        add free capping protein
        '''
        positions = np.random.uniform(size=(n,3)) * box_dims - box_dims * 0.5
        for p in range(len(positions)):
            simulation.add_topology("Cap", ["cap"], np.array([positions[p]]))

    @staticmethod
    def do_shrink(topology, recipe, barbed, ATP):
        '''
        remove an (ATP or ADP)-actin from the (barbed or pointed) end
        '''
        end_type = "actin#{}{}".format(
            "barbed" if barbed else "pointed", "_ATP" if ATP else "")
        v_end = ReaddyUtil.get_random_vertex_of_types(
            topology, get_types_of_all_numbers(end_type))
        if v_end is None:
            if VERBOSE:
                print("failed to find end actin to remove")
            return False

        v_arp = ReaddyUtil.get_neighbor_of_types(
            topology, v_end, ["arp2", "arp2#ATP", "arp3", "arp3#branched"], [])
        if v_arp is not None:
            if VERBOSE:
                print("failed to remove actin because a branch was attached")
            return False

        v_neighbor = ReaddyUtil.get_neighbor_of_types(
            topology, v_end, get_types_of_all_numbers("actin")
            + get_types_of_all_numbers("actin#ATP"), [])

        if v_neighbor is None:
            if VERBOSE:
                print("failed to find plain actin neighbor of actin to remove")
            return False

        recipe.remove_edge(v_end, v_neighbor)
        recipe.change_particle_type(
            v_end, "actin#free" if not ATP else "actin#free_ATP")
        ReaddyUtil.set_flags(topology, recipe, v_neighbor,
            ["barbed"] if barbed else ["pointed"], [], True)
        recipe.change_topology_type("Polymer#Shrinking")

        return True

    @staticmethod
    def do_arp23_unbind(topology, recipe, with_ATP):
        '''
        dissociate an arp2/3 from a mother filament
        '''
        v_arp3 = ActinUtil.get_random_arp3(topology, with_ATP, False)
        if v_arp3 is None:
            recipe.change_topology_type("Polymer#Fail-Arp-Unbind-{}".format(
                "ATP" if with_ATP else "ADP"))
            return recipe

        actin_types = (get_types_of_all_numbers("actin")
            + get_types_of_all_numbers("actin#ATP")
            + get_types_of_all_numbers("actin#pointed")
            + get_types_of_all_numbers("actin#pointed_ATP")
            + get_types_of_all_numbers("actin#barbed")
            + get_types_of_all_numbers("actin#barbed_ATP"))

        v_actin_arp3 = ReaddyUtil.get_neighbor_of_types(
            topology, v_arp3, actin_types, [])
        if v_actin_arp3 is None:
            if VERBOSE:
                print("failed to find actin_arp3")
            recipe.change_topology_type("Polymer#Fail-Arp-Unbind-{}".format(
                "ATP" if with_ATP else "ADP"))
            return recipe

        v_arp2 = ReaddyUtil.get_neighbor_of_types(
            topology, v_arp3, ["arp2", "arp2#ATP"], [])
        if v_arp2 is None:
            if VERBOSE:
                print("failed to find arp2")
            recipe.change_topology_type("Polymer#Fail-Arp-Unbind-{}".format(
                "ATP" if with_ATP else "ADP"))
            return recipe

        v_actin_arp2 = ReaddyUtil.get_neighbor_of_types(
            topology, v_arp2, actin_types, [])
        if v_actin_arp2 is None:
            if VERBOSE:
                print("failed to find actin_arp2")
            recipe.change_topology_type("Polymer#Fail-Arp-Unbind-{}".format(
                "ATP" if with_ATP else "ADP"))
            return recipe

        recipe.remove_edge(v_arp2, v_actin_arp2)
        recipe.remove_edge(v_arp3, v_actin_arp3)
        recipe.change_topology_type("Polymer#Shrinking")

    @staticmethod
    def do_debranching(topology, recipe, with_ATP):
        '''
        reaction function to detach a branch filament from arp2/3
        '''
        v_arp3 = ActinUtil.get_random_arp3(topology, with_ATP, True)
        if v_arp3 is None:
            recipe.change_topology_type("Polymer#Fail-Debranch-{}".format(
                "ATP" if with_ATP else "ADP"))
            return recipe

        actin_types = (get_types_of_all_numbers("actin#branch")
                    + get_types_of_all_numbers("actin#branch_ATP")
                    + get_types_of_all_numbers("actin#branch_barbed")
                    + get_types_of_all_numbers("actin#branch_barbed_ATP"))

        v_actin1 = ReaddyUtil.get_neighbor_of_types(
            topology, v_arp3, actin_types, [])
        if v_actin1 is None:
            if VERBOSE:
                print("failed to find first branch actin")
            recipe.change_topology_type("Polymer#Fail-Debranch-{}".format(
                "ATP" if with_ATP else "ADP"))
            return recipe

        recipe.remove_edge(v_arp3, v_actin1)
        ReaddyUtil.set_flags(topology, recipe, v_arp3, [], ["branched"], True)
        pt_actin1 = topology.particle_type_of_vertex(v_actin1)
        if "barbed" in pt_actin1: # branch is a monomer
            recipe.change_particle_type(v_actin1, "actin#free{}".format(
                "_ATP" if "ATP" in pt_actin1 else ""))
        else:
            ReaddyUtil.set_flags(topology, recipe, v_actin1,
                ["pointed"], ["branch"], True) # branch is a filament
        recipe.change_topology_type("Polymer#Shrinking")

    def add_bonds_between_actins(self, force_constant, system):
        '''
        add bonds between actins
        '''
        self.readdy_util.add_bond(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2"],
            force_constant, 4.27, system
        )
        self.readdy_util.add_bond(
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            force_constant, 4.27, system
        )
        self.readdy_util.add_bond(
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1"],
            force_constant, 4.27, system
        )

    def add_branch_bonds(self, force_constant, system):
        '''
        add bonds between arp2, arp3, and actins
        '''
        self.readdy_util.add_bond( # mother filament to arp2 bonds
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1",
            "actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2",
            "actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["arp2", "arp2#ATP", "arp2#new", "arp2#new_ATP"],
            force_constant, 4.17, system
        )
        self.readdy_util.add_bond( # mother filament to arp3 bonds
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1",
            "actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2",
            "actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            ["arp3", "arp3#branched"],
            force_constant, 8.22, system
        )
        self.readdy_util.add_bond( # arp2 to arp3 bonds
            ["arp2", "arp2#ATP", "arp2#new", "arp2#new_ATP"],
            ["arp3", "arp3#branched"],
            force_constant, 4.18, system
        )
        self.readdy_util.add_bond( # arp3 to daughter filament bonds
            ["arp3#branched"],
            ["actin#branch_1", "actin#branch_ATP_1",
            "actin#branch_barbed_1", "actin#branch_barbed_ATP_1"],
            force_constant, 4.19, system
        )

    def add_temporary_bonds(self, force_constant, system):
        '''
        add temporary bonds during growth reactions
        '''
        self.readdy_util.add_bond(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_barbed_1", "actin#branch_barbed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1",
            "actin#barbed_1", "actin#barbed_ATP_1",
            "actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2",
            "actin#barbed_2", "actin#barbed_ATP_2",
            "actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3",
            "actin#barbed_3", "actin#barbed_ATP_3", "arp3#branched"],
            ["actin#new", "actin#new_ATP"],
            force_constant, 4.27, system
        )
        self.readdy_util.add_bond(
            ["actin#1", "actin#ATP_1", "actin#2", "actin#ATP_2",
            "actin#3", "actin#ATP_3"],
            ["arp2", "arp2#ATP"],
            force_constant, 4.19, system
        )

    def add_filament_twist_angles(self, force_constant, system):
        '''
        add angles for filament twist and cohesiveness
        '''
        self.readdy_util.add_angle(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            force_constant, 1.48, system
        )
        self.readdy_util.add_angle(
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1"],
            force_constant, 1.48, system
        )
        self.readdy_util.add_angle(
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2"],
            force_constant, 1.48, system
        )

    def add_filament_twist_dihedrals(self, force_constant, system):
        '''
        add dihedrals for filament twist and cohesiveness
        '''
        self.readdy_util.add_dihedral(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1"],
            force_constant, 2.80, system
        )
        self.readdy_util.add_dihedral(
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2"],
            force_constant, 2.80, system
        )
        self.readdy_util.add_dihedral(
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            force_constant, 2.80, system
        )

    def add_branch_angles(self, force_constant, system):
        '''
        add angles for branching
        '''
        self.readdy_util.add_angle(
            ["arp2", "arp2#ATP"],
            ["arp3#branched"],
            ["actin#branch_1", "actin#branch_ATP_1",
            "actin#branch_barbed_1", "actin#branch_barbed_ATP_1"],
            force_constant, 1.43, system
        )
        self.readdy_util.add_angle(
            ["arp3", "arp3#branched"],
            ["actin#1", "actin#ATP_1"],
            ["actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2"],
            force_constant, 1.46, system
        )
        self.readdy_util.add_angle(
            ["arp3", "arp3#branched"],
            ["actin#2", "actin#ATP_2"],
            ["actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            force_constant, 1.46, system
        )
        self.readdy_util.add_angle(
            ["arp3", "arp3#branched"],
            ["actin#3", "actin#ATP_3"],
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1"],
            force_constant, 1.46, system
        )
        self.readdy_util.add_angle(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["arp2", "arp2#ATP"],
            force_constant, 2.15, system
        )
        self.readdy_util.add_angle(
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["arp2", "arp2#ATP"],
            force_constant, 2.15, system
        )
        self.readdy_util.add_angle(
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["arp2", "arp2#ATP"],
            force_constant, 2.15, system
        )
        self.readdy_util.add_angle(
            ["actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2"],
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1"],
            ["arp2", "arp2#ATP"],
            force_constant, 2.31, system
        )
        self.readdy_util.add_angle(
            ["actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["arp2", "arp2#ATP"],
            force_constant, 2.31, system
        )
        self.readdy_util.add_angle(
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1"],
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["arp2", "arp2#ATP"],
            force_constant, 2.31, system
        )
        self.readdy_util.add_angle(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1"],
            ["actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2"],
            ["arp3", "arp3#branched"],
            force_constant, 0.91, system
        )
        self.readdy_util.add_angle(
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            ["arp3", "arp3#branched"],
            force_constant, 0.91, system
        )
        self.readdy_util.add_angle(
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1"],
            ["arp3", "arp3#branched"],
            force_constant, 0.91, system
        )
        self.readdy_util.add_angle(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2",
            "actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["arp2", "arp2#ATP"],
            ["arp3", "arp3#branched"],
            force_constant, 1.80, system
        )
        self.readdy_util.add_angle(
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1",
            "actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2",
            "actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            ["arp3", "arp3#branched"],
            ["arp2", "arp2#ATP"],
            force_constant, 1.19, system
        )

    def add_branch_dihedrals(self, force_constant, system):
        '''
        add dihedrals for branching
        '''
        self.readdy_util.add_dihedral(
            ["arp2", "arp2#ATP"],
            ["arp3#branched"],
            ["actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2"],
            force_constant, 2.92, system
        )
        self.readdy_util.add_dihedral(
            ["arp3#branched"],
            ["actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            force_constant, 2.82, system
        )
        self.readdy_util.add_dihedral(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2",
            "actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3",],
            ["arp2", "arp2#ATP"],
            ["arp3#branched"],
            ["actin#branch_1", "actin#branch_ATP_1",
            "actin#branch_barbed_1", "actin#branch_barbed_ATP_1"],
            force_constant, 2.76, system
        )
        self.readdy_util.add_dihedral(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["arp2", "arp2#ATP"],
            ["arp3", "arp3#branched"],
            force_constant, 1.67, system
        )
        self.readdy_util.add_dihedral(
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["arp2", "arp2#ATP"],
            ["arp3", "arp3#branched"],
            force_constant, 1.67, system
        )
        self.readdy_util.add_dihedral(
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["arp2", "arp2#ATP"],
            ["arp3", "arp3#branched"],
            force_constant, 1.67, system
        )
        self.readdy_util.add_dihedral(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["arp2", "arp2#ATP"],
            force_constant, 0.60, system
        )
        self.readdy_util.add_dihedral(
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["arp2", "arp2#ATP"],
            force_constant, 0.60, system
        )
        self.readdy_util.add_dihedral(
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["arp2", "arp2#ATP"],
            force_constant, 0.60, system
        )
        self.readdy_util.add_dihedral(
            ["actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            ["actin#2", "actin#ATP_2"],
            ["arp3#branched"],
            ["actin#branch_1", "actin#branch_ATP_1",
            "actin#branch_barbed_1", "actin#branch_barbed_ATP_1"],
            force_constant, 1.25, system
        )
        self.readdy_util.add_dihedral(
            ["actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2"],
            ["actin#1", "actin#ATP_1"],
            ["arp3#branched"],
            ["actin#branch_1", "actin#branch_ATP_1",
            "actin#branch_barbed_1", "actin#branch_barbed_ATP_1"],
            force_constant, 1.25, system
        )
        self.readdy_util.add_dihedral(
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1"],
            ["actin#3", "actin#ATP_3"],
            ["arp3#branched"],
            ["actin#branch_1", "actin#branch_ATP_1",
            "actin#branch_barbed_1", "actin#branch_barbed_ATP_1"],
            force_constant, 1.25, system
        )
        self.readdy_util.add_dihedral(
            ["actin#3", "actin#ATP_3", "actin#barbed_3", "actin#barbed_ATP_3"],
            ["actin#2", "actin#ATP_2"],
            ["actin#1", "actin#ATP_1"],
            ["arp3", "arp3#branched"],
            force_constant, 1.89, system
        )
        self.readdy_util.add_dihedral(
            ["actin#2", "actin#ATP_2", "actin#barbed_2", "actin#barbed_ATP_2"],
            ["actin#1", "actin#ATP_1"],
            ["actin#3", "actin#ATP_3"],
            ["arp3", "arp3#branched"],
            force_constant, 1.89, system
        )
        self.readdy_util.add_dihedral(
            ["actin#1", "actin#ATP_1", "actin#barbed_1", "actin#barbed_ATP_1"],
            ["actin#3", "actin#ATP_3"],
            ["actin#2", "actin#ATP_2"],
            ["arp3", "arp3#branched"],
            force_constant, 1.89, system
        )

    def add_cap_bonds(self, force_constant, system):
        '''
        add capping protein to actin bonds
        '''
        self.readdy_util.add_bond(
            ["actin#1", "actin#ATP_1", "actin#2", "actin#ATP_2",
            "actin#3", "actin#ATP_3"],
            ["cap#bound", "cap#new"],
            force_constant, 4.27, system
        )

    def add_cap_angles(self, force_constant, system):
        '''
        add angles for capping protein
        '''
        self.readdy_util.add_angle(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["cap#bound"],
            force_constant, 1.48, system
        )
        self.readdy_util.add_angle(
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["cap#bound"],
            force_constant, 1.48, system
        )
        self.readdy_util.add_angle(
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["cap#bound"],
            force_constant, 1.48, system
        )

    def add_cap_dihedrals(self, force_constant, system):
        '''
        add dihedrals for capping protein
        '''
        self.readdy_util.add_dihedral(
            ["actin#1", "actin#ATP_1", "actin#pointed_1", "actin#pointed_ATP_1",
            "actin#branch_1", "actin#branch_ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["cap#bound"],
            force_constant, 2.80, system
        )
        self.readdy_util.add_dihedral(
            ["actin#2", "actin#ATP_2", "actin#pointed_2", "actin#pointed_ATP_2"],
            ["actin#3", "actin#ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["cap#bound"],
            force_constant, 2.80, system
        )
        self.readdy_util.add_dihedral(
            ["actin#3", "actin#ATP_3", "actin#pointed_3", "actin#pointed_ATP_3"],
            ["actin#1", "actin#ATP_1"],
            ["actin#2", "actin#ATP_2"],
            ["cap#bound"],
            force_constant, 2.80, system
        )
        self.readdy_util.add_dihedral(
            ["arp2", "arp2#ATP"],
            ["arp3#branched"],
            ["actin#branch_1", "actin#branch_ATP_1"],
            ["cap#bound"],
            force_constant, 2.92, system
        )

    def add_repulsions(self, force_constant, system):
        '''
        add repulsions
        '''
        self.readdy_util.add_repulsion(
            ["actin#pointed_1", "actin#pointed_ATP_1", "actin#pointed_2",
            "actin#pointed_ATP_2", "actin#pointed_3", "actin#pointed_ATP_3",
            "actin#1", "actin#ATP_1", "actin#2", "actin#ATP_2",
            "actin#3", "actin#ATP_3",
            "actin#branch_1", "actin#branch_ATP_1",
            "actin#branch_barbed_1", "actin#branch_barbed_ATP_1",
            "actin#barbed_1", "actin#barbed_ATP_1", "actin#barbed_2",
            "actin#barbed_ATP_2", "actin#barbed_3", "actin#barbed_ATP_3",
            "arp2", "arp2#ATP", "arp3", "arp3#branched", "cap", "cap#bound",
            "actin#free", "actin#free_ATP"],
            ["actin#pointed_1", "actin#pointed_ATP_1", "actin#pointed_2",
            "actin#pointed_ATP_2", "actin#pointed_3", "actin#pointed_ATP_3",
            "actin#1", "actin#ATP_1", "actin#2", "actin#ATP_2",
            "actin#3", "actin#ATP_3",
            "actin#branch_1", "actin#branch_ATP_1",
            "actin#branch_barbed_1", "actin#branch_barbed_ATP_1",
            "actin#barbed_1", "actin#barbed_ATP_1", "actin#barbed_2",
            "actin#barbed_ATP_2", "actin#barbed_3", "actin#barbed_ATP_3",
            "arp2", "arp2#ATP", "arp3", "arp3#branched", "cap", "cap#bound",
            "actin#free", "actin#free_ATP"],
            force_constant, 4., system
        )

    @staticmethod
    def add_dimerize_reaction(system, rate, reaction_distance):
        '''
        attach two monomers
        '''
        system.topologies.add_spatial_reaction(
            "Dimerize: Actin-Monomer(actin#free_ATP) + {}".format(
            "Actin-Monomer(actin#free_ATP) -> {}".format(
            "Actin-Dimer(actin#pointed_ATP_1--actin#barbed_ATP_2)")),
            rate=rate, radius=reaction_distance
        )

    @staticmethod
    def add_dimerize_reverse_reaction(system):
        '''
        detach two monomers
        '''
        system.topologies.add_structural_reaction(
            "Reverse_Dimerize",
            topology_type="Actin-Dimer",
            reaction_function=reaction_function_reverse_dimerize,
            rate_function=rate_function_reverse_dimerize
        )
        system.topologies.add_structural_reaction(
            "Fail_Reverse_Dimerize",
            topology_type="Actin-Dimer#Fail",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_trimerize_reaction(system, rate, reaction_distance):
        '''
        attach a monomer to a dimer
        '''
        for i in range(1,4):
            system.topologies.add_spatial_reaction(
                "Trimerize{}: Actin-Dimer(actin#barbed_ATP_{}) + {}".format(i, i,
                "Actin-Monomer(actin#free_ATP) -> {}".format(
                "Actin-Trimer#Growing(actin#ATP_{}--actin#new_ATP)".format(i))),
                rate=rate, radius=reaction_distance
            )

        system.topologies.add_structural_reaction(
            "Finish_Trimerize",
            topology_type="Actin-Trimer#Growing",
            reaction_function=reaction_function_finish_trimerize,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_trimerize_reverse_reaction(system):
        '''
        detach a monomer from a dimer
        '''
        system.topologies.add_structural_reaction(
            "Reverse_Trimerize",
            topology_type="Actin-Trimer",
            reaction_function=reaction_function_reverse_trimerize,
            rate_function=rate_function_reverse_trimerize
        )
        system.topologies.add_structural_reaction(
            "Fail_Reverse_Trimerize",
            topology_type="Actin-Trimer#Fail",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_nucleate_reaction(system, rate_ATP, rate_ADP, reaction_distance):
        '''
        reversibly attach a monomer to a trimer
        '''
        for i in range(1,4):
            system.topologies.add_spatial_reaction(
                "Barbed_Growth_Nucleate_ATP{}: {}".format(i,
                "Actin-Trimer(actin#barbed_ATP_{}) + {}".format(i,
                "Actin-Monomer(actin#free_ATP) -> {}".format(
                "Polymer#GrowingBarbed(actin#ATP_{}--actin#new_ATP)".format(i)))),
                rate=rate_ATP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Barbed_Growth_Nucleate_ADP{}: {}".format(i,
                "Actin-Trimer(actin#barbed_ATP_{}) + {}".format(i,
                "Actin-Monomer(actin#free) -> {}".format(
                "Polymer#GrowingBarbed(actin#ATP_{}--actin#new)".format(i)))),
                rate=rate_ADP, radius=reaction_distance
            )

    @staticmethod
    def add_pointed_growth_reaction(system, rate_ATP, rate_ADP, reaction_distance):
        '''
        attach a monomer to the pointed end of a filament
        '''
        for i in range(1,4):
            system.topologies.add_spatial_reaction(
                "Pointed_Growth_ATP1{}: Polymer(actin#pointed_{}) + {}".format(i, i,
                "Actin-Monomer(actin#free_ATP) -> {}".format(
                "Polymer#GrowingPointed(actin#{}--actin#new_ATP)".format(i))),
                rate=rate_ATP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Pointed_Growth_ATP2{}: Polymer(actin#pointed_ATP_{}) + {}".format(i, i,
                "Actin-Monomer(actin#free_ATP) -> {}".format(
                "Polymer#GrowingPointed(actin#ATP_{}--actin#new_ATP)".format(i))),
                rate=rate_ATP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Pointed_Growth_ADP1{}: Polymer(actin#pointed_{}) + {}".format(i, i,
                "Actin-Monomer(actin#free) -> {}".format(
                "Polymer#GrowingPointed(actin#{}--actin#new)".format(i))),
                rate=rate_ADP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Pointed_Growth_ADP2{}: Polymer(actin#pointed_ATP_{}) + {}".format(i, i,
                "Actin-Monomer(actin#free) -> {}".format(
                "Polymer#GrowingPointed(actin#ATP_{}--actin#new)".format(i))),
                rate=rate_ADP, radius=reaction_distance
            )

        system.topologies.add_structural_reaction(
            "Finish_Pointed_Growth",
            topology_type="Polymer#GrowingPointed",
            reaction_function=reaction_function_finish_pointed_grow,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_pointed_shrink_reaction(system):
        '''
        remove a monomer from the pointed end of a filament
        '''
        system.topologies.add_structural_reaction(
            "Pointed_Shrink_ATP",
            topology_type="Polymer",
            reaction_function=reaction_function_pointed_shrink_ATP,
            rate_function=rate_function_pointed_shrink_ATP
        )
        system.topologies.add_structural_reaction(
            "Pointed_Shrink_ADP",
            topology_type="Polymer",
            reaction_function=reaction_function_pointed_shrink_ADP,
            rate_function=rate_function_pointed_shrink_ADP
        )
        system.topologies.add_structural_reaction(
            "Fail_Pointed_Shrink_ATP",
            topology_type="Polymer#Fail-Pointed-Shrink-ATP",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )
        system.topologies.add_structural_reaction(
            "Fail_Pointed_Shrink_ADP",
            topology_type="Polymer#Fail-Pointed-Shrink-ADP",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )
        system.topologies.add_structural_reaction(
            "Cleanup_Shrink",
            topology_type="Polymer#Shrinking",
            reaction_function=reaction_function_cleanup_shrink,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_barbed_growth_reaction(system, rate_ATP, rate_ADP, reaction_distance):
        '''
        attach a monomer to the barbed end of a filament
        '''
        for i in range(1,4):
            system.topologies.add_spatial_reaction(
                "Barbed_Growth_ATP1{}: Polymer(actin#barbed_{}) + {}".format(i, i,
                "Actin-Monomer(actin#free_ATP) -> {}".format(
                "Polymer#GrowingBarbed(actin#{}--actin#new_ATP)".format(i))),
                rate=rate_ATP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Barbed_Growth_ATP2{}: Polymer(actin#barbed_ATP_{}) + {}".format(i, i,
                "Actin-Monomer(actin#free_ATP) -> {}".format(
                "Polymer#GrowingBarbed(actin#ATP_{}--actin#new_ATP)".format(i))),
                rate=rate_ATP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Barbed_Growth_ADP1{}: Polymer(actin#barbed_{}) + {}".format(i, i,
                "Actin-Monomer(actin#free) -> {}".format(
                "Polymer#GrowingBarbed(actin#{}--actin#new)".format(i))),
                rate=rate_ADP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Barbed_Growth_ADP2{}: Polymer(actin#barbed_ATP_{}) + {}".format(i, i,
                "Actin-Monomer(actin#free) -> {}".format(
                "Polymer#GrowingBarbed(actin#ATP_{}--actin#new)".format(i))),
                rate=rate_ADP, radius=reaction_distance
            )
        system.topologies.add_spatial_reaction(
            "Branch_Barbed_Growth_ATP1: Polymer(actin#branch_barbed_1) + {}".format(
            "Actin-Monomer(actin#free_ATP) -> {}".format(
            "Polymer#GrowingBarbed(actin#branch_1--actin#new_ATP)")),
            rate=rate_ATP, radius=reaction_distance
        )
        system.topologies.add_spatial_reaction(
            "Branch_Barbed_Growth_ATP2: Polymer(actin#branch_barbed_ATP_1) + {}".format(
            "Actin-Monomer(actin#free_ATP) -> {}".format(
            "Polymer#GrowingBarbed(actin#branch_ATP_1--actin#new_ATP)")),
            rate=rate_ATP, radius=reaction_distance
        )
        system.topologies.add_spatial_reaction(
            "Branch_Barbed_Growth_ADP1: Polymer(actin#branch_barbed_1) + {}".format(
            "Actin-Monomer(actin#free) -> {}".format(
            "Polymer#GrowingBarbed(actin#branch_1--actin#new)")),
            rate=rate_ADP, radius=reaction_distance
        )
        system.topologies.add_spatial_reaction(
            "Branch_Barbed_Growth_ADP2: Polymer(actin#branch_barbed_ATP_1) + {}".format(
            "Actin-Monomer(actin#free) -> {}".format(
            "Polymer#GrowingBarbed(actin#branch_ATP_1--actin#new)")),
            rate=rate_ADP, radius=reaction_distance
        )

        system.topologies.add_structural_reaction(
            "Finish_Barbed_growth",
            topology_type="Polymer#GrowingBarbed",
            reaction_function=reaction_function_finish_barbed_grow,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_barbed_shrink_reaction(system):
        '''
        remove a monomer from the barbed end of a filament
        '''
        system.topologies.add_structural_reaction(
            "Barbed_Shrink_ATP",
            topology_type="Polymer",
            reaction_function=reaction_function_barbed_shrink_ATP,
            rate_function=rate_function_barbed_shrink_ATP
        )
        system.topologies.add_structural_reaction(
            "Barbed_Shrink_ADP",
            topology_type="Polymer",
            reaction_function=reaction_function_barbed_shrink_ADP,
            rate_function=rate_function_barbed_shrink_ADP
        )
        system.topologies.add_structural_reaction(
            "Fail_Barbed_Shrink_ATP",
            topology_type="Polymer#Fail-Barbed-Shrink-ATP",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )
        system.topologies.add_structural_reaction(
            "Fail_Barbed_Shrink_ADP",
            topology_type="Polymer#Fail-Barbed-Shrink-ADP",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_hydrolyze_reaction(system):
        '''
        hydrolyze ATP
        '''
        system.topologies.add_structural_reaction(
            "Hydrolysis_Actin",
            topology_type="Polymer",
            reaction_function=reaction_function_hydrolyze_actin,
            rate_function=rate_function_hydrolyze_actin
        )
        system.topologies.add_structural_reaction(
            "Fail_Hydrolysis_Actin",
            topology_type="Polymer#Fail-Hydrolysis-Actin",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )
        system.topologies.add_structural_reaction(
            "Hydrolysis_Arp",
            topology_type="Polymer",
            reaction_function=reaction_function_hydrolyze_arp,
            rate_function=rate_function_hydrolyze_arp
        )
        system.topologies.add_structural_reaction(
            "Fail_Hydrolysis_Arp",
            topology_type="Polymer#Fail-Hydrolysis-Arp",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_actin_nucleotide_exchange_reaction(system):
        '''
        exchange ATP for ADP in free actin monomers
        '''
        system.topologies.add_structural_reaction(
            "Nucleotide_Exchange_Actin",
            topology_type="Actin-Monomer",
            reaction_function=reaction_function_nucleotide_exchange_actin,
            rate_function=rate_function_nucleotide_exchange_actin
        )

    @staticmethod
    def add_arp23_nucleotide_exchange_reaction(system):
        '''
        exchange ATP for ADP in free Arp2/3 dimers
        '''
        system.topologies.add_structural_reaction(
            "Nucleotide_Exchange_Arp",
            topology_type="Arp23-Dimer",
            reaction_function=reaction_function_nucleotide_exchange_arp,
            rate_function=rate_function_nucleotide_exchange_arp
        )

    @staticmethod
    def add_arp23_bind_reaction(system, rate_ATP, rate_ADP, reaction_distance):
        '''
        add arp2/3 along filament to start a branch
        '''
        for i in range(1,4):
            system.topologies.add_spatial_reaction(
                "Arp_Bind_ATP1{}: Polymer(actin#ATP_{}) + Arp23-Dimer(arp2) -> {}".format(i, i,
                "Polymer#Branching(actin#ATP_{}--arp2#new)".format(i)),
                rate=rate_ATP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Arp_Bind_ATP2{}: Polymer(actin#ATP_{}) + Arp23-Dimer(arp2#ATP) -> {}".format(i, i,
                "Polymer#Branching(actin#ATP_{}--arp2#new_ATP)".format(i)),
                rate=rate_ATP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Arp_Bind_ADP1{}: Polymer(actin#{}) + Arp23-Dimer(arp2) -> {}".format(i, i,
                "Polymer#Branching(actin#{}--arp2#new)".format(i)),
                rate=rate_ADP, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Arp_Bind_ADP2{}: Polymer(actin#{}) + Arp23-Dimer(arp2#ATP) -> {}".format(i, i,
                "Polymer#Branching(actin#{}--arp2#new_ATP)".format(i)),
                rate=rate_ADP, radius=reaction_distance
            )

        system.topologies.add_structural_reaction(
            "Finish_Arp_Bind",
            topology_type="Polymer#Branching",
            reaction_function=reaction_function_finish_arp_bind,
            rate_function=rate_function_infinity
        )
        system.topologies.add_structural_reaction(
            "Cleanup_Fail_Arp_Bind_ATP",
            topology_type="Polymer#Fail-Branch-ATP",
            reaction_function=reaction_function_cleanup_shrink,
            rate_function=rate_function_infinity
        )
        system.topologies.add_structural_reaction(
            "Cleanup_Fail_Arp_Bind_ADP",
            topology_type="Polymer#Fail-Branch-ADP",
            reaction_function=reaction_function_cleanup_shrink,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_arp23_unbind_reaction(system):
        '''
        remove an arp2/3 that is not nucleated
        '''
        system.topologies.add_structural_reaction(
            "Arp_Unbind_ATP",
            topology_type="Polymer",
            reaction_function=reaction_function_arp23_unbind_ATP,
            rate_function=rate_function_arp23_unbind_ATP
        )
        system.topologies.add_structural_reaction(
            "Arp_Unbind_ADP",
            topology_type="Polymer",
            reaction_function=reaction_function_arp23_unbind_ADP,
            rate_function=rate_function_arp23_unbind_ADP
        )
        system.topologies.add_structural_reaction(
            "Fail_Arp_Unbind_ATP",
            topology_type="Polymer#Fail-Arp-Unbind-ATP",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )
        system.topologies.add_structural_reaction(
            "Fail_Arp_Unbind_ADP",
            topology_type="Polymer#Fail-Arp-Unbind-ADP",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_nucleate_branch_reaction(system, rate_ATP, rate_ADP, reaction_distance):
        '''
        add actin to arp2/3 to begin a branch
        '''
        system.topologies.add_spatial_reaction(
            "Barbed_Growth_Branch_ATP: Polymer(arp3) + Actin-Monomer(actin#free_ATP) -> {}".format(
            "Polymer#Branch-Nucleating(arp3#branched--actin#new_ATP)"),
            rate=rate_ATP, radius=reaction_distance
        )
        system.topologies.add_spatial_reaction(
            "Barbed_Growth_Branch_ADP: Polymer(arp3) + Actin-Monomer(actin#free) -> {}".format(
            "Polymer#Branch-Nucleating(arp3#branched--actin#new)"),
            rate=rate_ADP, radius=reaction_distance
        )
        system.topologies.add_structural_reaction(
            "Nucleate_Branch",
            topology_type="Polymer#Branch-Nucleating",
            reaction_function=reaction_function_finish_start_branch,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_debranch_reaction(system):
        '''
        remove a branch
        '''
        system.topologies.add_structural_reaction(
            "Debranch_ATP",
            topology_type="Polymer",
            reaction_function=reaction_function_debranching_ATP,
            rate_function=rate_function_debranching_ATP
        )
        system.topologies.add_structural_reaction(
            "Debranch_ADP",
            topology_type="Polymer",
            reaction_function=reaction_function_debranching_ADP,
            rate_function=rate_function_debranching_ADP
        )
        system.topologies.add_structural_reaction(
            "Fail_Debranch_ATP",
            topology_type="Polymer#Fail-Debranch-ATP",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )
        system.topologies.add_structural_reaction(
            "Fail_Debranch_ADP",
            topology_type="Polymer#Fail-Debranch-ADP",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_cap_bind_reaction(system, rate, reaction_distance):
        '''
        add capping protein to a barbed end to stop growth
        '''
        for i in range(1,4):
            system.topologies.add_spatial_reaction(
                "Cap_Bind1{}: Polymer(actin#barbed_{}) + Cap(cap) -> {}".format(i, i,
                "Polymer#Capping(actin#{}--cap#new)".format(i)),
                rate=rate, radius=reaction_distance
            )
            system.topologies.add_spatial_reaction(
                "Cap_Bind2{}: Polymer(actin#barbed_ATP_{}) + Cap(cap) -> {}".format(i, i,
                "Polymer#Capping(actin#ATP_{}--cap#new)".format(i)),
                rate=rate, radius=reaction_distance
            )

        system.topologies.add_structural_reaction(
            "Finish_Cap-Bind",
            topology_type="Polymer#Capping",
            reaction_function=reaction_function_finish_cap_bind,
            rate_function=rate_function_infinity
        )

    @staticmethod
    def add_cap_unbind_reaction(system):
        '''
        remove capping protein
        '''
        system.topologies.add_structural_reaction(
            "Cap_Unbind",
            topology_type="Polymer",
            reaction_function=reaction_function_cap_unbind,
            rate_function=rate_function_cap_unbind
        )
        system.topologies.add_structural_reaction(
            "Fail_Cap_Unbind",
            topology_type="Polymer#Fail-Cap-Unbind",
            reaction_function=reaction_function_reset_state,
            rate_function=rate_function_infinity
        )

# **************************** Callbacks ****************************************

box_size = 0.
def set_box_size(b):
    global box_size
    box_size = b
    return b

def reaction_function_reverse_dimerize(topology):
    '''
    reaction function for a dimer falling apart
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Reverse Dimerize")

    v_barbed = ReaddyUtil.get_first_vertex_of_types(topology,
        ["actin#barbed_ATP_1", "actin#barbed_ATP_2", "actin#barbed_ATP_3",
        "actin#barbed_1", "actin#barbed_2", "actin#barbed_3"])
    if v_barbed is None:
        if VERBOSE:
            print("failed to find barbed end")
        recipe.change_topology_type("Actin-Dimer#Fail")
        return recipe

    v_pointed = ReaddyUtil.get_first_neighbor(topology, v_barbed, [])
    if v_pointed is None:
        if VERBOSE:
            print("failed to find pointed end")
        recipe.change_topology_type("Actin-Dimer#Fail")
        return recipe

    recipe.remove_edge(v_barbed, v_pointed)
    recipe.change_particle_type(v_barbed, "actin#free_ATP")
    recipe.change_particle_type(v_pointed, "actin#free_ATP")
    recipe.change_topology_type("Actin-Monomer")

    return recipe

def reaction_function_finish_trimerize(topology):
    '''
    reaction function for a trimer forming
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Trimerize")

    v_new = ReaddyUtil.get_first_vertex_of_types(topology,
        ["actin#new", "actin#new_ATP"])
    if v_new is None:
        for v in topology.graph.get_vertices():
            print(topology.particle_type_of_vertex(v))
        raise GrowError('failed to find new actin vertex {}'.format(
            len(topology.graph.get_vertices())))

    v_neighbor1 = ReaddyUtil.get_first_neighbor(topology, v_new, [])
    if v_neighbor1 is None:
        for v in topology.graph.get_vertices():
            print(topology.particle_type_of_vertex(v))
        raise GrowError('failed to find neighbor of new actin vertex {}'.format(
            topology.particle_type_of_vertex(v_new)))

    v_neighbor2 = ReaddyUtil.get_first_neighbor(topology, v_neighbor1, [v_new])
    if v_neighbor2 is None:
        for v in topology.graph.get_vertices():
            print(topology.particle_type_of_vertex(v))
        raise GrowError(
            'failed to find non-new neighbor of actin vertex {}'.format(
            topology.particle_type_of_vertex(v_neighbor1)))

    ReaddyUtil.set_flags(topology, recipe, v_new,
        ["barbed", str(ActinUtil.get_actin_number(topology, v_neighbor1, 1))],
        ["new"], True)
    ActinUtil.set_new_trimer_vertex_position(
        topology, recipe, v_new, v_neighbor2, v_neighbor1)

    recipe.change_topology_type("Actin-Trimer")

    return recipe

def reaction_function_reverse_trimerize(topology):
    '''
    reaction function for removing ATP-actin from a trimer
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Reverse Trimerize")

    v_barbed = ReaddyUtil.get_first_vertex_of_types(topology,
        ["actin#barbed_ATP_1", "actin#barbed_ATP_2", "actin#barbed_ATP_3",
        "actin#barbed_1", "actin#barbed_2", "actin#barbed_3"])
    if v_barbed is None:
        if VERBOSE:
            print("failed to find barbed end")
        recipe.change_topology_type("Actin-Trimer#Fail")
        return recipe

    v_neighbor = ReaddyUtil.get_first_neighbor(topology, v_barbed, [])
    if v_neighbor is None:
        if VERBOSE:
            print("failed to find neighbor of barbed end")
        recipe.change_topology_type("Actin-Trimer#Fail")
        return recipe

    recipe.remove_edge(v_barbed, v_neighbor)
    recipe.change_particle_type(v_barbed, "actin#free_ATP")
    ReaddyUtil.set_flags(topology, recipe, v_neighbor, ["barbed"], [], True)
    recipe.change_topology_type("Polymer#Shrinking")

    return recipe

def reaction_function_finish_pointed_grow(topology):
    '''
    reaction function for the pointed end growing
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Grow Pointed")

    v_new = ReaddyUtil.get_first_vertex_of_types(topology,
        ["actin#new", "actin#new_ATP"])
    if v_new is None:
        for v in topology.graph.get_vertices():
            print(topology.particle_type_of_vertex(v))
        raise GrowError('failed to find new actin vertex {}'.format(
            len(topology.graph.get_vertices())))

    v_neighbor = ReaddyUtil.get_first_neighbor(topology, v_new, [])
    if v_neighbor is None:
        for v in topology.graph.get_vertices():
            print(topology.particle_type_of_vertex(v))
        raise GrowError('failed to find neighbor of new actin vertex {}'.format(
            topology.particle_type_of_vertex(v_new)))

    ReaddyUtil.set_flags(topology, recipe, v_new,
        ["pointed", str(ActinUtil.get_actin_number(topology, v_neighbor, -1))],
        ["new"], True)
    recipe.change_topology_type("Polymer")

    ActinUtil.set_end_vertex_position(topology, recipe, v_new, False)

    return recipe

def reaction_function_finish_barbed_grow(topology):
    '''
    reaction function for the barbed end growing
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Grow Barbed")

    v_new = ReaddyUtil.get_first_vertex_of_types(topology,
        ["actin#new", "actin#new_ATP"])
    if v_new is None:
        for v in topology.graph.get_vertices():
            print(topology.particle_type_of_vertex(v))
        raise GrowError('failed to find new actin vertex {}'.format(
            len(topology.graph.get_vertices())))

    v_neighbor = ReaddyUtil.get_first_neighbor(topology, v_new, [])
    if v_neighbor is None:
        for v in topology.graph.get_vertices():
            print(topology.particle_type_of_vertex(v))
        raise GrowError('failed to find neighbor of new actin vertex {}'.format(
            topology.particle_type_of_vertex(v_new)))

    print("attach to {}".format(topology.particle_type_of_vertex(v_neighbor)))

    ReaddyUtil.set_flags(topology, recipe, v_new,
        ["barbed", str(ActinUtil.get_actin_number(topology, v_neighbor, 1))],
        ["new"], True)
    ActinUtil.set_end_vertex_position(topology, recipe, v_new, True)
    recipe.change_topology_type("Polymer")

    return recipe

def reaction_function_finish_arp_bind(topology):
    '''
    reaction function to finish a branching reaction
    (triggered by a spatial reaction)
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Try Branch")

    v_arp2, v_arp3 = ActinUtil.get_new_arp23(topology)
    if v_arp2 is None or v_arp3 is None:
        raise BranchError('failed to find new arp2/3 vertices')

    v_actin_pointed = ReaddyUtil.get_first_neighbor(topology, v_arp2, [v_arp3])
    if v_actin_pointed is None:
        raise BranchError("failed to find actin vertex binding to arp2")

    # make sure arp3 binds to the barbed end neighbor of the actin bound to arp2
    N_barbed = ActinUtil.get_actin_number(topology, v_actin_pointed, 1)
    v_actin_barbed = ReaddyUtil.get_neighbor_of_types(
        topology, v_actin_pointed,
        ["actin#ATP_{}".format(N_barbed), "actin#{}".format(N_barbed)], [])
    if v_actin_barbed is None:
        if VERBOSE:
            print("failed to find actin vertex to bind to arp3")
        ActinUtil.cancel_failed_branch_reaction(
            topology, recipe, v_actin_pointed, v_arp2)
        return recipe

    ReaddyUtil.set_flags(topology, recipe, v_arp2, [], ["new"], True)
    recipe.add_edge(v_actin_barbed, v_arp3)
    recipe.change_topology_type("Polymer")

    ActinUtil.set_arp23_vertex_position(
        topology, recipe, v_arp2, v_arp3, v_actin_pointed, v_actin_barbed)

    return recipe

def reaction_function_finish_start_branch(topology):
    '''
    reaction function for adding the first actin to an arp2/3 to start a branch
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Start Branch")

    v_new = ReaddyUtil.get_first_vertex_of_types(topology,
        ["actin#new", "actin#new_ATP"])
    if v_new is None:
        for v in topology.graph.get_vertices():
            print(topology.particle_type_of_vertex(v))
        raise GrowError('failed to find new actin vertex {}'.format(
            len(topology.graph.get_vertices())))

    ReaddyUtil.set_flags(
        topology, recipe, v_new, ["barbed", "1", "branch"], ["new"], True)
    recipe.change_topology_type("Polymer")

    ActinUtil.set_end_vertex_position(topology, recipe, v_new, True)

    return recipe

def reaction_function_pointed_shrink_ATP(topology):
    '''
    reaction function to remove an ATP-actin from the pointed end
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Shrink pointed ATP")

    if not ActinUtil.do_shrink(topology, recipe, False, True):
        recipe.change_topology_type("Polymer#Fail-Pointed-Shrink-ATP")

    return recipe

def reaction_function_pointed_shrink_ADP(topology):
    '''
    reaction function to remove an ADP-actin from the pointed end
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Shrink pointed ADP")

    if not ActinUtil.do_shrink(topology, recipe, False, False):
        recipe.change_topology_type("Polymer#Fail-Pointed-Shrink-ADP")

    return recipe

def reaction_function_barbed_shrink_ATP(topology):
    '''
    reaction function to remove an ATP-actin from the barbed end
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Shrink barbed ATP")

    if not ActinUtil.do_shrink(topology, recipe, True, True):
        recipe.change_topology_type("Polymer#Fail-Barbed-Shrink-ATP")

    return recipe

def reaction_function_barbed_shrink_ADP(topology):
    '''
    reaction function to remove an ADP-actin from the barbed end
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Shrink barbed ADP")

    if not ActinUtil.do_shrink(topology, recipe, True, False):
        recipe.change_topology_type("Polymer#Fail-Barbed-Shrink-ADP")

    return recipe

def reaction_function_cleanup_shrink(topology):
    '''
    reaction function for finishing a reverse polymerization reaction
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Cleanup Shrink")

    if len(topology.graph.get_vertices()) < 2:

        v_cap = ReaddyUtil.get_vertex_of_type(topology, "cap", True)
        if v_cap is not None:
            if VERBOSE:
                    print("cleaned up Capping Protein")
            recipe.change_topology_type("Cap")
            return recipe

        if VERBOSE:
            print("cleaned up Actin Monomer")
        recipe.change_topology_type("Actin-Monomer")
        return recipe


    if len(topology.graph.get_vertices()) < 3:

        v_arp3 = ReaddyUtil.get_vertex_of_type(topology, "arp3", True)
        if v_arp3 is not None:
            if VERBOSE:
                print("cleaned up Arp2/3 dimer")
            recipe.change_topology_type("Arp23-Dimer")
            return recipe

        if VERBOSE:
            print("cleaned up Actin Dimer")
        recipe.change_topology_type("Actin-Dimer")
        return recipe

    if len(topology.graph.get_vertices()) < 4:
        if VERBOSE:
            print("cleaned up Actin Trimer")
        recipe.change_topology_type("Actin-Trimer")
        return recipe

    if VERBOSE:
        print("cleaned up Actin Polymer")
    recipe.change_topology_type("Polymer")
    return recipe

def reaction_function_hydrolyze_actin(topology):
    '''
    reaction function to hydrolyze a filamentous ATP-actin to ADP-actin
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Hydrolyze Actin")

    v = ReaddyUtil.get_random_vertex_of_types(
        topology, get_types_of_all_numbers("actin#ATP")
        + get_types_of_all_numbers("actin#pointed_ATP")
        + get_types_of_all_numbers("actin#barbed_ATP"))

    if v is None:
        if VERBOSE:
            print("failed to find ATP-actin")
        recipe.change_topology_type("Polymer#Fail-Hydrolysis-Actin")
        return recipe

    ReaddyUtil.set_flags(topology, recipe, v, [], ["ATP"], True)

    return recipe

def reaction_function_hydrolyze_arp(topology):
    '''
    reaction function to hydrolyze a arp2/3 (hydrolyze the ATP in arp2,
    but for simplicity of the model,
    nucleotide state is tied to arp3 since it connects to the branch)
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Hydrolyze Arp2/3")

    v = ReaddyUtil.get_random_vertex_of_types(topology, ["arp2#ATP"])
    if v is None:
        if VERBOSE:
            print("failed to find arp2/3 with ATP")
        recipe.change_topology_type("Polymer#Fail-Hydrolysis-Arp")
        return recipe

    ReaddyUtil.set_flags(topology, recipe, v, [], ["ATP"], True)

    return recipe

def reaction_function_nucleotide_exchange_actin(topology):
    '''
    reaction function to exchange ATP for ADP in free actin
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Nucleotide Exchange Actin")

    v = ReaddyUtil.get_vertex_of_type(topology, "actin#free", True)

    if v is None:
        if VERBOSE:
            print("failed to find ADP-actin")
        recipe.change_topology_type("Polymer#Fail-Nucleotide-Exchange-Actin")
        return recipe

    ReaddyUtil.set_flags(topology, recipe, v, ["ATP"], [], True)

    return recipe

def reaction_function_nucleotide_exchange_arp(topology):
    '''
    reaction function to exchange ATP for ADP in free Arp2/3
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Nucleotide Exchange Arp2/3")

    v = ReaddyUtil.get_vertex_of_type(topology, "arp2", True)

    if v is None:
        if VERBOSE:
            print("failed to find ADP-arp2")
        recipe.change_topology_type("Polymer#Fail-Nucleotide-Exchange-Arp")
        return recipe

    ReaddyUtil.set_flags(topology, recipe, v, ["ATP"], [], True)

    return recipe

def reaction_function_arp23_unbind_ATP(topology):
    '''
    reaction function to dissociate an arp2/3 with ATP from a mother filament
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Remove Arp2/3 ATP")

    ActinUtil.do_arp23_unbind(topology, recipe, True)

    return recipe

def reaction_function_arp23_unbind_ADP(topology):
    '''
    reaction function to dissociate an arp2/3 with ADP from a mother filament
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Remove Arp2/3 ADP")

    ActinUtil.do_arp23_unbind(topology, recipe, False)

    return recipe

def reaction_function_debranching_ATP(topology):
    '''
    reaction function to detach a branch filament from arp2/3 with ATP
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Debranching ATP")

    ActinUtil.do_debranching(topology, recipe, True)

    return recipe

def reaction_function_debranching_ADP(topology):
    '''
    reaction function to detach a branch filament from arp2/3 with ADP
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Debranching ADP")

    ActinUtil.do_debranching(topology, recipe, False)

    return recipe

def reaction_function_finish_cap_bind(topology):
    '''
    reaction function for adding a capping protein
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Finish Cap Bind")

    v_new = ReaddyUtil.get_first_vertex_of_types(topology, ["cap#new"])
    if v_new is None:
        if VERBOSE:
            print("failed to find new capping protein")
        return recipe

    ReaddyUtil.set_flags(topology, recipe, v_new, ["bound"], ["new"], True)
    recipe.change_topology_type("Polymer")

    return recipe

def reaction_function_cap_unbind(topology):
    '''
    reaction function to detach capping protein from a barbed end
    '''
    recipe = readdy.StructuralReactionRecipe(topology)
    if VERBOSE:
        print("Remove Cap")

    v_cap = ReaddyUtil.get_random_vertex_of_types(topology, ["cap#bound"])
    if v_cap is None:
        if VERBOSE:
            print("failed to find cap")
        recipe.change_topology_type("Polymer#Fail-Cap-Unbind")
        return recipe

    v_actin = ReaddyUtil.get_neighbor_of_types(
        topology, v_cap, get_types_of_all_numbers("actin")
        + get_types_of_all_numbers("actin#ATP"), [])
    if v_actin is None:
        if VERBOSE:
            print("failed to find actin bound to cap")
        recipe.change_topology_type("Polymer#Fail-Cap-Unbind")
        return recipe

    recipe.remove_edge(v_cap, v_actin)
    ReaddyUtil.set_flags(topology, recipe, v_cap, [], ["bound"], True)
    ReaddyUtil.set_flags(topology, recipe, v_actin, ["barbed"], [], True)
    recipe.change_topology_type("Polymer#Shrinking")

    return recipe

def rate_function_reverse_dimerize(topology):
    '''
    rate function for a dimer falling apart
    '''
    return parameters["dimerize_reverse_rate_rate"]

def rate_function_reverse_trimerize(topology):
    '''
    rate function for removing ATP-actin from a trimer
    '''
    return parameters["trimerize_reverse_rate"]

def rate_function_barbed_shrink_ATP(topology):
    '''
    rate function for removing ATP-actin from the barbed end
    '''
    return parameters["barbed_shrink_ATP_rate"]

def rate_function_barbed_shrink_ADP(topology):
    '''
    rate function for removing ADP-actin from the barbed end
    '''
    return parameters["barbed_shrink_ADP_rate"]

def rate_function_pointed_shrink_ATP(topology):
    '''
    rate function for removing ATP-actin from the pointed end
    '''
    return parameters["pointed_shrink_ATP_rate"]

def rate_function_pointed_shrink_ADP(topology):
    '''
    rate function for removing ADP-actin from the pointed end
    '''
    return parameters["pointed_shrink_ADP_rate"]

def rate_function_hydrolyze_actin(topology):
    '''
    rate function for hydrolyzing filamentous ATP-actin to ADP-actin
    '''
    return parameters["hydrolysis_actin_rate"]

def rate_function_hydrolyze_arp(topology):
    '''
    rate function for hydrolyzing bound arp2/3
    '''
    return parameters["hydrolysis_arp_rate"]

def rate_function_nucleotide_exchange_actin(topology):
    '''
    rate function for exchanging an ATP for ADP in free actin
    '''
    return parameters["nucleotide_exchange_actin_rate"]

def rate_function_nucleotide_exchange_arp(topology):
    '''
    rate function for exchanging an ATP for ADP in free Arp2/3
    '''
    return parameters["nucleotide_exchange_arp_rate"]

def rate_function_debranching_ATP(topology):
    '''
    rate function for dissociation of a daughter filament from an arp2/3 with ATP
    '''
    return parameters["debranching_ATP_rate"]

def rate_function_debranching_ADP(topology):
    '''
    rate function for dissociation of a daughter filament from an arp2/3 with ADP
    '''
    return parameters["debranching_ADP_rate"]

def rate_function_arp23_unbind_ATP(topology):
    '''
    rate function for dissociation of bound arp2/3 with ATP
    '''
    return parameters["arp_unbind_ATP_rate"]

def rate_function_arp23_unbind_ADP(topology):
    '''
    rate function for dissociation of bound arp2/3 with ADP
    '''
    return parameters["arp_unbind_ADP_rate"]

def rate_function_cap_unbind(topology):
    '''
    rate function for dissociation of capping protein
    '''
    return parameters["cap_unbind_rate"]

class Error(Exception):
    pass

class GrowError(Error):
    pass

class BranchError(Error):
    pass
