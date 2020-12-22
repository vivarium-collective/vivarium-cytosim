#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

VERBOSE = False

# from PDB structure
initial_positions = [
    np.array([-0.1127203,  -2.675011,  -14.24995]), # actin6 (pointed)
    np.array([ 2.033443,   -0.3248997, -17.09995]), # actin7
    np.array([-0.5967388,  -2.116877,  -19.94994]), # actin8 (actin_arp2)
    np.array([ 2.375738,   -0.9795914, -22.79993]), # actin9 (actin_arp3)
    np.array([-0.7788677,  -1.400904,  -25.64992]), # actin10 (barbed)
    np.array([-2.600235,   -5.706778,  -19.27911]), # arp2
    np.array([ 0.1346893,  -8.574478,  -20.60971]), # arp3
    np.array([-2.736113,  -11.19299,   -19.05002]), # actin13 (first branch actin)
    np.array([ 0.3656578, -14.09693,   -19.49452]), # actin14 (second branch actin)
    np.array([-2.81676,   -16.93661,   -19.24974]), # actin15 (third branch actin)
]

VECTOR_TO_NEW_POINTED = initial_positions[0] - initial_positions[2]
VECTOR_TO_NEW_BARBED = initial_positions[4] - initial_positions[2]
VECTOR_TO_NEW_ARP2 = initial_positions[5] - initial_positions[2]
VECTOR_TO_NEW_ARP3 = initial_positions[6] - initial_positions[2]
VECTOR_TO_NEW_BRANCH_ACTIN = [
    initial_positions[7] - initial_positions[2],
    initial_positions[8] - initial_positions[2],
    initial_positions[9] - initial_positions[2]
]
INITIAL_ORIENTATION = get_actin_orientation(
    [initial_positions[1], initial_positions[2], initial_positions[3]])
ANGLE_BETWEEN_ACTINS = 1.475664 #radians
