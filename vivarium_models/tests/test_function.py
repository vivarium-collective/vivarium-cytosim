#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tests for actin ReaDDy models
"""

from vivarium_models import ReaddyActinProcess


def test_readdy_actin_process():
    """
    Test the initial ReaDDy actin process
    """
    output = ReaddyActinProcess.run_readdy_actin_process()

    found_monomer = False
    found_dimer = False
    assert len(output["topologies"]) == 2
    for t in output["topologies"]:
        top = output["topologies"][t]
        if top["type"][1] == "Actin-Monomer":
            found_monomer = True
            assert len(top["particles"]) == 1
            assert "actin#free" in top["particles"]["0"]["type"][1]
            assert len(top["particles"]["0"]["neighbors"][1]) == 0
        if top["type"][1] == "Arp23-Dimer":
            found_dimer = True
            assert len(top["particles"]) == 2
            assert "arp2" in top["particles"]["1"]["type"][1]
            assert top["particles"]["1"]["neighbors"][1] == [2]
            assert "arp3" in top["particles"]["2"]["type"][1]
            assert top["particles"]["2"]["neighbors"][1] == [1]
    assert found_monomer
    assert found_dimer
