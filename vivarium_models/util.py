def create_monomer_update(previous_monomers, new_monomers):
    topology_update = {}
    particle_update = {}
    if len(previous_monomers["topologies"]) > 0:
        topology_update["_delete"] = list(previous_monomers["topologies"].keys())
    if len(previous_monomers["particles"]) > 0:
        particle_update["_delete"] = list(previous_monomers["particles"].keys())

    add_topologies = []
    add_particles = []
    for fiber_id in new_monomers["topologies"]:
        particle_ids = new_monomers["topologies"][fiber_id]["particle_ids"]
        add_topologies.append(
            {
                "key": fiber_id,
                "state": {
                    "type_name": new_monomers["topologies"][fiber_id]["type_name"],
                    "particle_ids": particle_ids,
                },
            }
        )
        for monomer_id in particle_ids:
            monomer = new_monomers["particles"][monomer_id]
            add_particles.append(
                {
                    "key": monomer_id,
                    "state": {
                        "type_name": monomer["type_name"],
                        "position": monomer["position"],
                        "neighbor_ids": monomer["neighbor_ids"],
                    },
                }
            )
    topology_update["_add"] = add_topologies
    particle_update["_add"] = add_particles

    return {"monomers": {"topologies": topology_update, "particles": particle_update}}
