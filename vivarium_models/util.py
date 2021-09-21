def agents_update(existing, projected):
    update = {
        '_add': [],
        '_delete': []}

    for id, state in projected.items():
        if id in existing:
            update[id] = state
        else:
            update['_add'].append({
                "key": id,
                "state": state})

    for existing_id in existing.keys():
        if existing_id not in projected:
            update['_delete'].append(existing_id)

    return update


def create_monomer_update(previous_monomers, new_monomers):
    topologies_update = agents_update(
        previous_monomers['topologies'],
        new_monomers['topologies'])

    particles_update = agents_update(
        previous_monomers['particles'],
        new_monomers['particles'])

    return {
        "monomers": {
            "topologies": topologies_update,
            "particles": particles_update}}
