def backbone_hydrophobic(structure, non_bonded, nb_cutoff=0.5):
  # Backbone springs
  distances = md.compute_distances(structure, atom_pairs=non_bonded).ravel()
  filtered =  np.where(distances < nb_cutoff)[0]
  output = list(distances[filtered])
  # Hydrophobic springs
  atom_resid = [str(atom.residue)[:3] for atom in structure.top.atoms]
  hydrophobic = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'MET', 'TRP']
  for pair in non_bonded:
    if np.all(np.array([atom_resid[pair[0]] in hydrophobic, atom_resid[pair[1]] in hydrophobic])) and if pair not in output:
      output.append(pair)
  return set(output)
