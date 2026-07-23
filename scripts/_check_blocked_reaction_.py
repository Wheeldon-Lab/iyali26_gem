import csv
from cobra.io import read_sbml_model
from cobra.flux_analysis import find_blocked_reactions

m = read_sbml_model('/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/model.xml')
m.solver = 'glpk'
blocked = set(find_blocked_reactions(m, processes=1))

out_path = '/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/data/blocked_reactions.csv'
with open(out_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['id', 'name', 'reaction', 'gpr', 'subsystem', 'compartments', 'has_gpr'])
    for r in m.reactions:
        if r.id not in blocked:
            continue
        comps = ','.join(sorted({m.compartment for m in r.metabolites}))
        writer.writerow([
            r.id,
            r.name,
            r.reaction,
            r.gene_reaction_rule,
            r.subsystem,
            comps,
            bool(r.gene_reaction_rule),
        ])

print(f'Wrote {len(blocked)} blocked reactions to {out_path}')
