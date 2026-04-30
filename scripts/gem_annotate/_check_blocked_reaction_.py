from cobra.io import read_sbml_model
from cobra.flux_analysis import find_blocked_reactions
from collections import Counter
m = read_sbml_model('model.xml')
m.solver = 'glpk'
blocked = set(find_blocked_reactions(m, processes=1))
blocked_with_gpr = [r for r in m.reactions if r.id in blocked and r.gene_reaction_rule]
blocked_no_gpr = [r for r in m.reactions if r.id in blocked and not r.gene_reaction_rule]
print(f'Blocked with GPR: {len(blocked_with_gpr)}')
print(f'Blocked without GPR: {len(blocked_no_gpr)}')
comps = Counter()
for r in m.reactions:
    if r.id in blocked:
        c = set()
        for met in r.metabolites:
            c.add(met.compartment)
        comps[frozenset(c)] += 1
print(f'\nBlocked by compartment combination:')
for c, n in comps.most_common(10):
    print(f'  {set(c)}: {n}')
full_deadend = 0
partial_deadend = 0
for r in m.reactions:
    if r.id in blocked:
        mets = list(r.metabolites.keys())
        all_blocked = all(
            all(rxn.id in blocked for rxn in met.reactions)
            for met in mets
        )
        if all_blocked:
            full_deadend += 1
        else:
            partial_deadend += 1
print(f'\nFully isolated (all connected rxns blocked): {full_deadend}')
print(f'Partially connected (some rxns carry flux): {partial_deadend}')
