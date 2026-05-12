from cobra.io import read_sbml_model
from collections import Counter
m = read_sbml_model('model.xml')
m.solver = 'glpk'
total = len(m.metabolites)
has_ann = sum(1 for met in m.metabolites if isinstance(met.annotation, dict) and met.annotation)
has_formula = sum(1 for met in m.metabolites if met.formula and met.formula.strip())
has_charge = sum(1 for met in m.metabolites if met.charge is not None and met.charge != 0)
has_name = sum(1 for met in m.metabolites if met.name and met.name.strip())

# Per-database coverage
dbs = ['bigg.metabolite', 'kegg.compound', 'metanetx.chemical', 'chebi', 'hmdb', 'inchikey', 'pubchem.compound', 'seed.compound']
print(f'=== Metabolite overview ===')
print(f'Total: {total}')
print(f'Has any annotation: {has_ann} ({100*has_ann/total:.1f}%)')
print(f'Has formula: {has_formula} ({100*has_formula/total:.1f}%)')
print(f'Has charge: {has_charge} ({100*has_charge/total:.1f}%)')
print(f'Has name: {has_name} ({100*has_name/total:.1f}%)')

print(f'\n=== Per-database coverage ===')
for db in dbs:
    n = sum(1 for met in m.metabolites if isinstance(met.annotation, dict) and db in met.annotation)
    print(f'  {db:25s}: {n}/{total} ({100*n/total:.1f}%)')

# Unannotated metabolites
unann = [met for met in m.metabolites if not isinstance(met.annotation, dict) or not met.annotation]
print(f'\n=== Unannotated metabolites: {len(unann)} ===')
has_name_unann = sum(1 for met in unann if met.name and met.name.strip())
has_formula_unann = sum(1 for met in unann if met.formula and met.formula.strip())
print(f'  With name: {has_name_unann}')
print(f'  With formula: {has_formula_unann}')
print(f'  Samples:')
for met in unann[:15]:
    print(f'    {met.id:25s}  name={str(met.name)[:45]:45s}  formula={met.formula}  compartment={met.compartment}')

# Compartment distribution of unannotated
comp = Counter(met.compartment for met in unann)
print(f'\n  By compartment:')
for c, n in comp.most_common():
    print(f'    {c}: {n}')