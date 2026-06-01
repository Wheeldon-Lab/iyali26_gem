import csv
from cobra.io import read_sbml_model
from cobra.flux_analysis import find_blocked_reactions

m = read_sbml_model('/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/model.xml')
m.solver = 'glpk'

blocked_rxn_ids = set(find_blocked_reactions(m, processes=1))

dbs = ['bigg.metabolite', 'kegg.compound', 'metanetx.chemical', 'chebi', 'hmdb', 'inchikey', 'pubchem.compound', 'seed.compound']

out_path = '/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/data/lacking_metabolites.csv'
with open(out_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([
        'id', 'name', 'formula', 'charge', 'compartment',
        'has_annotation', 'missing_dbs',
        'missing_formula', 'missing_charge', 'missing_name',
        'n_reactions', 'n_blocked_reactions', 'all_reactions_blocked',
        'reaction_ids',
    ])
    for met in m.metabolites:
        ann = met.annotation if isinstance(met.annotation, dict) else {}
        has_annotation = bool(ann)
        missing_dbs = ';'.join(db for db in dbs if db not in ann)
        missing_formula = not (met.formula and met.formula.strip())
        missing_charge = met.charge is None
        missing_name = not (met.name and met.name.strip())

        rxn_ids = [r.id for r in met.reactions]
        n_blocked = sum(1 for rid in rxn_ids if rid in blocked_rxn_ids)
        all_blocked = len(rxn_ids) > 0 and n_blocked == len(rxn_ids)

        # only write rows that have at least one missing field or a blocked reaction
        if not (missing_formula or missing_charge or missing_name or not has_annotation or n_blocked > 0):
            continue

        writer.writerow([
            met.id,
            met.name,
            met.formula,
            met.charge,
            met.compartment,
            has_annotation,
            missing_dbs,
            missing_formula,
            missing_charge,
            missing_name,
            len(rxn_ids),
            n_blocked,
            all_blocked,
            ';'.join(rxn_ids),
        ])

print(f'Wrote to {out_path}')
