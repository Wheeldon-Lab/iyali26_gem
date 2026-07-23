import csv
from pathlib import Path
from cobra.io import read_sbml_model
from cobra.flux_analysis import find_blocked_reactions

if __name__ == '__main__':
    m = read_sbml_model('model.xml')
    m.solver = 'glpk'

    data_dir = Path(__file__).parent.parent.parent / 'data'

    # 当前 medium
    print('=== Current medium (open uptakes) ===')
    medium_rows = []
    for rxn_id, bound in m.medium.items():
        rxn = m.reactions.get_by_id(rxn_id)
        met = list(rxn.metabolites.keys())[0]
        print(f'  {rxn_id:20s}  bound={bound:8.1f}  met={met.name} ({met.id})')
        medium_rows.append({'reaction_id': rxn_id, 'bound': bound, 'met_name': met.name, 'met_id': met.id})

    with open(data_dir / 'minimum_medium.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['reaction_id', 'bound', 'met_name', 'met_id'])
        writer.writeheader()
        writer.writerows(medium_rows)

    # Blocked exchanges 的代谢物名称
    blocked = set(find_blocked_reactions(m))
    print(f'\n=== Blocked exchange metabolites (all 97) ===')
    blocked_ex = [(e, list(e.metabolites.keys())[0]) for e in m.exchanges if e.id in blocked]
    blocked_rows = []
    for e, met in sorted(blocked_ex, key=lambda x: x[1].name or ''):
        print(f'  {met.name:50s}  ({met.id})  exchange={e.id}')
        blocked_rows.append({'met_name': met.name, 'met_id': met.id, 'exchange_id': e.id})

    with open(data_dir / 'blocked_exchanges.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['met_name', 'met_id', 'exchange_id'])
        writer.writeheader()
        writer.writerows(blocked_rows)