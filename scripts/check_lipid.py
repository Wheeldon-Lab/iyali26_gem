import cobra
m = cobra.io.read_sbml_model('model.xml')

ceramide_names = set()
for met in m.metabolites:
    if 'ceramide' in (met.name or '').lower():
        # 用 "_" 之前的部分作为基础名
        base = (met.name or '').rsplit('_', 1)[0].strip()
        ceramide_names.add(base)

for n in sorted(ceramide_names):
    print(n)