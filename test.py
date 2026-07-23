from cobra.io import read_sbml_model

model = read_sbml_model("model.xml")

def is_unannotated(rxn):
    ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
    return not {k for k in ann if k != "sbo"}

unannotated = [rxn for rxn in model.reactions if is_unannotated(rxn)]

print(f"total {len(unannotated)} unannotated reacions（ {len(model.reactions)} ）")
for rxn in unannotated:
    print(f"  {rxn.id:15s}  {rxn.name}")
