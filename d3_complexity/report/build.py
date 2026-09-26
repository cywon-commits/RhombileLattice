"""Assemble d3_explainer.html from template.html and figs.json (run gen_svgs.py first)."""
import json, re
F = json.load(open("figs.json"))
t = open("template.html").read()
missing = set(re.findall(r"\{\{(\w+)\}\}", t)) - set(F)
assert not missing, missing
t = re.sub(r"\{\{(\w+)\}\}", lambda m: F[m.group(1)], t)
open("d3_explainer.html", "w").write(t)
print("wrote d3_explainer.html", len(t))
