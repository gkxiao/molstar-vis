"""Generate data.js: fully offline, inline structure data for the Mol* viewer."""
import json
import math
from collections import OrderedDict, defaultdict
from pathlib import Path

SRC = Path(r"C:/Users/gkxiao/Documents/blog/AIR/comfa_study/molstar")
DST = SRC / "molstar_viewer" / "data.js"

pdb_text = (SRC / "4ivc_receptor.pdb").read_text(encoding="utf-8", errors="replace")
sdf_text = (SRC / "test.sdf").read_text(encoding="utf-8", errors="replace")

# --- protein atoms (heavy only) for pocket detection -------------------------
protein_atoms = []
for line in pdb_text.splitlines():
    if line.startswith(("ATOM  ", "HETATM")):
        elem = line[76:78].strip() or line[12:16].strip()[0]
        if elem.upper() == "H":
            continue
        protein_atoms.append({
            "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            "chain": line[21].strip(),
            "seq": line[22:26].strip(),
            "name": line[17:20].strip(),
            "icode": line[26].strip(),
        })

# --- split SDF into individual pose records ---------------------------------
records = []
for chunk in sdf_text.split("$$$$"):
    lines = chunk.splitlines()
    # remove leading blank lines created by the split
    while lines and not lines[0].strip():
        lines.pop(0)
    if len(lines) < 4 or "V2000" not in lines[3]:
        continue
    title = lines[0].strip() or "unnamed"
    natoms = int(lines[3][0:3])
    atoms = []
    for line in lines[4:4 + natoms]:
        atoms.append({
            "xyz": (float(line[0:10]), float(line[10:20]), float(line[20:30])),
            "elem": line[31:34].strip() or "C",
        })
    heavy = [a for a in atoms if a["elem"].upper() != "H"]
    centroid = [round(sum(a["xyz"][i] for a in heavy) / max(1, len(heavy)), 3) for i in range(3)]
    formula = defaultdict(int)
    for a in atoms:
        formula[a["elem"].upper()] += 1
    formula_str = "".join(
        (k if v == 1 else f"{k}{v}") for k, v in sorted(formula.items(), key=lambda kv: ("C H".split().index(kv[0]) if kv[0] in ("C", "H") else 2, kv[0]))
    )

    # pocket residues within 5.5 A of any heavy ligand atom
    near = set()
    for pa in protein_atoms:
        for la in heavy:
            if abs(pa["xyz"][0] - la["xyz"][0]) > 5.5:
                continue
            if abs(pa["xyz"][1] - la["xyz"][1]) > 5.5:
                continue
            if abs(pa["xyz"][2] - la["xyz"][2]) > 5.5:
                continue
            if math.dist(pa["xyz"], la["xyz"]) <= 5.5:
                near.add((pa["chain"], pa["seq"], pa["icode"], pa["name"]))
                break
    pocket = []
    pocket_items = []
    for chain, seq, icode, name in sorted(near, key=lambda r: (r[0], int(r[1]) if r[1].isdigit() else 0)):
        pocket.append(f"{name} {seq}{icode} ({chain})" if icode else f"{name} {seq} ({chain})")
        item = {"auth_asym_id": chain, "auth_seq_id": int(seq)}
        if icode:
            item["pdbx_PDB_ins_code"] = icode
        pocket_items.append(item)

    # keep the original lines verbatim + the SDF record terminator
    sdf_record = "\n".join(lines).rstrip() + "\n$$$$\n"
    records.append({
        "name": title,
        "nAtoms": natoms,
        "nHeavy": len(heavy),
        "formula": formula_str,
        "centroid": centroid,
        "pocket": pocket,
        "pocketItems": pocket_items,
        "sdf": sdf_record,
    })

# --- group poses by molecule name -------------------------------------------
groups = OrderedDict()
for i, r in enumerate(records):
    group = groups.setdefault(r["name"], [])
    r["poseNo"] = len(group) + 1
    group.append(i)

data = {
    "protein": {
        "name": "JAK1 kinase (JH1 domain) - 4IVC receptor",
        "pdb": pdb_text,
        "nProteinAtoms": len(protein_atoms),
    },
    "poses": records,
    "groups": [{"name": k, "poseIds": v} for k, v in groups.items()],
}

DST.write_text(
    "// Auto-generated structure data (offline, no fetch required).\n"
    "// Source: 4ivc_receptor.pdb + test.sdf\n"
    "window.VIEWER_DATA = " + json.dumps(data, ensure_ascii=False, separators=(",", ":")) + ";\n",
    encoding="utf-8",
)
print(f"poses: {len(records)}, groups: {len(groups)}, data.js bytes: {DST.stat().st_size}")
print("largest pose group:", max((len(v), k) for k, v in groups.items()))
print("pocket residues in first pose:", len(records[0]["pocket"]), records[0]["pocket"][:8])
