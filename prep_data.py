#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prep_data.py — 为离线 Mol* 查看器生成 data.js（结构数据全部内嵌，无网络请求）

最简用法：把本脚本、受体 PDB、配体 SDF 放同一个目录，然后

    python prep_data.py

脚本自动找到唯一的 *.pdb 与 *.sdf，结果写到 molstar_viewer/data.js。
也可显式指定：

    python prep_data.py --pdb 7s8l.pdb --sdf linker-design.sdf
    python prep_data.py --pdb receptor.pdb --sdf poses.sdf \
                        --out molstar_viewer/data.js \
                        --name "MRGPRX2 / cortistatin-14 (7S8L)" \
                        --cutoff 6.0

参数
    --pdb / --sdf / --out   输入与输出（相对路径基于本脚本所在目录）
    --name                  蛋白显示名；不写则从 PDB 的 TITLE 行自动提取
    --cutoff                口袋残基距离阈值，Å（默认 5.5，只按重原子）
    --score-field NAME      指定哪个 SDF 属性作为列表里显示的打分
    --no-props              不提取 SDF 属性（打分等），data.js 更小
    --group-regex REGEX     含 1 个捕获组的正则，从分子名生成分组键；
                            例 '^(\\w+?\\d+)' 可把 Linker16 / Linker16a 归为一组
    --no-obabel             禁用 Open Babel，只用内置 V3000→V2000 转换器
    --strict                无法解析的 SDF 记录直接报错（默认跳过并警告）

────────────────────────────────────────────────────────────────────────
关于 V3000 SDF（重要）
    Mol* 5.13.1 的 SDF 解析器明确「不支持 V3000 的形式电荷」，所有原子电荷
    一律置 0（lib/mol-io/reader/sdf/parser.js）。这会让离子/电荷类相互作用
    判据失效。因此必须先把 V3000 转成 V2000：

        obabel -isdf linker-design.sdf -osdf -O linker-design-v2000.sdf

    Open Babel 写的 V2000 会生成标准 `M  CHG` 行，Mol* 能正确读取
    （lib/mol-model-formats/structure/mol.js）。本脚本检测到 V3000 时会
    自动调用本机 obabel 做这一步（生成 <原名>-v2000.sdf）；若未安装
    Open Babel，则使用内置的纯 Python 转换器兜底（同样写出 M  CHG）。

    已验证：Open Babel 3.1.1 转换 31 个分子后，坐标零偏差、元素/键数一致、
    CHG=+1/-1 → `M  CHG  2  31 1 34 -1`。

    V3000 也可以手动先转好，再把这个 V2000 文件传给 --sdf。
────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import shutil
import subprocess
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent

SCORE_CANDIDATES = [
    "r_i_docking_score", "r_i_glide_gscore", "r_i_glide_emodel",
    "r_psp_IFDScore", "r_i_Potential_Energy-S-OPLS", "r_psp_Prime_Energy",
    "docking_score", "glide_gscore", "Score", "score", "Affinity",
    # Cresset Spark / linker-design 常见字段
    "BIF%", "FScore+P", "FScore", "SScore+P", "SScore", "Sim", "2D Sim",
    "FSim", "SSim", "Field Penalty", "Excl Vol Pen",
]
# 不作为打分自动选取的字段
SCORE_BLACKLIST = re.compile(
    r"(formal charge|cluster|frequency|entry|subgroup|resolution|expdta|\bNC\b|prepared|pdb_)", re.I)
PROP_SKIP = re.compile(r"(_cresset|fieldpoints|fieldconstraints|SMILES|entry_id|subgroup)", re.I)
# V2000 原子行电荷码：1=+3 2=+2 3=+1 5=-1 6=-2 7=-3
CHARGE_CODE = {3: 1, 2: 2, 1: 3, -1: 5, -2: 6, -3: 7}


# ---------------------------------------------------------------- 参数
def parse_args():
    ap = argparse.ArgumentParser(
        description="为离线 Mol* 查看器生成 data.js",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例: python prep_data.py --pdb 7s8l.pdb --sdf linker-design.sdf",
    )
    ap.add_argument("--pdb", help="受体 PDB（默认自动探测）")
    ap.add_argument("--sdf", help="配体 SDF（默认自动探测）")
    ap.add_argument("--out", help="输出 data.js（默认 molstar_viewer/data.js）")
    ap.add_argument("--name", help="蛋白显示名（默认从 PDB TITLE 提取）")
    ap.add_argument("--cutoff", type=float, default=5.5, help="口袋残基阈值 Å（默认 5.5）")
    ap.add_argument("--score-field", help="强制指定打分属性名")
    ap.add_argument("--no-props", action="store_true", help="不提取 SDF 属性")
    ap.add_argument("--group-regex", help="用正则捕获组从分子名生成分组键")
    ap.add_argument("--no-obabel", action="store_true", help="只用内置 V3000 转换器")
    ap.add_argument("--strict", action="store_true", help="解析失败直接报错")
    return ap.parse_args()


def _abs(p):
    p = Path(p)
    return p if p.is_absolute() else BASE / p


def resolve(path_arg, patterns, label):
    if path_arg:
        p = _abs(path_arg)
        if not p.is_file():
            sys.exit(f"[错误] {label} 不存在：{p}")
        return p
    hits = []
    for pat in patterns:
        hits.extend(BASE.glob(pat))
    hits = sorted({h for h in hits if h.is_file()})
    if not hits:
        sys.exit(f"[错误] 在 {BASE} 找不到 {label} 文件，请用命令行参数显式指定。")
    # 优先「原始」文件：若它是 V3000，ensure_v2000 会复用同名的 -v2000.sdf
    orig = [h for h in hits if "v2000" not in h.name.lower()]
    v2000 = [h for h in hits if "v2000" in h.name.lower()]
    if len(orig) == 1:
        return orig[0]
    if not orig and len(v2000) == 1:
        return v2000[0]
    if len(hits) == 1:
        return hits[0]
    sys.exit("[错误] 发现多个候选，请显式指定：\n  " + "\n  ".join(h.name for h in hits))


def find_obabel():
    exe = shutil.which("obabel") or shutil.which("obabel.exe")
    if exe:
        return exe
    for pat in (r"C:/Program Files/OpenBabel*/obabel.exe",
                r"C:/Program Files (x86)/OpenBabel*/obabel.exe",
                r"C:/OpenBabel*/obabel.exe"):
        found = glob.glob(pat)
        if found:
            return found[0]
    return None


def ensure_v2000(sdf_path, use_obabel):
    """若是 V3000，转成 V2000 并返回 (实际使用的路径, 说明文字)。"""
    text = sdf_path.read_text(encoding="utf-8", errors="replace")
    if "V3000" not in text:
        return sdf_path, "V2000"
    out = sdf_path.with_name(sdf_path.stem + "-v2000.sdf")
    if out.is_file() and "V3000" not in out.read_text(encoding="utf-8", errors="replace"):
        print(f"[信息] 检测到 V3000，复用已转换文件 {out.name}")
        return out, "V2000（Open Babel 预转换）"
    if use_obabel:
        exe = find_obabel()
        if exe:
            cmd = [exe, "-isdf", sdf_path.name, "-osdf", "-O", out.name]
            print("[信息] 检测到 V3000，调用 Open Babel 转换：" + " ".join(cmd))
            try:
                r = subprocess.run(cmd, cwd=str(sdf_path.parent), capture_output=True,
                                   text=True, timeout=300)
                if r.returncode == 0 and out.is_file():
                    print("        " + (r.stderr or r.stdout).strip().splitlines()[-1])
                    return out, "V3000→V2000（Open Babel）"
                print(f"[警告] Open Babel 失败：{(r.stderr or r.stdout).strip()[:300]}", file=sys.stderr)
            except Exception as exc:
                print(f"[警告] Open Babel 调用失败：{exc}", file=sys.stderr)
        else:
            print("[警告] 未找到 Open Babel，改用内置转换器（建议安装 Open Babel 3.1.1）",
                  file=sys.stderr)
    return sdf_path, "V3000→V2000（内置转换器）"


# ---------------------------------------------------------------- PDB
def read_protein(path, name_override):
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.endswith("\n"):
        text += "\n"
    atoms = []
    n_water = 0
    for line in text.splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == "HOH":
            n_water += 1
        if not line.startswith("ATOM  "):   # 口袋残基只取聚合物原子，排除 HETATM 水/配体
            continue
        elem = line[76:78].strip() or line[12:16].strip()[:1]
        if elem.upper() in ("H", "D"):
            continue
        atoms.append({
            "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            "chain": line[21].strip(), "seq": line[22:26].strip(),
            "name": line[17:20].strip(), "icode": line[26].strip(),
        })
    if not atoms:
        sys.exit("[错误] PDB 里没有可用的蛋白重原子（ATOM 记录）。")
    if name_override:
        title = name_override
    else:
        raw = " ".join(l[10:].strip() for l in text.splitlines() if l.startswith("TITLE"))
        raw = re.sub(r"\s+", " ", raw).strip()
        title = f"{path.stem.upper()} · {raw}" if raw else path.stem.upper()
        if len(title) > 110:
            title = title[:107] + "…"
    return {"name": title, "pdb": text, "nProteinAtoms": len(atoms),
            "nWater": n_water, "source": path.name}, atoms


def pocket_of(heavy, protein_atoms, cutoff):
    near, c2 = set(), cutoff * cutoff
    for pa in protein_atoms:
        px, py, pz = pa["xyz"]
        for lx, ly, lz in heavy:
            dx = px - lx
            if abs(dx) > cutoff:
                continue
            dy = py - ly
            if abs(dy) > cutoff:
                continue
            dz = pz - lz
            if abs(dz) > cutoff:
                continue
            if dx * dx + dy * dy + dz * dz <= c2:
                near.add((pa["chain"], pa["seq"], pa["icode"], pa["name"]))
                break

    def key(r):
        m = re.match(r"(-?\d+)", r[1])
        return (r[0], int(m.group(1)) if m else 10 ** 9, r[1], r[2])

    display, items = [], []
    for chain, seq, icode, resname in sorted(near, key=key):
        display.append(f"{resname} {seq}{icode} ({chain})" if icode
                       else f"{resname} {seq} ({chain})")
        m = re.fullmatch(r"(-?\d+)", seq)
        if not m:
            continue
        item = {"auth_asym_id": chain, "auth_seq_id": int(m.group(1))}
        if icode:
            item["pdbx_PDB_ins_code"] = icode
        items.append(item)
    return display, items


# ---------------------------------------------------------------- SDF
def split_records(text):
    for chunk in text.split("$$$$"):
        lines = chunk.splitlines()
        while lines and not lines[0].strip():
            lines.pop(0)
        while lines and not lines[-1].strip():
            lines.pop()
        if len(lines) >= 4:
            yield lines


def parse_data_fields(lines):
    """读取 `>  <name>` 之后的单行属性值（跳过多行的 fieldpoints 大块）。"""
    props, i = {}, 0
    while i < len(lines):
        if lines[i].startswith("> "):
            m = re.search(r"<\s*([^>]+?)\s*>", lines[i])
            key = m.group(1).strip() if m else ""
            vals = []
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith("> "):
                vals.append(lines[i].strip())
                i += 1
            if key and not PROP_SKIP.search(key) and len(vals) == 1 and len(vals[0]) <= 120:
                props[key] = vals[0]
            continue
        i += 1
    return props


def convert_v3000_record(lines):
    """纯 Python 兜底：把一条 V3000 记录转成等价 V2000 molfile。"""
    m = re.search(r"M\s+V30\s+COUNTS\s+(\d+)\s+(\d+)", "\n".join(lines))
    if not m:
        raise ValueError("V3000 缺少 COUNTS 行")
    na, nb = int(m.group(1)), int(m.group(2))
    atoms, bonds, section = [], [], None
    for line in lines:
        s = line.strip()
        if not s.startswith("M  V30"):
            continue
        body = s[6:].strip()
        if body.startswith("BEGIN ATOM"):
            section = "atom"; continue
        if body.startswith("END ATOM"):
            section = None; continue
        if body.startswith("BEGIN BOND"):
            section = "bond"; continue
        if body.startswith("END BOND"):
            section = None; continue
        if body.startswith(("BEGIN ", "END ", "COUNTS")):
            section = None; continue
        parts = body.split()
        # 普通原子行只有 6 个 token：idx sym x y z aamap；带 CHG=/CFG= 时更多
        if section == "atom" and len(parts) >= 6 and parts[0].isdigit():
            sym = parts[1]
            try:
                xyz = (float(parts[2]), float(parts[3]), float(parts[4]))
            except ValueError:
                continue
            chg, parity = 0, 0
            for tok in parts[6:]:          # parts[5] 是 aamap
                if tok.upper().startswith("CHG="):
                    chg = int(tok.split("=", 1)[1])
                elif tok.upper().startswith("CFG="):
                    v = int(tok.split("=", 1)[1])
                    parity = v if v in (1, 2) else 0
            atoms.append((int(parts[0]), sym, xyz, chg, parity))
        elif section == "bond" and len(parts) >= 4 and parts[0].isdigit():
            try:
                bonds.append((int(parts[1]), int(parts[2]), int(parts[3])))
            except ValueError:
                pass
    if len(atoms) != na or len(bonds) != nb:
        raise ValueError(f"V3000 原子/键数量不符：声明 {na}/{nb}，实读 {len(atoms)}/{len(bonds)}")

    out = [lines[0].strip(), lines[1].strip(), "",
           f"{na:3d}{nb:3d}  0  0  0  0  0  0  0  0999 V2000"]
    for _, sym, (x, y, z), chg, parity in atoms:
        cc = CHARGE_CODE.get(chg, 0)
        out.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {sym:<3s}" +
                   "".join(f"{v:3d}" for v in (0, cc, parity, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
    for order, a, b in bonds:
        out.append(f"{a:3d}{b:3d}{order:3d}{0:3d}{0:3d}{0:3d}{0:3d}")
    charges = [(i + 1, a[3]) for i, a in enumerate(atoms) if a[3]]
    for start in range(0, len(charges), 8):
        grp = charges[start:start + 8]
        out.append("M  CHG" + f"{len(grp):3d}" +
                   "".join(f"{idx:4d}{chg:4d}" for idx, chg in grp))
    out.append("M  END")
    # 保留少量有用的单行属性，丢弃 Cresset fieldpoints 大块
    props = parse_data_fields(lines)
    for k, v in props.items():
        out += [f">  <{k}>", v, ""]
    parsed = [{"xyz": xyz, "elem": sym, "chg": chg}
              for _, sym, xyz, chg, _ in atoms]
    return out, parsed


def parse_record(lines):
    """返回 (单条 SDF 文本, 原子列表, 来源说明)。"""
    ver = lines[3].strip()
    if ver.endswith("V3000"):
        body, atoms = convert_v3000_record(lines)
        return "\n".join(body) + "\n$$$$\n", atoms, "内置转换器"
    if "V2000" not in ver:
        raise ValueError(f"既不是 V2000 也不是 V3000：{ver[:40]!r}")
    na, nb = int(lines[3][0:3]), int(lines[3][3:6])
    if len(lines) < 4 + na + nb:
        raise ValueError("V2000 记录行数不足")
    atoms = [{"xyz": (float(l[0:10]), float(l[10:20]), float(l[20:30])),
              "elem": l[31:34].strip() or "C", "chg": 0}
             for l in lines[4:4 + na]]
    return "\n".join(lines).rstrip() + "\n$$$$\n", atoms, "V2000"


def formula_of(atoms):
    counts = defaultdict(int)
    for a in atoms:
        counts[a["elem"].upper()] += 1
    order = {"C": 0, "H": 1}
    return "".join(k if v == 1 else f"{k}{v}"
                   for k, v in sorted(counts.items(),
                                      key=lambda kv: (order.get(kv[0], 2), kv[0])))


# ---------------------------------------------------------------- 主流程
def main():
    args = parse_args()
    pdb_path = resolve(args.pdb, ["*.pdb", "*.ent"], "受体")
    sdf_input = resolve(args.sdf, ["*.sdf", "*.sd", "*.mol"], "配体")
    sdf_path, ver_note = ensure_v2000(sdf_input, not args.no_obabel)
    out_path = _abs(args.out) if args.out else BASE / "molstar_viewer" / "data.js"

    protein, protein_atoms = read_protein(pdb_path, args.name)
    sdf_text = sdf_path.read_text(encoding="utf-8", errors="replace")

    records = []
    skipped = 0
    for idx, lines in enumerate(split_records(sdf_text), 1):
        title = lines[0].strip() or f"mol_{idx}"
        try:
            sdf, atoms, _ = parse_record(lines)
        except Exception as exc:
            if args.strict:
                raise
            print(f"[警告] 跳过第 {idx} 条 {title!r}：{exc}", file=sys.stderr)
            skipped += 1
            continue
        heavy = [a["xyz"] for a in atoms if a["elem"].upper() not in ("H", "D")]
        centroid = [round(sum(p[i] for p in heavy) / max(1, len(heavy)), 3) for i in range(3)]
        pocket, pocket_items = pocket_of(heavy, protein_atoms, args.cutoff)
        group_key = title
        if args.group_regex:
            m = re.search(args.group_regex, title)
            group_key = m.group(1) if m else title
        records.append({
            "name": title, "groupKey": group_key,
            "nAtoms": len(atoms), "nHeavy": len(heavy),
            "formula": formula_of(atoms), "centroid": centroid,
            "pocket": pocket, "pocketItems": pocket_items,
            "props": {} if args.no_props else parse_data_fields(lines),
            "sdf": sdf,
        })
    if not records:
        sys.exit("[错误] SDF 中没有解析出任何有效记录。")

    groups = OrderedDict()
    for i, r in enumerate(records):
        g = groups.setdefault(r["groupKey"], [])
        r["poseNo"] = len(g) + 1
        g.append(i)

    score_label = ""
    if not args.no_props:
        def numeric(key):
            vals = {}
            for i, r in enumerate(records):
                try:
                    vals[i] = float(str(r["props"][key]))
                except (KeyError, TypeError, ValueError):
                    pass
            return vals

        if args.score_field:
            vals = numeric(args.score_field)
            if vals:
                score_label = args.score_field
                for i, v in vals.items():
                    records[i]["score"] = v
            else:
                print(f"[警告] 打分属性 {args.score_field!r} 不是数值或不存在，已忽略。",
                      file=sys.stderr)
        else:
            # 按覆盖率选最有用的打分；混合数据集里部分记录缺字段是允许的
            all_keys = list(SCORE_CANDIDATES)
            for r in records:
                for k in r["props"]:
                    if k not in all_keys and not SCORE_BLACKLIST.search(k):
                        all_keys.append(k)
            best_key, best_vals, best_cov = "", {}, 0
            for key in all_keys:            # 覆盖率优先；覆盖率相同则取优先级更高的字段
                vals = numeric(key)
                if len(vals) > best_cov:
                    best_key, best_vals, best_cov = key, vals, len(vals)
            if best_vals:
                score_label = best_key
                for i, v in best_vals.items():
                    records[i]["score"] = v
                print(f"[信息] 打分字段 {best_key!r}：{best_cov}/{len(records)} 个 pose 有值")

    data = {
        "protein": protein, "ligandSource": sdf_input.name,
        "sdfVersion": ver_note, "cutoff": args.cutoff, "scoreName": score_label,
        "poses": records,
        "groups": [{"name": k, "poseIds": v} for k, v in groups.items()],
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        "// 由 prep_data.py 自动生成；结构数据全部内嵌，无需联网。\n"
        f"// 受体 {pdb_path.name} | 配体 {sdf_input.name} | "
        f"{len(records)} poses / {len(groups)} 分子\n"
        "window.VIEWER_DATA = " +
        json.dumps(data, ensure_ascii=False, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )

    print(f"\n[完成] {out_path}")
    print(f"  受体     : {pdb_path.name} · {protein['nProteinAtoms']} 个蛋白重原子")
    print(f"            {protein['name']}")
    print(f"  配体     : {sdf_input.name} → {sdf_path.name}（{ver_note}）")
    print(f"            {len(records)} 个 pose / {len(groups)} 个分子"
          + (f"，跳过 {skipped} 条" if skipped else ""))
    print(f"  口袋阈值 : {args.cutoff} Å；首个 pose 口袋残基 {len(records[0]['pocket'])} 个")
    print(f"  打分字段 : {score_label or '（无）'}")
    print(f"  data.js  : {out_path.stat().st_size:,} bytes")
    print("  下一步   : 打开 molstar_viewer/index.html（或打包整个 molstar_viewer 文件夹）")


if __name__ == "__main__":
    main()
