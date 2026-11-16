"""Draw actual products with colours traced through their parent atom mappings."""

import colorsys
import html
import json
import math
from pathlib import Path

import cairosvg
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D

from fragment import FragmentConfig
from fragment.modern.chemistry.input import read_parent
from fragment.modern.chemistry.identity import identify
from fragment.modern.chemistry.products import cap

COLORS = [
    (0.73, 0.85, 0.98),
    (0.99, 0.80, 0.63),
    (0.73, 0.88, 0.75),
    (0.87, 0.77, 0.94),
    (0.98, 0.90, 0.60),
    (0.68, 0.88, 0.88),
]
INK = "#253247"
WIDTH = 1600


def block_label(index):
    label = ""
    while index >= 0:
        label = chr(65 + index % 26) + label
        index = index // 26 - 1
    return label


def draw_figure(record, result, occurrences, output):
    parent, _ = read_parent(record, FragmentConfig(**result["config"]))
    if identify(parent)["identity_key"] != result["parent"]["identity_key"]:
        raise ValueError("Drawing input differs from the processed parent.")
    blocks = [set(b["components"]) for b in result["elementary_fragments"]]
    colors = list(COLORS)
    while len(colors) < len(blocks):
        colors.append(colorsys.hsv_to_rgb((len(colors) * 0.618) % 1, 0.25, 0.96))
    owner = {a: i for i, block in enumerate(blocks) for a in block}
    elementary, combinations, skipped = [], [], []
    seen = set()
    all_atoms = set(range(parent.GetNumAtoms()))
    for item in occurrences:
        construction = item["construction"]
        members = frozenset(construction["parent_atom_indices"])
        if members == all_atoms:
            continue  # The parent is already shown at the top.
        indices = [i for i, block in enumerate(blocks) if block <= members]
        whole = set().union(*(blocks[i] for i in indices))
        # The diagram illustrates whole elementary blocks joined along original
        # connections. Ring-specific internal edits remain in the full results.
        if whole != set(members) or construction["internal_bond_edits"]:
            skipped.append(item["identity_key"])
            continue
        key = (members, item["identity_key"])
        if key in seen:
            continue
        seen.add(key)
        product, rebuilt = cap(parent, members)
        if identify(product)["identity_key"] != item["identity_key"]:
            raise ValueError("Drawn product does not match the exported identity.")
        entry = dict(
            label=" + ".join(block_label(i) for i in indices),
            indices=indices,
            mol=product,
            mapping=rebuilt["product_to_parent"],
            identity_key=item["identity_key"],
            smiles=item["smiles"],
        )
        (elementary if len(indices) == 1 else combinations).append(entry)
    # A single elementary block can equal the entire parent; show it in panel b.
    if len(blocks) == 1:
        elementary = [
            dict(
                label="A",
                indices=[0],
                mol=parent,
                mapping=list(range(parent.GetNumAtoms())),
                identity_key=result["parent"]["identity_key"],
                smiles=result["parent"]["smiles"],
            )
        ]
    elementary.sort(key=lambda item: item["indices"])
    combinations.sort(key=lambda item: (len(item["indices"]), item["indices"]))
    if len(elementary) != len(blocks):
        raise ValueError(
            "Not all elementary fragments can be displayed; inspect the full results."
        )

    parent_height = max(320, 100 + len(blocks) * 38)
    first_arrow = 165 + parent_height + 10
    elem_title = first_arrow + 50
    elem_top = elem_title + 45
    elem_bottom = elem_top + math.ceil(len(elementary) / 4) * 260
    second_arrow = elem_bottom + 10
    combination_title = second_arrow + 50
    combination_top = combination_title + 45
    height = combination_top + max(1, math.ceil(len(combinations) / 3)) * 280 + 70
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" viewBox="0 0 {WIDTH} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]

    def text(x, y, value, size=24, bold=False, anchor="start", color=INK):
        svg.append(
            f'<text x="{x}" y="{y + size * .8}" font-family="sans-serif" font-size="{size}" '
            f'font-weight="{700 if bold else 400}" text-anchor="{anchor}" fill="{color}">{html.escape(str(value))}</text>'
        )

    def line(x1, y1, x2, y2, color="#CBD2DA", width=2):
        svg.append(
            f'<path d="M{x1},{y1} L{x2},{y2}" stroke="{color}" stroke-width="{width}"/>'
        )

    def arrow(y, caption):
        line(600, y, 600, y + 30, INK)
        line(600, y + 30, 593, y + 22, INK)
        line(600, y + 30, 607, y + 22, INK)
        text(625, y + 6, caption, 20)

    def molecule(mol, mapping, x, y, width, cell_height, cuts=False):
        mol = Chem.Mol(mol)
        for atom, origin in zip(mol.GetAtoms(), mapping):
            atom.SetIntProp("_figure_parent", -1 if origin is None else origin)
        mol = Chem.RemoveHs(mol)
        rdDepictor.Compute2DCoords(mol)
        atom_colors = {
            a.GetIdx(): colors[owner[a.GetIntProp("_figure_parent")]]
            for a in mol.GetAtoms()
            if a.GetIntProp("_figure_parent") >= 0
        }
        bond_colors, boundaries = {}, []
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if atom_colors.get(a) == atom_colors.get(b) and a in atom_colors:
                bond_colors[bond.GetIdx()] = atom_colors[a]
            elif cuts:
                boundaries.append((a, b))
        draw = rdMolDraw2D.MolDraw2DSVG(width, cell_height)
        opts = draw.drawOptions()
        opts.useBWAtomPalette()
        opts.padding = 0.12
        opts.fixedBondLength = 34
        opts.minFontSize = 16
        opts.maxFontSize = 24
        opts.highlightBondWidthMultiplier = 16
        draw.DrawMolecule(
            mol,
            highlightAtoms=list(atom_colors),
            highlightBonds=list(bond_colors),
            highlightAtomColors=atom_colors,
            highlightBondColors=bond_colors,
            highlightAtomRadii={a: 0.32 for a in atom_colors},
        )
        draw.SetColour((0.78, 0.22, 0.20))
        draw.SetLineWidth(2.5)
        for a, b in boundaries:
            pa, pb = draw.GetDrawCoords(a), draw.GetDrawCoords(b)
            dx, dy = pb.x - pa.x, pb.y - pa.y
            length = math.hypot(dx, dy)
            cx, cy = (pa.x + pb.x) / 2, (pa.y + pb.y) / 2
            draw.DrawLine(
                Point2D(cx - dy / length * 10, cy + dx / length * 10),
                Point2D(cx + dy / length * 10, cy - dx / length * 10),
                True,
            )
        draw.FinishDrawing()
        drawing = draw.GetDrawingText()
        drawing = drawing[drawing.index("<svg") :]
        svg.append(f'<g transform="translate({x},{y})">{drawing}</g>')

    text(55, 28, "Fragmentation and fragment combinations", 32, True)
    text(
        55,
        78,
        f'{result["id"]}  |  {result["backend"]}  |  {Chem.RemoveHs(parent).GetNumAtoms()} heavy atoms',
        22,
    )
    line(55, 115, 1545, 115)
    text(55, 138, "a  Parent structure", 27, True)
    molecule(
        parent,
        list(range(parent.GetNumAtoms())),
        440,
        175,
        640,
        parent_height - 10,
        True,
    )
    text(1150, 200, "Elementary regions", 21, True)
    for i in range(len(blocks)):
        color = "#" + "".join(f"{round(v * 255):02x}" for v in colors[i])
        svg.append(
            f'<rect x="1150" y="{240 + 38 * i}" width="22" height="22" rx="5" fill="{color}"/>'
        )
        text(1186, 239 + 38 * i, block_label(i), 22, True)
    text(1150, 255 + 38 * len(blocks), "Red marks: cut bonds", 19, color="#C73833")
    arrow(first_arrow, "Partition and cap exposed boundaries")
    text(55, elem_title, "b  Elementary fragments", 27, True)
    for i, item in enumerate(elementary):
        row, col = divmod(i, 4)
        n = min(4, len(elementary) - row * 4)
        x, y = (WIDTH - n * 380) // 2 + col * 380, elem_top + row * 260
        molecule(item["mol"], item["mapping"], x, y, 360, 210)
        text(x + 180, y + 218, item["label"], 25, True, "middle")
    arrow(second_arrow, "Combine adjacent regions; cap outer boundaries")
    text(55, combination_title, "c  Combined fragments", 27, True)
    if not combinations:
        text(
            800,
            combination_top + 100,
            "No additional non-parent combinations.",
            24,
            anchor="middle",
        )
    for i, item in enumerate(combinations):
        row, col = divmod(i, 3)
        n = min(3, len(combinations) - row * 3)
        x, y = (WIDTH - n * 510) // 2 + col * 510, combination_top + row * 280
        molecule(item["mol"], item["mapping"], x, y, 490, 230)
        text(x + 245, y + 238, item["label"], 25, True, "middle")
    text(
        800,
        height - 42,
        "Colours track parent regions; exposed boundaries are hydrogen-capped.",
        20,
        anchor="middle",
    )
    svg.append("</svg>")
    data = "\n".join(svg).encode("utf-8")
    output = Path(output)
    (output / "figure2.svg").write_bytes(data)
    # PDF is vector graphics, not a bitmap embedded in a PDF page.
    cairosvg.svg2pdf(bytestring=data, write_to=str(output / "figure2.pdf"))
    cairosvg.svg2png(bytestring=data, write_to=str(output / "figure2.png"))
    metadata = dict(
        backend=result["backend"],
        parent_id=result["id"],
        elementary=[
            {k: v for k, v in x.items() if k not in ("mol", "mapping")}
            for x in elementary
        ],
        combinations=[
            {k: v for k, v in x.items() if k not in ("mol", "mapping")}
            for x in combinations
        ],
        excluded_from_figure=sorted(set(skipped)),
        note="Parent-equivalent products appear only in the parent panel. Repeated families sharing the same atom membership and identity are drawn once. Full exports are unchanged.",
    )
    (output / "figure2.json").write_text(json.dumps(metadata, indent=2) + "\n")
