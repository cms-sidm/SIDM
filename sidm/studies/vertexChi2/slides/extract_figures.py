#!/usr/bin/env python3
"""Write the PNG outputs of the write-up's code cells to figures/, named by section.

Each image is the notebook's own rendered output (outputs[].data["image/png"]), decoded
unchanged; nothing is re-plotted. A code cell belongs to the section of the last markdown
heading above it ("## 8." or "### 10.4"). Files are named sec<section>_cell<index>[_<n>].png,
e.g. sec10.4_cell70.png.
"""
import argparse
import base64
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
HEADING = re.compile(r"^#{2,3}\s+(\d+(?:\.\d+)?)\.?\s", re.M)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--notebook", default=os.path.join(HERE, "..", "vertexChi2CosmicVetoStudy.ipynb"))
    ap.add_argument("--outdir", default=os.path.join(HERE, "figures"))
    ap.add_argument("--sections", nargs="+", default=["3", "7", "8", "9", "10"],
                    help="top-level sections to extract")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    with open(args.notebook) as fh:
        nb = json.load(fh)
    section = None
    for idx, cell in enumerate(nb["cells"]):
        src = "".join(cell["source"])
        if cell["cell_type"] == "markdown":
            found = HEADING.findall(src)
            if found:
                section = found[-1]
            continue
        if section is None or section.split(".")[0] not in args.sections:
            continue
        pngs = [o["data"]["image/png"] for o in cell.get("outputs", []) if "image/png" in o.get("data", {})]
        for n, png in enumerate(pngs):
            suffix = f"_{n}" if len(pngs) > 1 else ""
            path = os.path.join(args.outdir, f"sec{section}_cell{idx}{suffix}.png")
            data = png if isinstance(png, str) else "".join(png)
            with open(path, "wb") as fh:
                fh.write(base64.b64decode(data))
            print(path)


if __name__ == "__main__":
    main()
