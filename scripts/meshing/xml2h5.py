#!/usr/bin/env python3
"""
xml2h5.py -- convert a cardiax mesh XML into an HDF5 file containing ONLY the
mesh (geometry + topology) and the fibers.

Output layout
--------------
/mesh
    /geometry        (n_points, 3)  float64   node coordinates
    /topology        (n_elements, nodes_per_element) int32   connectivity
    /fiber           (n_elements, 3) float64   fiber direction  (if present)
    /sheet           (n_elements, 3) float64   sheet/transverse (if present)
    /normal          (n_elements, 3) float64   normal           (if present)
    attrs: celltype, dim, fiber_model

/cell_data           empty group, ready to be filled (per-element fields)
/point_data          empty group, ready to be filled (per-node fields)

Usage
-----
    python xml2h5.py mesh.xml [mesh.h5]

If the output name is omitted, the input name with a .h5 suffix is used.
"""

import sys
import os
import xml.etree.ElementTree as ET

import numpy as np
import h5py


def _parse_vec(text):
    """'1.0,0.0,0.0' -> np.array([1.0, 0.0, 0.0])"""
    return np.fromstring(text.strip(), sep=",")


def read_xml(path):
    tree = ET.parse(path)
    mesh = tree.getroot()
    if mesh.tag != "mesh":
        raise ValueError(f"root element is <{mesh.tag}>, expected <mesh>")

    celltype = mesh.get("celltype", "unknown")
    dim = int(mesh.get("dim", "3"))

    # ---- nodes ------------------------------------------------------------
    nodes_el = mesh.find("nodes")
    n_points = int(nodes_el.get("size"))
    geometry = np.zeros((n_points, 3), dtype=np.float64)
    for node in nodes_el.findall("node"):
        i = int(node.get("id"))
        geometry[i, 0] = float(node.get("x"))
        geometry[i, 1] = float(node.get("y"))
        geometry[i, 2] = float(node.get("z"))

    # ---- elements ---------------------------------------------------------
    elems_el = mesh.find("elements")
    n_elements = int(elems_el.get("size"))
    elem_list = elems_el.findall("element")

    # infer nodes-per-element from the vN attributes of the first element
    npe = sum(1 for k in elem_list[0].keys() if k.startswith("v"))
    topology = np.zeros((n_elements, npe), dtype=np.int32)
    for el in elem_list:
        i = int(el.get("id"))
        for j in range(npe):
            topology[i, j] = int(el.get(f"v{j}"))

    # ---- fibers (per element) --------------------------------------------
    fiber_model = None
    fiber = sheet = normal = None
    edata = mesh.find("element_data")
    if edata is not None:
        fiber_model = edata.get("type")
        fiber = np.zeros((n_elements, 3), dtype=np.float64)
        has_sheet = has_normal = False
        sheet = np.zeros((n_elements, 3), dtype=np.float64)
        normal = np.zeros((n_elements, 3), dtype=np.float64)
        for el in edata.findall("element"):
            i = int(el.get("id"))
            f = el.find("fiber")
            if f is not None and f.text:
                fiber[i] = _parse_vec(f.text)
            s = el.find("sheet")
            if s is not None and s.text:
                sheet[i] = _parse_vec(s.text)
                has_sheet = True
            n = el.find("normal")
            if n is not None and n.text:
                normal[i] = _parse_vec(n.text)
                has_normal = True
        if not has_sheet:
            sheet = None
        if not has_normal:
            normal = None

    return dict(
        celltype=celltype,
        dim=dim,
        fiber_model=fiber_model,
        geometry=geometry,
        topology=topology,
        fiber=fiber,
        sheet=sheet,
        normal=normal,
    )


def write_h5(data, path):
    with h5py.File(path, "w") as h5:
        mesh = h5.create_group("mesh")
        mesh.attrs["celltype"] = data["celltype"]
        mesh.attrs["dim"] = data["dim"]
        if data["fiber_model"] is not None:
            mesh.attrs["fiber_model"] = data["fiber_model"]

        mesh.create_dataset("geometry", data=data["geometry"],
                            compression="gzip")
        mesh.create_dataset("topology", data=data["topology"],
                            compression="gzip")

        if data["fiber"] is not None:
            mesh.create_dataset("fiber", data=data["fiber"],
                                compression="gzip")
        if data["sheet"] is not None:
            mesh.create_dataset("sheet", data=data["sheet"],
                                compression="gzip")
        if data["normal"] is not None:
            mesh.create_dataset("normal", data=data["normal"],
                                compression="gzip")

        # empty groups, to be filled by the solver later
        h5.create_group("cell_data")
        h5.create_group("point_data")


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        sys.exit(1)

    in_path = argv[1]
    out_path = argv[2] if len(argv) > 2 else os.path.splitext(in_path)[0] + ".h5"

    data = read_xml(in_path)
    write_h5(data, out_path)

    print(f" wrote {out_path}")
    print(f"   celltype    : {data['celltype']}")
    print(f"   dim         : {data['dim']}")
    print(f"   points      : {data['geometry'].shape[0]}")
    print(f"   elements    : {data['topology'].shape[0]} "
          f"({data['topology'].shape[1]} nodes/element)")
    print(f"   fiber_model : {data['fiber_model']}")
    print(f"   fiber       : {'yes' if data['fiber']  is not None else 'no'}")
    print(f"   sheet       : {'yes' if data['sheet']  is not None else 'no'}")
    print(f"   normal      : {'yes' if data['normal'] is not None else 'no'}")


if __name__ == "__main__":
    main(sys.argv)
