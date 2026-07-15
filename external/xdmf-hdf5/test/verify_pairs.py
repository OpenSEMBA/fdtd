#!/usr/bin/env python3
"""Black-box validation for generated XDMF/HDF5 pairs."""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import numpy as np


PAIR_NAMES = (
    "uniform",
    "volume",
    "rectilinear",
    "curvilinear",
    "unstructured",
    "mixed",
    "time-series",
    "frequency-series",
    "escaped-&-pair",
    "validation-static",
    "validation-series",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def text(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def dimensions(item: ET.Element) -> tuple[int, ...]:
    return tuple(int(value) for value in item.attrib["Dimensions"].split())


def validate_pair(root: Path, name: str) -> tuple[ET.Element, h5py.File]:
    xdmf_path = root / f"{name}.xdmf"
    hdf5_path = root / f"{name}.h5"
    require(xdmf_path.is_file(), f"missing {xdmf_path.name}")
    require(hdf5_path.is_file(), f"missing {hdf5_path.name}")

    document = ET.parse(xdmf_path).getroot()
    require(document.tag == "Xdmf", f"{name}: invalid XML root")
    require(document.attrib.get("Version") == "3.0", f"{name}: wrong XDMF version")
    hdf5 = h5py.File(hdf5_path, "r")
    require(text(hdf5.attrs["schema_name"]) == "OpenSEMBA XDMF-HDF5", f"{name}: schema name")
    require(text(hdf5.attrs["schema_version"]) == "1.0", f"{name}: schema version")

    for item in document.findall(".//DataItem[@Format='HDF']"):
        reference = (item.text or "").strip()
        require(":" in reference, f"{name}: malformed HDF reference {reference!r}")
        referenced_file, dataset_path = reference.split(":", 1)
        require(referenced_file == hdf5_path.name, f"{name}: external HDF reference")
        require(dataset_path in hdf5, f"{name}: missing dataset {dataset_path}")
        dataset = hdf5[dataset_path]
        require(dataset.shape == dimensions(item), f"{name}:{dataset_path} shape {dataset.shape}")

        number_type = item.attrib.get("NumberType")
        precision = int(item.attrib.get("Precision", "0"))
        if dataset.dtype.kind == "f":
            require(number_type == "Float", f"{name}:{dataset_path} NumberType")
        elif dataset.dtype.kind in "iu":
            require(number_type == "Int", f"{name}:{dataset_path} NumberType")
        else:
            raise AssertionError(f"{name}:{dataset_path} unsupported dtype {dataset.dtype}")
        require(dataset.dtype.itemsize == precision, f"{name}:{dataset_path} precision")

    for slab in document.findall(".//DataItem[@ItemType='HyperSlab']"):
        children = slab.findall("./DataItem")
        require(len(children) == 2, f"{name}: HyperSlab child count")
        selector, source = children
        source_shape = dimensions(source)
        selector_shape = dimensions(selector)
        require(selector_shape == (3, len(source_shape)), f"{name}: HyperSlab selector shape")
        selector_values = [int(value) for value in (selector.text or "").split()]
        require(len(selector_values) == 3 * len(source_shape), f"{name}: HyperSlab selector size")
        rank = len(source_shape)
        start = selector_values[:rank]
        stride = selector_values[rank : 2 * rank]
        count = selector_values[2 * rank :]
        require(all(value >= 0 for value in start), f"{name}: HyperSlab start")
        require(all(value > 0 for value in stride), f"{name}: HyperSlab stride")
        require(all(value > 0 for value in count), f"{name}: HyperSlab count")
        require(
            all(start[i] + (count[i] - 1) * stride[i] < source_shape[i] for i in range(rank)),
            f"{name}: HyperSlab outside source",
        )
        require(tuple(count[1:]) == dimensions(slab), f"{name}: HyperSlab output shape")

    for topology in document.findall(".//Topology"):
        data_item = topology.find("./DataItem[@Format='HDF']")
        if data_item is None or topology.attrib.get("TopologyType") == "Mixed":
            continue
        element_count = int(topology.attrib["NumberOfElements"])
        nodes_per_element = int(topology.attrib["NodesPerElement"])
        require(
            int(np.prod(dimensions(data_item))) == element_count * nodes_per_element,
            f"{name}: topology connectivity cardinality",
        )

    return document, hdf5


def verify_uniform(document: ET.Element, hdf5: h5py.File) -> None:
    topology = document.find(".//Topology")
    require(topology is not None, "uniform: missing topology")
    require(topology.attrib["Dimensions"] == "2 3", "uniform: topology dimensions")
    values = hdf5["/attributes/a0001/values"][:]
    require(values.shape == (2, 3), "uniform: external axis order")
    np.testing.assert_allclose(values.ravel(), np.arange(1, 7) * 0.25)
    require(values.dtype == np.dtype("float32"), "uniform: real32 data")
    require(hdf5["/grids/g0001/origin"].dtype == np.dtype("float32"), "uniform: real32 geometry")


def verify_volume(document: ET.Element, hdf5: h5py.File) -> None:
    topology = document.find(".//Topology")
    require(topology is not None, "volume: missing topology")
    require(topology.attrib["Dimensions"] == "21 21 21", "volume: topology dimensions")
    attributes = {element.attrib["Name"] for element in document.findall(".//Attribute")}
    require(attributes == {"electric-field-magnitude"}, "volume: scalar attribute")
    values = hdf5["/attributes/a0001/values"]
    require(values.shape == (21, 21, 21), "volume: field shape")
    require(values.dtype == np.dtype("float64"), "volume: field dtype")
    require(float(values[:].min()) >= 0.0, "volume: negative field magnitude")
    require(float(values[:].max()) > 1.0, "volume: field maximum")


def verify_rectilinear(document: ET.Element, hdf5: h5py.File) -> None:
    names = {element.attrib.get("Name") for element in document.findall(".//Grid")}
    require("collection & <escaped>" in names, "rectilinear: escaped collection name")
    require('grid & "quoted"' in names, "rectilinear: escaped grid name")
    attributes = {element.attrib["Name"] for element in document.findall(".//Attribute")}
    require("velocity <xyz>" in attributes, "rectilinear: escaped attribute name")
    np.testing.assert_allclose(hdf5["/grids/g0001/axis_x"][:], [-2.0, 0.5])
    np.testing.assert_allclose(hdf5["/grids/g0001/axis_y"][:], [0.0, 1.0, 4.0])
    np.testing.assert_allclose(hdf5["/grids/g0001/axis_z"][:], [10.0, 11.0, 13.0, 17.0])

    vector = hdf5["/attributes/a0001/values"]
    require(vector.shape == (4, 3, 2, 3), "rectilinear: vector axis order")
    require(vector.dtype == np.dtype("float64"), "rectilinear: vector dtype")
    np.testing.assert_allclose(vector[:].ravel(), np.arange(1, 73))
    material = hdf5["/attributes/a0002/values"]
    require(material.shape == (3, 2, 1), "rectilinear: cell shape")
    require(material.dtype == np.dtype("int32"), "rectilinear: int32 dtype")
    np.testing.assert_array_equal(material[:].ravel(), np.arange(11, 17))
    identifier = hdf5["/attributes/a0003/values"]
    require(identifier.dtype == np.dtype("int64"), "rectilinear: int64 dtype")
    require(int(identifier[0]) == 9876543210, "rectilinear: int64 value")
    tensor = hdf5["/attributes/a0004/values"]
    require(tensor.shape == (4, 3, 2, 3, 3), "rectilinear: tensor shape")
    tensor6 = hdf5["/attributes/a0005/values"]
    require(tensor6.shape == (3, 2, 1, 6), "rectilinear: tensor6 shape")
    matrix = hdf5["/attributes/a0006/values"]
    require(matrix.shape == (1, 2, 3), "rectilinear: matrix component order")
    require(matrix.dtype == np.dtype("float32"), "rectilinear: matrix dtype")
    np.testing.assert_allclose(
        matrix[0],
        [[4001.0, 4003.0, 4005.0], [4002.0, 4004.0, 4006.0]],
    )


def verify_curvilinear(document: ET.Element, hdf5: h5py.File) -> None:
    topology = document.find(".//Topology")
    geometry = document.find(".//Geometry")
    require(topology is not None and topology.attrib["TopologyType"] == "2DSMesh", "curvilinear: topology")
    require(geometry is not None and geometry.attrib["GeometryType"] == "XY", "curvilinear: geometry")
    points = hdf5["/grids/g0001/points"]
    require(points.shape == (3, 2, 2), "curvilinear: structured point shape")
    require(points.dtype == np.dtype("float32"), "curvilinear: point dtype")


def verify_unstructured(document: ET.Element, hdf5: h5py.File) -> None:
    expected = {
        "Polyvertex",
        "Polyline",
        "Polygon",
        "Triangle",
        "Quadrilateral",
        "Tetrahedron",
        "Pyramid",
        "Wedge",
        "Hexahedron",
        "Edge_3",
        "Quadrilateral_9",
        "Triangle_6",
        "Quadrilateral_8",
        "Tetrahedron_10",
        "Pyramid_13",
        "Wedge_15",
        "Wedge_18",
        "Hexahedron_20",
        "Hexahedron_24",
        "Hexahedron_27",
    }
    actual = {topology.attrib["TopologyType"] for topology in document.findall(".//Topology")}
    require(actual == expected, f"unstructured: topology suite mismatch {actual ^ expected}")
    require(len(hdf5["/grids"]) == len(expected), "unstructured: grid count")
    for grid in hdf5["/grids"].values():
        connectivity = grid["connectivity"][:]
        np.testing.assert_array_equal(connectivity.ravel(), np.arange(connectivity.size))
        require(grid["points"].shape == (27, 3), "unstructured: point shape")


def verify_mixed(document: ET.Element, hdf5: h5py.File) -> None:
    topology = document.find(".//Topology")
    require(topology is not None and topology.attrib["TopologyType"] == "Mixed", "mixed: topology")
    require(topology.attrib["NumberOfElements"] == "8", "mixed: element count")
    expected = np.array(
        [
            4, 0, 1, 2,
            5, 0, 1, 3, 2,
            3, 5, 0, 1, 5, 6, 2,
            6, 0, 1, 2, 4,
            1, 2, 6, 7,
            2, 3, 0, 1, 3,
            16, 4,
            3, 0, 1, 2,
            3, 0, 1, 4,
            3, 1, 2, 4,
            3, 2, 0, 4,
            36, 0, 1, 2, 3, 4, 5,
        ],
        dtype=np.int64,
    )
    connectivity = hdf5["/grids/g0001/connectivity"][:]
    require(connectivity.dtype == np.dtype("int64"), "mixed: connectivity dtype")
    np.testing.assert_array_equal(connectivity, expected)


def verify_time_series(document: ET.Element, hdf5: h5py.File) -> None:
    temporal = document.find(".//Grid[@CollectionType='Temporal']")
    require(temporal is not None, "time-series: temporal collection")
    times = [float(element.attrib["Value"]) for element in temporal.findall("./Grid/Time")]
    expected_times = np.linspace(0.0, 1.0, 20)
    np.testing.assert_allclose(times, expected_times, rtol=1.0e-6)
    np.testing.assert_allclose(hdf5["/series/values"][:], expected_times, rtol=1.0e-6)
    attributes = {element.attrib["Name"] for element in document.findall(".//Attribute")}
    require(
        attributes == {"electric-field-magnitude", "electric-field"},
        "time-series: field names",
    )
    magnitude = hdf5["/attributes/a0001/values"]
    electric_field = hdf5["/attributes/a0002/values"]
    require(magnitude.shape == (20, 25, 25, 25), "time-series: scalar shape")
    require(electric_field.shape == (20, 25, 25, 25, 3), "time-series: vector shape")
    require(magnitude.dtype == np.dtype("float32"), "time-series: scalar dtype")
    require(electric_field.dtype == np.dtype("float32"), "time-series: vector dtype")
    magnitude_values = magnitude[:]
    electric_values = electric_field[:]
    require(float(magnitude_values.min()) >= 0.0, "time-series: negative magnitude")
    require(float(magnitude_values.max()) > 0.99, "time-series: pulse maximum")
    np.testing.assert_allclose(magnitude_values, np.abs(electric_values[..., 1]))
    np.testing.assert_array_equal(electric_values[..., 0], 0.0)
    np.testing.assert_array_equal(electric_values[..., 2], 0.0)
    early_peak = int(np.argmax(magnitude_values[0, 12, 12, :]))
    late_peak = int(np.argmax(magnitude_values[-1, 12, 12, :]))
    require(early_peak < late_peak, "time-series: pulse does not move along +x")
    require(magnitude.chunks is not None and magnitude.chunks[0] == 1, "time-series: chunking")


def verify_frequency_series(document: ET.Element, hdf5: h5py.File) -> None:
    frequencies = [
        float(element.attrib["Value"])
        for element in document.findall(".//Information[@Name='Frequency']")
    ]
    np.testing.assert_allclose(frequencies, [1.0e6, 2.0e6])
    np.testing.assert_allclose(hdf5["/series/values"][:], frequencies)
    real_part = hdf5["/attributes/a0001/values"][:]
    imaginary_part = hdf5["/attributes/a0002/values"][:]
    np.testing.assert_allclose(real_part[0].ravel(), [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(real_part[1].ravel(), [2.0, 4.0, 6.0, 8.0])
    np.testing.assert_allclose(imaginary_part, -real_part)


def verify_validation_series(hdf5: h5py.File) -> None:
    require(hdf5["/series/values"].shape == (1,), "validation-series: rollback extent")
    np.testing.assert_allclose(hdf5["/series/values"][:], [1.0])
    np.testing.assert_allclose(hdf5["/attributes/a0001/values"][0].ravel(), [1, 2, 3, 4])
    np.testing.assert_allclose(hdf5["/attributes/a0002/values"][0].ravel(), [-1, -2, -3, -4])


def main() -> int:
    root = Path(sys.argv[1]).resolve()
    require(root.is_dir(), f"generated output directory does not exist: {root}")
    requested_names = tuple(sys.argv[2:]) or PAIR_NAMES
    for name in requested_names:
        require(name in PAIR_NAMES, f"unknown generated pair: {name}")

    opened: list[h5py.File] = []
    try:
        pairs: dict[str, tuple[ET.Element, h5py.File]] = {}
        for name in requested_names:
            document, hdf5 = validate_pair(root, name)
            opened.append(hdf5)
            pairs[name] = (document, hdf5)

        if "uniform" in pairs:
            verify_uniform(*pairs["uniform"])
        if "volume" in pairs:
            verify_volume(*pairs["volume"])
        if "rectilinear" in pairs:
            verify_rectilinear(*pairs["rectilinear"])
        if "curvilinear" in pairs:
            verify_curvilinear(*pairs["curvilinear"])
        if "unstructured" in pairs:
            verify_unstructured(*pairs["unstructured"])
        if "mixed" in pairs:
            verify_mixed(*pairs["mixed"])
        if "time-series" in pairs:
            verify_time_series(*pairs["time-series"])
        if "frequency-series" in pairs:
            verify_frequency_series(*pairs["frequency-series"])
        if "validation-series" in pairs:
            verify_validation_series(pairs["validation-series"][1])
    finally:
        for hdf5 in opened:
            hdf5.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
