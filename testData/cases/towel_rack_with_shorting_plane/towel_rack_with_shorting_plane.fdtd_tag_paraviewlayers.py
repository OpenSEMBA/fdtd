"""
ParaView macro: create separate threshold layers from ugrfdtd VTK output.

Usage
-----
1. Open ParaView and load your VTK file (e.g. the _MAP_*.vtk or current_t/current_f file).
2. Select the VTK source in the Pipeline Browser.
3. Run this macro: Tools > Macros > Add Macro… then click Run Macro.

Two array modes are handled automatically:
  - 'tagnumber'  → used in SLICE CURRENT VTK PROBES (current_t / current_f files)
  - 'mediatype'  → used in MAP VTK PROBES
Both arrays are processed if present.
"""

import colorsys
from paraview.simple import (
    GetActiveSource,
    GetActiveViewOrCreate,
    Threshold,
    Show,
    Hide,
    Render,
    ResetCamera,
)

# ─────────────────────────────────────────────────────────────────────────────
# Layer definitions – (name, lower_threshold, upper_threshold, description)
# ─────────────────────────────────────────────────────────────────────────────

TAGNUMBER_LAYERS = [
    # name                        low       high   description
    ("Tag_FreeSpace_Undesired",  -1e21,    -1e-3,  "Candidates for undesired free-space slots"),
    ("Tag_NodalSources",          0,        63,    "Nodal sources, etc."),
    ("Tag_PEC_layer1",            64,      127,    "PEC@layer1"),
    ("Tag_PEC_layer2",           128,      191,    "PEC@layer2"),
    ("Tag_PEC_layer3",           192,      255,    "PEC@layer3"),
]

MEDIATYPE_LAYERS = [
    # name                                      low      high    description
    ("Media_FreeSpace_Undesired",              -100,    -100,    "Candidates for undesired free-space slots (Surface)"),
    ("Media_PEC_Surface",                       0.0,     0.0,    "PEC (Surface)"),
    ("Media_PEC_Line",                          0.5,     0.5,    "PEC (Line)"),
    ("Media_Dispersive_Line",                   1.5,     1.5,    "Dispersive electric/magnetic isotropic/anisotropic (Line)"),
    ("Media_Dispersive_Surface",              100,     199,      "Dispersive electric/magnetic isotropic/anisotropic (Surface)"),
    ("Media_Dielectric_Line",                   2.5,     2.5,    "Dielectric isotropic or anisotropic (Line)"),
    ("Media_Dielectric_Surface",              200,     299,      "Dielectric isotropic or anisotropic (Surface)"),
    ("Media_SGBC_MIBC_Line",                    3.5,     3.5,    "sgbc/mibc Isotropic/anisotropic Multiport (Line)"),
    ("Media_SGBC_MIBC_Surface",               300,     399,      "sgbc/mibc Isotropic/anisotropic Multiport (Surface)"),
    ("Media_ThinSlot_Line",                     4.5,     4.5,    "Thin slot (Line)"),
    ("Media_YEEadvanced_conformal_Surface",     5.0,     5.0,    "Already_YEEadvanced_byconformal (Surface)"),
    ("Media_YEEadvanced_conformal_Line",        5.5,     5.5,    "Already_YEEadvanced_byconformal (Line)"),
    ("Media_Split_useless_Surface",             6.0,     6.0,    "Split_and_useless (Surface)"),
    ("Media_Split_useless_Line",                6.5,     6.5,    "Split_and_useless (Line)"),
    ("Media_Edge_NotColliding_ThinWires",       7.0,     7.0,    "Edge Not colliding thin wires (Line)"),
    ("Media_ThinWire_Colliding",                8.0,     8.0,    "Thin wire segments colliding with structure (Line)"),
    ("Media_NodalSource_E",                     8.5,     8.5,    "Soft/Hard Nodal CURRENT/FIELD ELECTRIC DENSITY SOURCE (Line)"),
    ("Media_NodalSource_H",                     9.0,     9.0,    "Soft/Hard Nodal CURRENT/FIELD MAGNETIC DENSITY SOURCE (Line)"),
    ("Media_Wire_Ending",                      10,      11,      "LeftEnd/RightEnd/Ending wire segment (Wire)"),
    ("Media_Wire_Intermediate",                20,      20,      "Intermediate wire segment (Wire)"),
    ("Media_Multiwire_Edge_NotColliding",      12,      12,      "Edge Not colliding multiwires (Multiwire)"),
    ("Media_Multiwire_Colliding",              13,      13,      "Multiwire segments colliding with structure (Multiwire)"),
    ("Media_Multiwire_Ending",                 14,      15,      "LeftEnd/RightEnd/Ending multiwire segment (Multiwire)"),
    ("Media_Multiwire_Intermediate",           60,      60,      "Intermediate multiwire segment (Multiwire)"),
    ("Media_ThinSlot_Surface",                400,     499,      "Thin slot (+indexmedium) (Surface)"),
    ("Media_ConformalPEC_Surface",           1000,    1999,      "Conformal Volume PEC (+indexmedium) (Surface)"),
    ("Media_ConformalPEC_Line",              2000,    2999,      "Conformal Volume PEC (+indexmedium) (Line)"),
    ("Media_Other_Line",                      -0.5,    -0.5,     "Other types of media (Line)"),
    ("Media_Other_Surface",                   -1.0,    -1.0,     "Other types of media (Surface)"),
]

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _distinct_colors(n, saturation=0.75, value=0.88):
    """Return *n* perceptually distinct RGB triples spread across the hue wheel."""
    return [list(colorsys.hsv_to_rgb(i / n, saturation, value)) for i in range(n)]


def _cell_array_names(source):
    """Return the list of cell-data array names present in *source*."""
    source.UpdatePipeline()
    cell_info = source.GetDataInformation().GetAttributeInformation(1)  # 1 = cell data
    return [cell_info.GetArrayInformation(i).GetName()
            for i in range(cell_info.GetNumberOfArrays())]


def _set_threshold_range(thresh, low, high):
    """Apply ThresholdBetween in a version-agnostic way."""
    try:
        # ParaView 5.10+
        thresh.ThresholdMethod = "Between"
        thresh.LowerThreshold = low
        thresh.UpperThreshold = high
    except AttributeError:
        # ParaView < 5.10
        thresh.ThresholdRange = [low, high]


def _create_layers(source, array_name, layer_defs, render_view):
    """Create one Threshold filter per entry in *layer_defs* and show each."""
    colors = _distinct_colors(len(layer_defs))
    created = []
    for i, (name, low, high, desc) in enumerate(layer_defs):
        thresh = Threshold(registrationName=name, Input=source)
        thresh.Scalars = ["CELLS", array_name]
        _set_threshold_range(thresh, low, high)

        display = Show(thresh, render_view)
        display.ColorArrayName = [None, ""]   # solid color, no scalar LUT
        display.DiffuseColor = colors[i]
        display.Opacity = 0.85

        created.append(thresh)
        print(f"  [OK] {name:50s}  ({desc})")

    return created


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    source = GetActiveSource()
    if source is None:
        raise RuntimeError(
            "No active source found.\n"
            "Please load a VTK file and select it in the Pipeline Browser first."
        )

    available = _cell_array_names(source)
    print(f"Available cell arrays in source: {available}")

    render_view = GetActiveViewOrCreate("RenderView")
    any_created = False

    if "tagnumber" in available:
        print("\n── tagnumber layers ──────────────────────────────────────────────")
        _create_layers(source, "tagnumber", TAGNUMBER_LAYERS, render_view)
        any_created = True
    else:
        print("Array 'tagnumber' not found – skipping tagnumber layers.")

    if "mediatype" in available:
        print("\n── mediatype layers ──────────────────────────────────────────────")
        _create_layers(source, "mediatype", MEDIATYPE_LAYERS, render_view)
        any_created = True
    else:
        print("Array 'mediatype' not found – skipping mediatype layers.")

    if not any_created:
        print(
            "\nNeither 'tagnumber' nor 'mediatype' arrays were found.\n"
            "Make sure you have selected the correct VTK source."
        )
        return

    # Hide the raw source so only the labelled layers are visible
    Hide(source, render_view)

    Render()
    ResetCamera()
    print("\nDone. All layers are visible in the Pipeline Browser and the 3-D view.")


main()
