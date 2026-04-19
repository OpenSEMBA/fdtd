from __future__ import annotations

project = "semba-fdtd"
author = "OpenSEMBA"
copyright = "OpenSEMBA"

extensions = [
    "myst_parser",
    "sphinx.ext.githubpages",
]

source_suffix = {
    ".md": "markdown",
}

root_doc = "index"
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "templates",
]

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
]
myst_heading_anchors = 3

html_theme = "sphinx_rtd_theme"
html_title = "semba-fdtd documentation"
html_static_path: list[str] = []
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
}
html_context = {
    "display_github": True,
    "github_user": "OpenSEMBA",
    "github_repo": "fdtd",
    "github_version": "main",
    "conf_py_path": "/doc/",
}
