project = "AIM"
author = "Chung Research Group"
copyright = "2025–2026, Muhammad Hassan and Chung Research Group"

extensions = [
    "sphinx.ext.mathjax",
    "sphinx_design",
]
templates_path = ["_templates"]
exclude_patterns = []

html_theme = "pydata_sphinx_theme"
html_title = "AIM documentation"
html_static_path = ["_static"]
html_css_files = ["aim.css"]
html_logo = "_static/AIM_logo.png"
html_favicon = "_static/AIM_logo.png"
html_theme_options = {
    "show_toc_level": 2,
    "navigation_with_keys": True,
    "navbar_align": "left",
    "header_links_before_dropdown": 6,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/Chung-Research-Group/AIM",
            "icon": "fa-brands fa-github",
        }
    ],
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
}

master_doc = "index"
