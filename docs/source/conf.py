
# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------
project = 'Quarto Project'
copyright = '2024'
author = 'User'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx_book_theme'
]
templates_path = ['_templates']
exclude_patterns = ['_build','**/.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_book_theme'
html_theme_options = {
    "repository_url": "https://github.com/your-username/your-repo",
    "use_repository_button": True,
}
html_static_path = ['_static']
