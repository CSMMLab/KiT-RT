import subprocess, os

project = 'KiT-RT'
copyright = '2022, Karlsruhe Institute of Technology'
author = 'Jonas Kusch, Steffen Schotthöfer, Pia Stammer, Jannick Wolters, Tianbai Xiao'

version = 'v0.1'
release = 'v0.1'


extensions = [
    'breathe',
    'sphinx_rtd_theme'
]

#templates_path = ['_templates']
master_doc = 'index'
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'display_version': False
}
html_logo = 'images/KiT-RT_logo_small.png'
html_static_path = ['_static']
def setup(app):
    app.add_css_file('theme_overrides.css')

breathe_projects = {
    "KiT-RT": "../build/docs/doxygen/xml",
}

breathe_default_project = "KiT-RT"

read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
if read_the_docs_build:
    inputDir = '../'
    outputDir = '../build/docs/doxygen'
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    with open("Doxyfile.in", "rt") as fin:
        with open("Doxyfile", "wt") as fout:
            for line in fin:
                line = line.replace('GENERATE_HTML          = YES', 'GENERATE_HTML          = NO')
                line = line.replace('SEARCH_INCLUDES        = YES', 'SEARCH_INCLUDES        = NO')
                line = line.replace('CLASS_DIAGRAMS         = YES', 'CLASS_DIAGRAMS         = NO')
                line = line.replace('@DOXYGEN_OUTPUT_DIR@', outputDir)
                line = line.replace('@DOXYGEN_INPUT_DIR@', inputDir)
                fout.write(line)
    subprocess.call('doxygen', shell=True)
