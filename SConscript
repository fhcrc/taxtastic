import os

def add_src(p):
    return(os.path.join("/home/matsen/wopi/repo/src/mokaphy/", p))

env = Environment(ENV=os.environ)
inkscape = Builder(action = 'inkscape --without-gui --export-pdf=$TARGET $SOURCE')
env['BUILDERS']['Inkscape'] = inkscape
env['BUILDERS']['Cp'] = Builder(action = 'cp $SOURCE $TARGET')

env.Cp(target="library.bib", source="/home/matsen/Mendeley/library.bib")
final=env.PDF(target='algotax.pdf',source='algotax.tex')
Depends(final, "algotax.bib")

