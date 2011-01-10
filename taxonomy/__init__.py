try:
    __version__ = '0.1.'+'$Rev: 10 $'.split()[1]
except IndexError:
    __version__ = '0.1'

__version_info__ = tuple([ int(num) for num in __version__.split('.')])

import ncbi
import utils
import package
from alignment import Alignment
from taxonomy import Taxonomy

