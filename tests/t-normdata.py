import os
_BE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys
sys.path.append(_BE)

from bioelectronica.normdata.bacstr import norm_bacstr

print(norm_bacstr())
