import os
_BE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys
sys.path.append(_BE)

from bioelectronica.normdata.bacstr import norm_bacstr

_txt = [
  "Sporosarcina psychrophila DSM 3",
  "Helicobacter pylori NCTC 11637 = CCUG  17874 = ATCC 43504 = JCM 12093",
  "Deinococus radiodurans R1 = ATCC 13939 = DSM 20539",
  "Deinococus radioduns R1 = ATCC 13939 = DSM 20539"
]

_input = "; ".join(_txt)

_re = norm_bacstr(_input)
print(_re)
