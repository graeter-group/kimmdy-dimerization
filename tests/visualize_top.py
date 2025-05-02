from kimmdy.parsing import read_top
from kimmdy.topology.topology import Topology
from pathlib import Path
from kimmdy.tools import write_top_as_dot
import copy

path_a = read_top(Path('Ala_out.top'), use_gmx_dir=False)
top = Topology(path_a)
top_a = copy.deepcopy(top)
top_b = top

top_b.break_bond(('7', '8'))
top_b.bind_bond(('8', '9'))

write_top_as_dot(top_a, "ala-top-a.dot")
write_top_as_dot(top_b, "ala-top-b.dot")