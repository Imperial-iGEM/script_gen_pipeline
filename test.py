# Testing exception warning

import sys
# Yes this is awful but it lets modules from sibling directories be imported  https://docs.python.org/3/tutorial/modules.html#the-module-search-path
sys.path.insert(0,'../') # print('sys.path', sys.path)

from typing import List, Tuple, Union
from designs.construct import Construct, Variant
from uuid import uuid4

from script_gen_pipeline.protocol.mix_copy import Mix
from script_gen_pipeline.labware.containers_copy import Content


# Make toy construct
real_construct: Construct = Construct([[1,2,3,4,5], 'a','b','c','d', [1,2,3,4,5],'e'])

for mod_idx, module in enumerate(real_construct.modules):
    for part_idx, part in enumerate(module.parts):
        if mod_idx % 2 == 0:
            real_construct.modules[mod_idx].parts[part_idx].set_role('Linker')

real_construct = real_construct.update_construct()
clip_components: List[List[Variant]] = [[]]

constructs: List[Tuple[Variant]] = real_construct.get_unique_constructs()
for i, construct in enumerate(constructs):
    print(i, 'construct', construct)
construct = constructs[0]
# for mod_idx, module in enumerate(real_construct.modules):
#     for part in module.parts:
#         if not part.is_linker():
#             variants = []
#             for variant in part.variants:
#                 if not variant.is_linker():  # double check
#                     prefix = real_construct.modules[variant.prefix]
#                     suffix = real_construct.modules[variant.suffix]
#                     variants.append(variant)
#             clip_components.append([prefix, variant, suffix])
#         else:
#             print('part',part)
# clip_components.pop(0)

#####################
# 27.07.20
def get_variant(variant_list: List[Variant], linker_id):
    for var_idx, variant in enumerate(variant_list):
        if linker_id == (variant.prefix):
            return variant_list[var_idx-1]
        if linker_id == (variant.suffix):
            return variant_list[var_idx+1]
    return 0

print('construct',construct)

for variant in construct:
    if not variant.is_linker():
        prefix = get_variant(construct, variant.prefix)
        suffix = get_variant(construct, variant.suffix)
        clip_components.append([prefix, variant, suffix])
clip_components.pop(0)
print('clip_components', clip_components)
######################


######################
"""Reagents."""


class Reagent:
    """A Reagent. Ex T4 Ligase, Buffer, etc.

    Keyword Args:
        name: the reagent's name (default: {""})
    """

    def __init__(self, name: str = ""):
        assert name

        self.name = name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other) -> bool:
        return hash(self) == hash(other)


"""Lab species and organisms."""


class Species:
    """A Species. Ex E coli.

    Keyword Args:
        name: the species's name (default: {""})
    """

    def __init__(self, name: str = ""):
        assert name

        self.name = name

    def __hash__(self):
        """Hash species."""

        return hash(self.name)

    def __eq__(self, other) -> bool:
        return hash(self) == hash(other)


class Container:
    """A container with contents.

    TODO: make contents a set for uniqueness

    Keyword Args:
        contents: the contents of the container (default: {None})
        volumes: volumes of each content (default: {None})
    """

    volume_dead = -1
    """Volume that should be left unused at bottom of container."""

    volume_max = -1
    """Max volume within each container."""

    rows = 1
    """Rows in the container (or its parents, as with Wells and their Plate)"""

    cols = 1
    """Cols in the container (or its parents, as with Wells and their Plate)"""

    def __init__(
        self,
        contents: Union[Content, List[Content]] = None,
        volumes: List[float] = None,
    ):
        self.id = uuid4()

        if not contents:
            self.contents: List[Content] = []
        elif not isinstance(contents, list):
            self.contents = [contents]
        else:
            self.contents = contents

        self.volumes = volumes if volumes else [-1] * len(self.contents)
        self.withdrawn = 0.0  # volume with withdraws during pipette sim


class Well(Container):
    """A single well in a plate."""

    volume_max = 200
    volume_dead = 15
    rows = 8
    cols = 12


GOLDEN_GATE_MIX = Mix(
    {Reagent("master mix"): 4.0, Variant: 3.0},
    fill_with=Reagent("water"),
    fill_to=20.0,
)

mix = GOLDEN_GATE_MIX

# TEST
all_wells: List[List[Well]] = [[]]
for con_i, construct in enumerate(constructs):
    for variant in construct:
        if not variant.is_linker():
            prefix = get_variant(construct, variant.prefix)
            suffix = get_variant(construct, variant.suffix)
            clip_components.append([prefix, variant, suffix])
    clip_components.pop(0)
    print(con_i, 'clip_components', clip_components)

    wells = []
    for clip_part in clip_components:
        # add reaction mix and water
        # TODO: define Clip reaction mix for each well
        well_contents, well_volumes = mix(clip_part)

        # create a well that mixes the assembly mix, plasmids, and reagents
        well = Well(contents=well_contents, volumes=well_volumes)

        wells.append(well)
    all_wells.append(wells)
all_wells.pop(0)

# TEST 
n = [(0, 'e'), (3, 'e'), 'a', 'e', 't', 's', 'n']
for k in n:
    e = n.pop(-1)
    print('e', e)
    print('k', k)

shape = (8, 12)
layout_wells: List[List[Container]] = [[None] * shape[1]] * shape[0]
row_count = 0
col_count = 0
remaining_wells = all_wells
for i in range(len(remaining_wells)):
    added_well = remaining_wells[i].pop(0)
    if not layout_wells[row_count][col_count]:
        layout_wells[row_count][col_count] = added_well
        col_count += 1
        if col_count >= shape[1]:
            col_count = 0
            row_count += 1
            if row_count >= shape[0]:
                break
print('remaining_wells:', len(remaining_wells), remaining_wells)          
print(layout_wells)
