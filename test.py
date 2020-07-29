# Testing exception warning

import sys
# Yes this is awful but it lets modules from sibling directories be imported  https://docs.python.org/3/tutorial/modules.html#the-module-search-path
sys.path.insert(0,'../') # print('sys.path', sys.path)

from typing import List, Tuple, Union
from uuid import uuid4

from script_gen_pipeline.designs.construct import Construct, Variant
from script_gen_pipeline.labware.containers import Content, Container, Well
from script_gen_pipeline.labware.mix import Mix
from script_gen_pipeline.protocol.steps import Setup, Transfer
from script_gen_pipeline.protocol.protocol import Basic, Plate
from script_gen_pipeline.protocol.biochem_utils import Reagent, Species


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

for variant in construct:
    if not variant.is_linker():
        prefix = get_variant(construct, variant.prefix)
        suffix = get_variant(construct, variant.suffix)
        clip_components.append([prefix, variant, suffix])
clip_components.pop(0)
print('clip_components', clip_components)
######################


######################

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
    # print(con_i, 'clip_components', clip_components)

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

# TEST pop() in a loop
# n = [[0, 'e'], (3, 'e'), 'a', 'e', 't', 's', 'n']
# for k in n:
#     e = n.pop(-1)
#     print('e', e)
#     print('k', k)
# a = 1
# c = a or (n[0][0] if isinstance(n[0], list) else n[0])
# print('c', c)

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


def recursive_len(item):
    if type(item) == list:
        return sum(recursive_len(subitem) for subitem in item)
    else:
        return 1


print('remaining constructs:', len(remaining_wells))
print('remaining wells:', recursive_len(remaining_wells))

plate_count = 0
parameters = None
num_slots = parameters if parameters else 12  # OT default slots
slots: List[Container] = [[None]]*num_slots
slots[0] = Plate()
slots[1] = Plate()
slots[2] = Plate()

plate_count += sum(isinstance(p, Plate) for p in slots)
print('plate_count', plate_count)

# Test that Setup works
protocol = Basic(real_construct)

