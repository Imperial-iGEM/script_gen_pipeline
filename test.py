# Testing exception warning

import pandas as pd 
from typing import List, Tuple
from designs.construct import Construct, Variant


f = {'prefixes': [1,1,1,1,1], 'parts': [2,2,2,2,2], 'suffixes': [3,3,3,3,3]}

real_construct: Construct = Construct([[1,2,3,4,5], 'a','b','c','d', [1,2,3,4,5],'e'])
# for i in range(5):
#     construct.append(Variant(i))

for mod_idx, module in enumerate(real_construct.modules):
    for part_idx, part in enumerate(module.parts):
        if mod_idx % 2 == 0:
            real_construct.modules[mod_idx].parts[part_idx].set_role('Linker')
        
real_construct = real_construct.update_construct()
print('len(real_construct.modules)',len(real_construct.modules))
print('real_construct.modules[0].parts[0].role',real_construct.modules[0].parts[0].role)

clip_components: List[List[Variant]] = [[]]

constructs: List[Tuple[Variant]] = real_construct.get_unique_constructs()
for i, construct in enumerate(constructs):
    print(i, 'construct',construct)
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


def get_variant(variant_list: List[Variant], linker_id):
    print('[get_variant]')
    for var_idx, variant in enumerate(variant_list):
        print(variant, variant.module_order_idx, 'variant.prefix',variant.prefix)
        print(variant, variant.module_order_idx, 'variant.suffix',variant.suffix)
        # if linker_id == (variant.prefix or variant.suffix):
        if linker_id == (variant.prefix):
            print('PREF returned;',variant)
            return variant_list[var_idx-1]
        if linker_id == (variant.suffix):
            print('SUFF returned;',variant)
            return variant_list[var_idx+1]
    return 0


print('construct',construct)

for variant in construct:
    if not variant.is_linker():
        print('variant.role',variant.role)
        print('giving variant.prefix', variant.prefix)
        prefix = get_variant(construct, variant.prefix)
        print('giving variant.suffix', variant.suffix)
        suffix = get_variant(construct, variant.suffix)
        clip_components.append([prefix, variant, suffix])
clip_components.pop(0)
print('clip_components',clip_components)


