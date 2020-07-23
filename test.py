# Testing exception warning

import pandas as pd 
from typing import List
from designs.construct import Construct, Variant


f = {'prefixes': [1,1,1,1,1], 'parts': [2,2,2,2,2], 'suffixes': [3,3,3,3,3]}

real_construct = Construct(['a','b','c','d','e'])
# for i in range(5):
#     construct.append(Variant(i))

for mod_idx, module in enumerate(real_construct.modules):
    for part_idx, part in enumerate(module.parts):
        if mod_idx % 2 == 0:
            real_construct.modules[mod_idx].parts[part_idx].role = 'Linker'
        
real_construct = real_construct.update_construct()
print(len(real_construct.modules))

clip_components: List[List[Variant]] = [[]]

for mod_idx, module in enumerate(real_construct.modules):
    for part in module.parts:
        if not part.is_linker():
            variants = []
            for variant in part.variants:
                if not variant.is_linker():  # double check
                    prefix = real_construct.modules[variant.prefix]
                    suffix = real_construct.modules[variant.suffix]
                    variants.append(variant)
            clip_components.append([prefix, variant, suffix])
        else:
            print(part)
clip_components.pop(0)

print(clip_components)


