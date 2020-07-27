# Testing exception warning

from typing import Any, Dict, Iterable, List, Set, Tuple, Union
from designs.construct import Construct, Variant
from uuid import uuid4
import inspect
import logging

# Make toy construct
real_construct: Construct = Construct([[1,2,3,4,5], 'a','b','c','d', [1,2,3,4,5],'e'])

for mod_idx, module in enumerate(real_construct.modules):
    for part_idx, part in enumerate(module.parts):
        if mod_idx % 2 == 0:
            real_construct.modules[mod_idx].parts[part_idx].set_role('Linker')

real_construct = real_construct.update_construct()
print('len(real_construct.modules)', len(real_construct.modules))
print('real_construct.modules[0].parts[0].role', real_construct.modules[0].parts[0].role)

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

#####################
# 27.07.20
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


Content = Union[Construct, Reagent]


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


class Mix:

    def __init__(
        self,
        mix: Dict[Any, float] = None,
        fill_with: Content = None,
        fill_to: float = 0,
    ):
        if fill_with and not fill_to:
            raise ValueError(
                f"Cannot specify what to 'fill_with' without specifying container's 'fill_to'"
            )

        if not fill_with and fill_to:
            raise ValueError(
                f"Cannot specify what to 'fill_to' without 'fill_with' Content (Reagent, Species or SeqRecord)"
            )

        self.mix = mix or {}
        self.fill_with = fill_with
        self.fill_to = fill_to
    
    def __call__(
        self, contents: Iterable[Content]
    ) -> Tuple[List[Content], List[float]]:

        contents_out: List[Content] = []
        volumes: List[float] = []

        seen: Set[Any] = set()

        for content in contents:
            if (
                isinstance(content, Reagent) or isinstance(content, Species)
            ) and content in self.mix:
                # reagent was explicitly specified
                contents_out.append(content)
                volumes.append(self.mix[content])
                seen.add(content)
            elif type(content) in self.mix:
                # example is SeqRecord
                contents_out.append(content)
                volumes.append(self.mix[type(content)])
                seen.add(type(content))
            elif any(
                isinstance(content, t) for t in self.mix.keys() if inspect.isclass(t)
            ):
                # example is RestrictionType class
                class_type = next(
                    t
                    for t in self.mix.keys()
                    if inspect.isclass(t) and isinstance(content, t)
                )
                contents_out.append(content)
                volumes.append(self.mix[class_type])
                seen.add(class_type)
            else:
                logging.warning(f"Content {content} not found in mix")

        for content, volume in self.mix.items():
            if content in seen or inspect.isclass(content):
                continue

            contents_out.append(content)
            volumes.append(volume)

        if self.fill_to and self.fill_with:
            volume_total = sum(volumes)
            volume_remaining = max([self.fill_to - volume_total, 0.0])
            contents_out.append(self.fill_with)
            volumes.append(volume_remaining)

        return (contents_out, volumes)


GOLDEN_GATE_MIX = Mix(
    {Reagent("master mix"): 4.0, Variant: 2.0},
    fill_with=Reagent("water"),
    fill_to=20.0,
)

mix = GOLDEN_GATE_MIX

# TEST
clip_component = clip_components[0]
well_contents, well_volumes = mix(clip_component)

# create a well that mixes the assembly mix, plasmids, and reagents
well = Well(contents=well_contents, volumes=well_volumes)

print('well_contents', well_contents)
print('well_volumes', well_volumes)
print('well.volumes', well.volumes)
