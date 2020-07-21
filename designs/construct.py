# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 11:18:14
#  * @modify date 2020-07-17 11:18:14
#  * @desc [description]
#  */

""" Parsing of SBOL input happens through Construct class. 
Hierarchy of parts here is Construct > Module > Part > Variant == BuildIn 
"""

from typing import List
from warnings import warn
import numpy as np
from uuid import uuid4


class Construct():
    """ Highest level of design. A construct is the full plasmid
    including all parts and variations. Takes in the output of the 
    SBOL Designer editing / processing step, possibly JSON (TBD)

    Attributes:
        id: unique internal id for Construct, unrelated to SBOL
        parts: List of all parts in construct without Module hierarchy
        modules: List of modules making up construct, way of grouping
            parts and translating them into buildeable units
    """
    def __init__(self, sbol_input):

        self.id = uuid4()
        self.parts: List[Part] = self.make_parts(sbol_input)
        self.modules: List[Module] = self.make_modules()
        self.simp_modules: List[str] = self.simplify_modules()
        """ Modules in order of construct assembly """

    # The following couple of fxns may be unnecessary if parsing happens in java
    def add_part(self, part: Part, rel_location):
        """ Add a part to the construct. Automatically make module
        from part and update flat list of parts.
        Args:
            part: New Part to be added
            rel_location: Model order index int of new part location
        """
        self.modules.insert(rel_location, Module(part))
        self.update_parts()
        return self

    def rm_part(self, bad_part: Part):
        """ Remove a part from the construct """
        for part in self.parts:
            if bad_part.id == part.id:
                self.parts.remove(part)
        return self

    def make_parts(self, sbol_input):
        parts: List[Part] = []
        components = self.get_components(sbol_input)

        for component in components:
            part = Part(component)
            parts.append(part)
        return parts

    def update_parts(self):
        self.parts: List[Part] = []
        for module in self.modules:
            self.parts.append(module.parts)

    def get_components(self, sbol_input):
        warn("NotImplem: get_components() should return the root ComponentDefs of an SBOL input")
        return 0

    def make_modules(self):
        """ Return all the parts within the final construct
        that are building blocks for assembly """
        modules: List[Module] = []

        for i, part in enumerate(self.parts):
            module = Module(i, part)
            modules.append(module)
        return modules

    def fuse_modules(self, module_1: Module, module_2: Module):
        assert self.check_adjacent(module_1.order_idx, module_2.order_idx), f"Please pick adjacent modules to fuse instead of {module_1.order_idx} and {module_2.order_idx}"

        order_idx = min(module_1.order_idx, module_2.order_idx)
        return Module(order_idx, [module_1.parts, module_2.parts])

    def check_adjacent(self, num1, num2):
        num1 = int(np.floor(num1))
        num2 = int(np.floor(num2))
        is_adjacent = num1 < (num2+1) and num1 > (num2-1)
        return is_adjacent

    def simplify_modules(self):
        """ Returns list of module names """
        simplified_modules: List[str] = ''
        for module in self.modules:
            simplified_modules.append(module.name)
        return simplified_modules


class Module():
    """ A Module is a unit of assembly. Way of grouping parts """
    def __init__(self, order_idx, parts):
        self.parts: List[Part] = parts if len(parts) > 1 else [parts]
        self.id = uuid4() 
        self.order_idx = order_idx          # use integers that reflect module order
        self.name = f'Module {self.order_idx}'


class Part:
    """ Equivalent to an SBOL root component. Store part 
    information as in a ComponentDefinition but with more
    explicit combinatorial variation handling. All parts
    consist of Variants.

    Attributes:
        comb_variants: all combinatorial variations on this part
        role: same as SBOL role (promoter, DNA, RBS...) but 
            includes linkers and assembly specific parts 
        module_id: the unique id of a Module that the Part is in
        range: sequence information (relation)
    """
    def __init__(self, component):
        print("Parse component into constituent Variant")
        self.comb_variants: List(Variant) = self.make_variants(component)
        self.role = self.get_role(component)
        self.module_id = 0  # Set once Modules are made
        self.id = uuid4()

        self.buildins: List[BuildIn] = []
        self.prefix: BuildIn = 0
        self.suffix: BuildIn = 0

        self.range = self.get_range()

    def get_module_id(self):
        pass

    def make_variants(self, component):
        variants: List(Variant) = []
        comb_ders = self.unpack_comb_ders(component)

        for single_comb_def in comb_ders:
            variant = Variant(single_comb_def)
            variants.append(variant)
        return variants

    def unpack_comb_ders(self, component):
        """ Enumerate each combinatorial design in this part; 
        refer to Ming's combinatorial derivation code"""
        comb_ders = 0
        warn(f"NotImplem: from the root component {component} get child variant components")
        return comb_ders

    def get_role(self, component):
        warn("NotImplem: get component role SBOL style")
        role = ''
        return role

    def __len__(self):
        return 1


class Variant:
    """ Equivalent to a part on SynBioHub or Parts Registry.
    The lowest level in combinatorial derivation.
    Attributes:
        name: Actual part name (eg. BBa_K10002)
        id: Unique id for internal referencing, SBOL unrelated
        uri: Equivalent to SBOL URI
        sequence: Equivalent to SBOL sequence
        annotattion: Equivalent to SBOL annotations + scars
    """
    def __init__(self, component):
        self.component = component
        self.name = self.get_name()
        self.id = uuid4()
        self.uri = self.get_uri()
        self.sequence = self.get_seq()
        self.annotations = self.get_annotations()
        self.subpart = "DO WE WANT THESE"

    def get_name(self):
        warn("NotImplem: get SBOL part name")
        # name = pysbol.get_part(self.component)
        return 'BBA_fake'

    def get_uri(self):
        warn("NotImplem: get SBOL uri")
        # uri = pysbol.get_uri(self.component)
        return 0

    def get_seq(self):
        warn("NotImplem: get SBOL DNA sequence. Watch out for internal references.")
        return 0

    def get_annotations(self):
        warn("NotImplem: get SBOL annotations")
        return 0


class BuildIn(Variant):
    """ Edits that need to integrate with part, like iP+iS.
    Might be unnecessary.
    """

    def __init__(self):
        self.overlaps = ''
        self.part_ref: float = 0
        # """Reference to relative position within the part"""