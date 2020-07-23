# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 11:18:14
#  * @modify date 2020-07-17 11:18:14
#  * @desc [description]
#  */

""" Parsing of SBOL input happens through Construct class.
Hierarchy of parts here is Construct > Module > Part > Variant
Expect construct to come with all necessary parts (eg. prefix / 
suffix linkers, scars)
"""

from typing import List
from warnings import warn
import numpy as np
from uuid import uuid4


class Variant:
    """ Equivalent to a part on SynBioHub or Parts Registry.
    The lowest level in combinatorial derivation.
    Attributes:
        annotations: Equivalent to SBOL annotations + scars
        id: Unique id for internal referencing, SBOL unrelated
        module_id: Unique id of the module this variant is in
        name: Actual part name (eg. BBa_K10002)
        prefix: BASIC linker that comes before a part (Module) in construct
        role: same as Part.role, equivalent to SBOL roles + custom linker
        sequence: Equivalent to SBOL sequence
        suffix: BASIC linker that comes after a part (Module) in construct
        uri: Equivalent to SBOL URI
    """
    def __init__(self, component):
        self.component = component

        self.id = uuid4()
        self.name = self.get_name()
        self.module_id = None
        self.module_order_idx = None

        self.uri = self.get_uri()
        self.sequence = self.get_seq()
        self.annotations = self.get_annotations()
        self.role = None

        self.prefix = None
        self.suffix = None
        
    def get_name(self):
        warn("NotImplem: get SBOL part name")
        # name = pysbol.get_part(self.component)
        return f'BBA_fake_{self.component}'  # str(self.id)[0:3]

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

    def is_linker(self):
        return (self.role == 'Linker')

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name


class Part:
    """ Equivalent to an SBOL root component. Store part 
    information as in a ComponentDefinition but with more
    explicit combinatorial variation handling. All parts
    consist of Variants.

    Attributes:
        comb_variants: all combinatorial variations on this part
        role: same as SBOL role (promoter, DNA, RBS...), but
            includes linkers (and assembly specific parts)
        module_id: the unique id of a Module that the Part is in
        range: sequence information (relation)
    """
    def __init__(self, component):
        print("Parse component into constituent Variant")
        self.role = self.get_role(component)
        self.variants: List(Variant) = self.make_variants(component)
        self.module_id = None  # Set once Modules are made
        self.id = uuid4()

        print("NotImplem: define roles ('Linker') through ids not strs")
        if self.role != "Linker":
            self.prefix = None  # self.module_id
            self.suffix = None

        # self.range = self.get_range()

    def get_role(self, component):
        warn("NotImplem: get component role SBOL style")
        role = 'Yuh'
        return role

    def is_linker(self):
        """ Check if this part's role is 'Linker' """
        warn("make sure is_linker() matches the same fxn in Variant")
        return (self.role == 'Linker')

    def get_module_id(self):
        pass

    def make_variants(self, component):
        variants: List(Variant) = []
        comb_ders = self.unpack_comb_ders(component)

        for single_comb_def in comb_ders:
            variant = Variant(single_comb_def)
            variant.role = self.role
            variants.append(variant)
        return variants

    def unpack_comb_ders(self, component):
        """ Enumerate each combinatorial design in this part; 
        refer to Ming's combinatorial derivation code"""
        comb_ders = component
        warn(f"NotImplem: from the root component {component} get child variant components")
        return [comb_ders]

    def set_module_info(self, module_id, **kwargs):
        """ Set the module id of the current part and its variants """
        self.module_id = module_id
        # propagate to variants
        temp_variants: List[Variant] = []
        for variant in self.variants:
            variant.module_id = module_id
            for key, value in kwargs.items():
                variant.key = value
                # variant.module_order_idx = module_order_idx
            temp_variants.append(variant)
        self.variants = temp_variants

    def __len__(self):
        print(f"Using the length of Part {self.role}: {self.id}")
        return 1

    def __repr__(self):
        """ Print the object representation as its role and id """
        return f"{self.role}: {self.id}"

    


class Module():
    """ A Module is a unit of assembly. Way of grouping parts """
    def __init__(self, order_idx, parts):
        self.id = uuid4()
        self.order_idx = order_idx          # use integers that reflect module order
        self.parts: List[Part] = self.make_parts_list(parts)
        self.name = f'Module {self.order_idx}'

    def make_parts_list(self, parts):
        """ Make parts input list type and propagate module info to parts """
        parts = parts if isinstance(parts, list) else [parts]

        for part in parts:
            part.set_module_info(module_id=self.id, 
                module_order_idx=self.order_idx)
        return parts


class Construct():
    """ Highest level of design. A construct is the full plasmid
    including all parts and variations. Takes in the output of the 
    SBOL Designer editing / processing step, possibly JSON (TBD)

    Attributes:
        id: unique internal id for Construct, unrelated to SBOL
        modules: List of modules making up construct, way of grouping
            parts and translating them into buildeable units. Ordered
            by index of creation starting at 0.
        parts: List of all parts in construct without Module hierarchy
        unique_constructs: construct hierarchy of modules > parts > variants
            is flattened to the Variant level. Each entry in list corresponds
            to one unique construct.
    """
    def __init__(self, sbol_input):

        self.id = uuid4()
        # self.sbol_input = sbol_input
        print("\n\nshould be making modules now")
        self.modules: List[Module] = self.make_modules(sbol_input)
        """ Might be able to make modules right away from sbol_input
        depending on output of SBOL Designer """
        self.parts: List[Part] = self.make_parts()
        self.simp_modules: List[str] = self._simplify_modules()
        """ Modules in order of construct assembly """
        self.unique_constructs: List[List[Variant]] = None        

    def make_modules(self, sbol_input):
        """ Return all the parts within the final construct
        that are building blocks for assembly """

        parts: List[Part] = []
        components = self.get_components(sbol_input)

        for component in components:
            part = Part(component)
            parts.append(part)

        modules: List[Module] = []

        for i, part in enumerate(parts):
            module = Module(i, part)
            modules.append(module)
        
        return modules

    def get_components(self, sbol_input):
        warn("NotImplem: get_components() should return the root ComponentDefs of an SBOL input")
        sbol_input = sbol_input if isinstance(sbol_input, list) else [sbol_input]
        return sbol_input

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

    def make_parts(self) -> List[Part]:
        """ Create list of parts """
        parts: List[Part] = []
        for module in self.modules:
            for part in module.parts:
                parts.append(part)
        return parts

    def update_parts(self):
        self.parts: List[Part] = []
        for module in self.modules:
            self.parts.append(module.parts)

    def fuse_modules(self, module_1: Module, module_2: Module):
        assert self.check_adjacent(module_1.order_idx, module_2.order_idx), f"Please pick adjacent modules to fuse instead of {module_1.order_idx} and {module_2.order_idx}"

        order_idx = min(module_1.order_idx, module_2.order_idx)
        return Module(order_idx, [module_1.parts, module_2.parts])

    def check_adjacent(self, num1, num2):
        num1 = int(np.floor(num1))
        num2 = int(np.floor(num2))
        is_adjacent = num1 < (num2+1) and num1 > (num2-1)
        return is_adjacent

    def check_module_order(self):
        """ Check that the module order makes sense """
        last_order_idx = self.modules[0].order_idx
        for module in self.modules:
            if last_order_idx != module.order_idx:
                raise "The modules are not listed in the correct order."
            last_order_idx = module.order_idx

    def get_unique_constructs(self, 
                remove_modules: List = None) -> List[List[Variant]]:
        """ List each unique, full construct possible by flattening the
        construct by Variants in order of assembly. Optionally use the 
        remove_modules argument to remove specific modules.
        """
        unique_constructs: List[List[Variant]] = None
        print("NotImplem: loop the modules list, within that loop each \
            part, within that each Variant. Maybe set self.unique_constructs\
            too?")
        return unique_constructs

    def _set_pref_suff(self):
        """ Set the prefix and suffix of each variant as the module id.
        Propagate the module id of linker prefix and suffixes to 
        the parts they are flanking """
        print('\n\n\nset_pref_suff: type(self.modules[0]) \n\n\n', type(self.modules[0]))

        for mod_idx, module in enumerate(self.modules):  # keep things diagonal :)  
            for part_idx, part in enumerate(module.parts):
                for var_idx, variant in enumerate(part.variants):
                    if variant.is_linker():
                        continue
                    else:
                        if mod_idx != 0:
                            if self.modules[mod_idx-1].parts[0].is_linker():
                                variant.prefix = self.modules[mod_idx-1].order_idx
                                print(variant.prefix)
                            if mod_idx != (len(self.modules)-1):
                                if self.modules[mod_idx+1].parts[0].is_linker():
                                    variant.suffix = self.modules[mod_idx+1].order_idx
                                    print(variant.suffix)
                    # Update
                    self.modules[mod_idx].parts[part_idx].variants[var_idx] = variant

        print('variant.prefix',variant.prefix, 'variant.suffix', variant.suffix)
        return self

    def _simplify_modules(self):
        """ Returns list of module names """
        simplified_modules: List[str] = []
        for module in self.modules:
            simplified_modules.append(module.name)
        return simplified_modules

    def update_construct(self):
        """ Perform various validations and propagate updated 
        values to child lists (eg new Module info) """

        self._set_pref_suff()
        return self
