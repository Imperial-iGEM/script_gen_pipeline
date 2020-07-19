# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 11:18:14
#  * @modify date 2020-07-17 11:18:14
#  * @desc [description]
#  */

from typing import List

class Construct():
    def __init__(self, sbol_input):

        self.id = ""
        self.parts: List[Part] = make_parts(sbol_input)
        """ List of all parts in construct """
        self.modules: List[Module] = []
        """ Modules in order of construct assembly """

        self.input_construct_path = input_construct_path
        self.output_sources_paths = output_sources_paths

        pass

    def add_part(self, part, rel_location):
        return self.parts.insert(rel_location, part)

    def make_parts(self, sbol_input):
        parts: List[Part] = []

        for component in sbol_input:
            part = Path(component)
            parts.append(part)
        return parts

    def get_all_modules(self):
        """ Return all the parts within the final construct
        that are building blocks for assembly """
        raise NotImplementedError


class Module():
    """ A Module is a unit of assemble """
    def __init__(self):
        self.parts: List[Part] = []
        self.id = uuid4() 
        self.name = get_module_name() # use integers that reflect module order

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
        self.comb_variants: List(Variant) = make_variants(component)
        self.role = get_role()
        self.module_id = 0  # Set once Modules are made
        self.buildins: List[BuildIn] = []

        self.range = get_range()

    def get_module_id(self):
        pass

    def make_variants(self, component):
        variants: List(Variant) = []
        comb_ders = unpack_comb_ders(component)

        for single_comb_def in component:
            variant = Variant(single_comb_def)
            variants.append(variant)
        return variants


class Variant:
    def __init__(self):
        self.name = get_name()
        self.id = make_id()
        self.uri = get_uri()
        self.sequence = get_seq()
        self.subpart = "DO WE WANT THESE"


class BuildIn(Variant):
    """ Edits that need to integrate with part, like iP+iS.
    Might be unnecessary.
    """

    def __init__(self):
        self.overlaps = ''
        self.part_ref: float = 0
        # """Reference to relative position within the part"""