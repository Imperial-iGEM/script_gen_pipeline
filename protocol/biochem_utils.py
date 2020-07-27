# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-27 11:58:04
#  * @modify date 2020-07-27 11:58:04
#  * @desc [description]
#  */
""" Contains classes supporting definitions of chemicals """
# TODO: Add Primers?


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


class Species:
    """Lab species and organisms. A Species. Ex E coli. 

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