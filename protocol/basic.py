"""BASIC assembly design process and steps."""

from typing import List

from Bio.Restriction.Restriction import RestrictionType
from Bio.Restriction import BsaI
from Bio.SeqRecord import SeqRecord

from synbio.containers import Container, Well
from synbio.designs import Design
from synbio.instructions import Temperature
from synbio.mix import Mix
from synbio.reagents import Reagent
from synbio.steps import Setup, Pipette, ThermoCycle, HeatShock
from synbio.protocol import Protocol

#BASIC_MIX = Mix()

class Basic(Protocol):
    """ BASIC assembly

        See: https://pubs.acs.org/doi/pdf/10.1021/sb500356d

        Takes in BASIC design containing parts, linkers, and backbone to produce a
        set of instructions. Checks to ensure existence of compatible linkers between
        each part and a backbone. Sites are cut using BsaI restriction enzyme. 

        Inspired by DNABot. 
    """
    
    def __init__(self, 
        design: Design = Design(),
        name: str = "",
        enzymes: List[RestrictionType] = [BsaI],
        include: List[str] = None,
        separate_reagents: bool = False,
        #mix = BASIC_MIX, 
    ):
        super().__init__(name=name, design=design, separate_reagents=separate_reagents)
        #self.mix = mix
        self.enzymes = enzymes

    def run(self):

    def _create_mixed_wells(self):


