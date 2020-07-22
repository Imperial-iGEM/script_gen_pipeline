"""BASIC assembly design process and steps."""

from collections import defaultdict
from typing import Dict, List, Set, Tuple, Iterable, Optional
import pandas as pd

from Bio.Restriction.Restriction import RestrictionType
from Bio.Restriction import BsaI
from Bio.SeqRecord import SeqRecord

#from synbio.containers import Container, Well
#from synbio.designs import Design
#from synbio.instructions import Temperature
#from synbio.mix import Mix
#from synbio.reagents import Reagent
#from synbio.steps import Setup, Pipette, ThermoCycle, HeatShock
#from synbio.protocol import Protocol

from script_gen_pipeline.protocol.instructions import Instruction, instr_to_txt, Temperature
from script_gen_pipeline.labware.containers_copy import Container, Fridge, Well
from script_gen_pipeline.designs.construct import Construct, Module, Part
from script_gen_pipeline.protocol.protocol import Protocol, Step

#BASIC_MIX = Mix()

SOURCE_DECK_POS = ['2', '5', '8', '7', '10', '11'] # example pos, used in DNABot

class Basic(Protocol):
    """ BASIC assembly

        See: https://pubs.acs.org/doi/pdf/10.1021/sb500356d

        Takes in BASIC construct containing parts, linkers, and backbone to produce a
        set of instructions. Checks to ensure existence of compatible linkers between
        each part and a backbone. Sites are cut using BsaI restriction enzyme. 

        Inspired by DNABot. 
    """
    
    def __init__(self, 
        construct: Construct = Construct(),
        name: str = "",
        source_wells: Dict[str, str] = []
    ):
        super().__init__(name=name, construct=construct)
        #self.mix = mix
        self.source_wells = source_wells
        # source well keys in form ['well', 'deck/plate']
        # source well values in form [Part, Part, ...] or [Part]

    def run(self):
        self.clip_df = self._create_clip_df()
        self.source_df = self._check_source_wells()

    def _create_clip_df(self):
        clips_info = {'prefixes': [], 'parts': [],
	              'suffixes': []}
        for index, module in enumerate(self.construct.modules):
            if index % 2 != 0:
                clips_info['parts'].append(module.parts)
                prefix_linker = self.construct.modules[index - 1].parts[0]
                clips_info['prefixes'].append(prefix_linker)
                if index == len(self.construct.modules) - 1:
                    suffix_linker = self.construct.modules[0].parts[0]
                    clips_info['suffixes'].append(suffix_linker)
                else:
                    suffix_linker = self.construct.modules[index + 1].parts[0]
                    clips_info['suffixes'].append(suffix_linker)
        clips_df = pd.DataFrame.from_dict(clips_info)
        return clips_df
        
    def _check_source_wells(self):
        if not source_wells:
            # generate well positions for parts, include warning
        else:
            # check to ensure all parts included in dict values
            # gen well positions for parts if any missing wells, include warning

        # return dataframe for easier lookup



            


