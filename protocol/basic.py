"""BASIC assembly design process and steps."""

from collections import defaultdict
from typing import Dict, List, Set, Tuple, Iterable, Optional
import pandas as pd

#from Bio.Restriction.Restriction import RestrictionType
#from Bio.Restriction import BsaI
#from Bio.SeqRecord import SeqRecord

#from synbio.containers import Container, Well
#from synbio.designs import Design
#from synbio.instructions import Temperature
from synbio.reagents import Reagent
#from synbio.steps import Setup, Pipette, ThermoCycle, HeatShock
#from synbio.protocol import Protocol

from script_gen_pipeline.protocol.instructions import Instruction, instr_to_txt, Temperature
from script_gen_pipeline.labware.containers import Container, Fridge, Well, Plate
from script_gen_pipeline.labware.mix import Mix
from script_gen_pipeline.designs.construct import Construct, Module, Part
from script_gen_pipeline.protocol.protocol import Protocol, Step, Subprotocol

#BASIC_MIX = Mix()

# Constant floats/ints - from DNABot - move to parameters?
CLIP_DEAD_VOL = 60
CLIP_VOL = 30
T4_BUFF_VOL = 3
BSAI_VOL = 1
T4_LIG_VOL = 0.5
CLIP_MAST_WATER = 15.5
PART_PER_CLIP = 200
MIN_VOL = 1
MAX_CONSTRUCTS = 96
MAX_CLIPS = 48
FINAL_ASSEMBLIES_PER_CLIP = 15
DEFAULT_PART_VOL = 1
MAX_SOURCE_PLATES = 6

CLIP_OUT_PATH = '1_clip.ot2.py'
MAGBEAD_OUT_PATH = '2_purification.ot2.py'
F_ASSEMBLY_OUT_PATH = '3_assembly.ot2.py'
TRANS_SPOT_OUT_PATH = '4_transformation.ot2.py'
basic_steps = [CLIP_OUT_PATH, MAGBEAD_OUT_PATH, F_ASSEMBLY_OUT_PATH, TRANS_SPOT_OUT_PATH]

basic_mix = Mix(
    {Reagent("Promega T4 DNA Ligase buffer, 10X"): T4_BUFF_VOL, 
    Reagent("NEB BsaI-HFv2"): BSAI_VOL, 
    Reagent("Promega T4 DNA Ligase"): T4_LIG_VOL, 
    Part: DEFAULT_PART_VOL, Module: DEFAULT_PART_VOL, 
    Part: DEFAULT_PART_VOL}, 
    fill_with=Reagent("water"), fill_to=CLIP_VOL,
)

class Basic(Protocol):
    """ BASIC assembly

        See: https://pubs.acs.org/doi/pdf/10.1021/sb500356d

        Takes in BASIC construct to produce a set of instructions. 
        Checks to ensure existence of compatible linkers between
        each part and a backbone. Sites are cut using BsaI restriction enzyme. 

        Inspired by synbio and DNABot. 
    """
    
    def __init__(self, 
        construct: Construct = Construct(),
        name: str = "",
        #source_wells: Dict[str] = [], 
    ):
        super().__init__(name=name, construct=construct)
        self.mix = basic_mix
        #self.source_wells = source_wells
        self.parameters = {
            'SPOTTING_VOLS_DICT': {2: 5, 3: 5, 4: 5, 5: 5, 6: 5, 7: 5},
            'SOURCE_DECK_POS': ['2', '5', '8', '7', '10', '11'],
            'ethanol_well_for_stage_2': "A11"
        }
        self.scripts = [CLIP_OUT_PATH, MAGBEAD_OUT_PATH, F_ASSEMBLY_OUT_PATH, TRANS_SPOT_OUT_PATH]
        self.subprotocols = [Subprotocol(str(script), script, self.construct, self.parameters) for script in self.scripts]

    def run(self):
        self.clip_df = self._create_clip_df()

        # clip reaction subprotocol
        # purification subprotocol
        # assembly subprotocol
        # transformation subprotocol
        for subprotocol in self.subprotocols:
            self = subprotocol(self)
            self.history.append(subprotocol)
            self.generate_ot_script(self, assay, template_script)
            raise NotImplementedError

    def _create_clip_df(self):
        clips_info = {'prefixes': [], 'parts': [],
	              'suffixes': []}
        for index, module in enumerate(self.construct.modules):
            if index % 2 != 0:
                clips_info['parts'].append(module)
                prefix_linker = self.construct.modules[index - 1].parts[0]
                clips_info['prefixes'].append(prefix_linker)
                if index == len(self.construct.modules) - 1:
                    suffix_linker = self.construct.modules[0].parts[0]
                    clips_info['suffixes'].append(suffix_linker)
                else:
                    suffix_linker = self.construct.modules[index + 1].parts[0]
                    clips_info['suffixes'].append(suffix_linker)
        clips_df = pd.DataFrame.from_dict(clips_info)
        # add a 'number' column for future?
        return clips_df
    
    def _create_source_plate(self):
        # list of modules -> plate
        source_plate = Plate()
        return source_plate

    def _create_mixed_wells(self):
        # clips_df -> Plate, use basic_mix
        mixed_wells = Plate()
        for clip_info in self.clip_df.iterrows():
            modules = []
            modules.append(clip_info['prefixes'])
            modules.append(clip_info['parts'])
            modules.append(clip_info['suffixes'])
            # need to make mix compatible
            #well_contents, well_volumes = self.mix(modules)
            #mixed_wells.add_wells(well_contents, well_volumes)
        return mixed_wells



            


