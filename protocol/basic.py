"""BASIC assembly design process and steps."""

from collections import defaultdict
from typing import Dict, List, Set, Tuple, Iterable, Optional
import pandas as pd
import numpy as np

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
SOURCE_VOL = 15 # dead vol of 10-15 uL recommended for each part/linker

CLIP_OUT_PATH = '1_clip.ot2.py'
MAGBEAD_OUT_PATH = '2_purification.ot2.py'
F_ASSEMBLY_OUT_PATH = '3_assembly.ot2.py'
TRANS_SPOT_OUT_PATH = '4_transformation.ot2.py'
basic_steps = [CLIP_OUT_PATH, MAGBEAD_OUT_PATH, F_ASSEMBLY_OUT_PATH, TRANS_SPOT_OUT_PATH]

basic_mix = Mix(
    {Reagent("Promega T4 DNA Ligase buffer, 10X"): T4_BUFF_VOL, 
    Reagent("NEB BsaI-HFv2"): BSAI_VOL, 
    Reagent("Promega T4 DNA Ligase"): T4_LIG_VOL, 
    Module: DEFAULT_PART_VOL, Module: DEFAULT_PART_VOL, 
    Module: DEFAULT_PART_VOL}, 
    fill_with=Reagent("water"), fill_to=CLIP_VOL,
)

source_mix = Mix({Module: SOURCE_VOL}) # not sure if this is right vol

class Basic(Protocol):
    """ BASIC assembly

        See: https://pubs.acs.org/doi/pdf/10.1021/sb500356d

        Takes in BASIC construct to produce a set of instructions. 
        Checks to ensure existence of compatible linkers between
        each part and a backbone. Sites are cut using BsaI restriction enzyme. 

        Inspired by synbio and DNABot. 
    """
    
    def __init__(self, 
        constructs: List[Construct] = [Construct()],
        name: str = "",
        #source_wells: Dict[str] = [], 
    ):
        super().__init__(name=name, constructs=constructs)
        self.mix = basic_mix
        #self.source_wells = source_wells
        self.parameters = {
            'SPOTTING_VOLS_DICT': {2: 5, 3: 5, 4: 5, 5: 5, 6: 5, 7: 5},
            'SOURCE_DECK_POS': ['2', '5', '8', '7', '10', '11'],
            'ethanol_well_for_stage_2': "A11"
        }
        self.scripts = [CLIP_OUT_PATH, MAGBEAD_OUT_PATH, F_ASSEMBLY_OUT_PATH, TRANS_SPOT_OUT_PATH]
        self.subprotocols = [Subprotocol(str(script), script, self.constructs, self.parameters) for script in self.scripts]

    def run(self):
        self.clip_df, self.master_mix = self._create_clip_df()
        self.source_plate, self.source_info = self._create_source_plate()
        self.mixed_wells = self._create_mixed_wells()

        # clip reaction subprotocol
        # purification subprotocol
        # assembly subprotocol
        # transformation subprotocol
        for subprotocol in self.subprotocols:
            self = subprotocol(self)
            self.history.append(subprotocol)
            self.generate_ot_script(self, assay, template_script)
            raise NotImplementedError

    def _get_construct_modules(self, construct):
        clips_info = {'prefixes': [], 'parts': [],
	              'suffixes': []}
        for index, module in enumerate(construct.modules):
            if index % 2 != 0:
                clips_info['parts'].append(module)
                prefix_linker = construct.modules[index - 1]
                clips_info['prefixes'].append(prefix_linker)
                if index == len(construct.modules) - 1:
                    suffix_linker = construct.modules[0]
                    clips_info['suffixes'].append(suffix_linker)
                else:
                    suffix_linker = construct.modules[index + 1]
                    clips_info['suffixes'].append(suffix_linker)
        clips_info_df = pd.DataFrame.from_dict(clips_info)
        return clips_info_df
    
    def _get_final_well(self, sample_number):
        """Determines well containing the final sample from sample number.
        """
        letter = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        final_well_column = sample_number // 8 + \
            (1 if sample_number % 8 > 0 else 0)
        final_well_row = letter[sample_number - (final_well_column - 1) * 8 - 1]
        return final_well_row + str(final_well_column)

    def _create_clip_df(self):
        constructs_list = []
        for construct in self.constructs:
            constructs_list.append(self._get_construct_modules(construct))

        merged_construct_dfs = pd.concat(constructs_list, ignore_index=True)
        unique_clips_df = merged_construct_dfs.drop_duplicates()
        unique_clips_df = unique_clips_df.reset_index(drop=True)
        clips_df = unique_clips_df.copy()

        # Error
        if len(unique_clips_df.index) > MAX_CLIPS:
            raise ValueError(
                'Number of CLIP reactions exceeds 48.')

        # Count number of each CLIP reaction
        clip_count = np.zeros(len(clips_df.index))
        for i, unique_clip in unique_clips_df.iterrows():
            for _, clip in merged_construct_dfs.iterrows():
                if unique_clip.equals(clip):
                    clip_count[i] = clip_count[i] + 1
        clip_count = clip_count // FINAL_ASSEMBLIES_PER_CLIP + 1
        clips_df['number'] = [int(i) for i in clip_count.tolist()]

        # Associate well/s for each CLIP reaction
        clips_df['mag_well'] = pd.Series(['0'] * len(clips_df.index),
                                        index=clips_df.index)
        for index, number in clips_df['number'].iteritems():
            if index == 0:
                mag_wells = []
                for x in range(number):
                    mag_wells.append(self._get_final_well(x + 1 + 48))
                clips_df.at[index, 'mag_well'] = tuple(mag_wells)
            else:
                mag_wells = []
                for x in range(number):
                    well_count = clips_df.loc[
                        :index - 1, 'number'].sum() + x + 1 + 48
                    mag_wells.append(self._get_final_well(well_count))
                clips_df.at[index, 'mag_well'] = tuple(mag_wells)

        multiple = (clips_df['number'].sum())*CLIP_DEAD_VOL/CLIP_VOL
        # in future mutiple = (clips_df['number'].sum())*CLIP_DEAD_VOL/CLIP_VOL
        master_mix = Mix({Reagent("Promega T4 DNA Ligase buffer, 10X"): multiple*T4_BUFF_VOL, 
        Reagent("NEB BsaI-HFv2"): multiple*BSAI_VOL, 
        Reagent("Promega T4 DNA Ligase"): multiple*T4_LIG_VOL, 
        Reagent("water"): multiple*CLIP_MAST_WATER})
        return clips_df, master_mix
    
    def _create_source_plate(self):
        # list of modules -> plate
        # is this necessary?
        # what is the mix?
        source_plate = Plate()
        source_info = {'modules': [], 'well_index': [], 'well': []}
        modules_list = []
        for construct in self.constructs:
            for module in construct.modules:
                # need to create new df saying module with corresponding index
                if module not in modules_list:
                    modules_list.append(module)
                    well_contents, well_volumes = source_mix(module)
                    well = Well(well_contents, well_volumes)
                    indx = source_plate.add_wells(well)
                    source_info['modules'].append(module)
                    source_info['well_index'].append(indx)
                    source_info['well'].append(well)
        return source_plate, source_info

    def _create_mixed_wells(self):
        # clips_df -> Plate, use basic_mix
        mixed_wells = Plate()
        for clip_index, clip_info in self.clip_df.iterrows():
            modules = []
            modules.append(clip_info['prefixes'])
            modules.append(clip_info['parts'])
            modules.append(clip_info['suffixes'])
            well_contents, well_volumes = self.mix(modules)
            well = Well(well_contents, well_volumes)
            indx = mixed_wells.add_wells(well)
            self.clip_df.insert(clip_index, 'well', well)
            self.clip_df.insert(clip_index, 'well_index', indx)

        return mixed_wells



            


