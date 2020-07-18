# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 14:10:20
#  * @modify date 2020-07-17 14:10:20
#  * @desc [description]
#  */

from typing import List, Iterable
import csv
import pandas as pd
import numpy as np
import os

import Instruction
from script_gen_pipeline.protocol.instructions import instr_to_txt
import Container
import Fridge 

import Step
import Construct


class Protocol:
    def __init__(self, construct: Construct = Construct()):
        self.construct = construct  # the final construct to be built
        self.steps: List[Step] = []  # list of steps for this assembly

        # each of the below is set during a 'run()'
        self.containers: List[Container] = []
        self.instructions: List[Instruction] = []  # assembly instructions

        raise NotImplementedError

    def generate_ot_script(self, assay, template_script):
        """ Build a python script for Opentrons. Based on DNAbot

        Args:
            assay: name of the type of protocol to be run on liquid handler
            template_script: pathname to the python script used as template
        """

        raise NotImplementedError

    def add_step(self, step: Step) -> "Protocol":
        """Add an instruction step to the protocol for documentation.

        Args:
            step: The Step to add to this protocol
        """

        if not isinstance(step, Step):
            raise TypeError

        self.steps.append(step)

        return self

    def run(self) -> "Protocol":
        """ Core running of a protocol """

        # all input records, each start out assigned to a single Fridge source
        records = self.construct.get_all_modules()

        # update containers
        self.containers = [Fridge(r) for r in records]  # everything comes from fridge

        for step in self.steps:
            step(self)

        return self

    # TODO: this was originally decorated with @property
    def accum_content(self) -> List[Construct]:
        """ Based on 'output' def in Lattice synbio Protocol code. 
        Gather the output contents from the final containers after protocol
        has been run. 

        Returns: 
            Constructs or collection of parts that were created from the assembly
        """

        accumulator: List[Construct] = []
        for container in self.containers:
            for content in container:     # custom iter in container (can also try [c for c in gibson_well if isinstance(c, SeqRecord)])
                if isinstance(content, Construct):
                    accumulator.append(content)
        return accumulator

    # TODO: unpack the DNA sequence for each Construct component
    def to_fasta(self, filename: str = "") -> int:
        """Write each output record to a FASTA file.

        Uses `SeqIO.write(records, filename, "fasta")`.

        Args:
            filename: The filename to write the FASTA file to

        Returns:
            The number of records that were written
        """

        if not filename:
            filename = self._filename() + ".fasta"

        self._check_output()
        return SeqIO.write(self.output, filename, "fasta")

    # TODO: unpack the DNA sequence for each Construct component
    def to_genbank(self, filename: str = "", split: bool = False) -> int:
        """Write each output record to a Genbank file.

        Uses `SeqIO.write(records, filename, "genbank")`.

        Args:
            filename: The filename to write the Genbanks to
            split: Write a separate Genbank for each SeqRecord

        Returns:
            The number of records that were written
        """

        self._check_output()

        # limit for id is 16 characters https://github.com/biopython/biopython/issues/747
        # shorten ids of those that exceed the limit
        output_renamed: List[SeqRecord] = []
        for i, output in enumerate(self.output):
            if len(output.id) > 16:
                output = output.upper()
                output.description = output.id
                output.id = "Seq" + str(i + 1)
                output_renamed.append(output)
            else:
                output_renamed.append(output)

        if split:
            write_dir = os.path.dirname(filename)
            for record, renamed_record in zip(self.output, output_renamed):
                filename = os.path.join(write_dir, record.id + ".gb")
                SeqIO.write(renamed_record, filename, "genbank")
            return len(output_renamed)

        if not filename:
            filename = self._filename() + ".gb"

        return SeqIO.write(output_renamed, filename, "genbank")

    def to_txt(self, filename: str = ""):
        """Write the protocol's instructions to a text file.

        Args:
            filename: the filename of the instructin file
        """

        if not filename:
            filename = self._filename() + ".csv"

        protocol_txt = instr_to_txt(self.name, self.instructions)
        with open(filename, "w") as instruction_file:
            instruction_file.write(protocol_txt)

    # TODO: check this 
    def to_csv(self, filename: str = "") -> str:
        """Write CSV file(s) describing the containers/Layout after each step.

        Rows of Layout/containers are written in CSV format with each step's name as its heading.

        Args:
            filename: The name of the CSV's name
        """

        if not filename:
            filename = self._filename() + ".csv"

        csv = ""
        row = 0
        for instruction in self.instructions:
            if not instruction.transfers:
                continue

            name = instruction.name
            if not name and instruction.instructions:
                name = instruction.instructions[0]
            row += 1
            csv += f"{name}:\n" if name else f"Setup step {row}:\n"
            csv += Layout.from_instruction(
                instruction,
                existing_plates=self.instruction_to_plate_count[instruction],
                log_volume=row == 1,
                separate_reagents=self.separate_reagents,
            ).to_csv()

        with open(filename, "w") as csvfile:
            csvfile.write(csv)

        return csv

    def to_picklists(self, filename: str = "", platform: str = "tecan"):
        """Create picklists for robotic pipetting.

        Supported platforms are `tecan`, `hamilton`, and `labcyte`.

        For each step where there's plate to plate pipetting, create a
        robotic picklists. Steps where reagents or samples come from the Fridge
        are not written to a picklist right now.

        If there are multiple passable steps, each are saved with their
        index in their filename. Ex: picklist.1.gwl, picklist.2.gwl

        Keyword Args:
            filename: Name of picklist file (default: {self.name})
            platform: Picklist platform (default: {"tecan"})
        """

        picklist_generators = {
            "tecan": to_tecan,
            "hamilton": to_hamilton,
            "labcyte": to_labcyte,
        }
        if platform not in picklist_generators:
            picklist_platforms = ", ".join(picklist_generators.keys())
            raise ValueError(
                f"'{platform}' is an unrecognized platform. Choose from: {picklist_platforms}"
            )

        self._check_output()

        if not filename:
            filename = self._filename()

            if platform == "tecan":
                filename += ".gwl"
            elif platform == "labcyte":
                filename += ".csv"

        # accumulate instructions from the protocol that are from plate to plate
        picklist_instructions: List[Instruction] = []
        for instruction in self.instructions:
            if not instruction.transfers:
                continue

            srcs = {t.src for t in instruction.transfers}
            dests = {t.dest for t in instruction.transfers}

            def no_fridge(containers: Iterable[Container]) -> bool:
                return all(not isinstance(s, Fridge) for s in containers)

            if no_fridge(srcs) and no_fridge(dests):
                picklist_instructions.append(instruction)

        if not picklist_instructions:
            raise RuntimeWarning(f"no picklist-capable steps in protocol")

        def picklist_filename(index: int) -> str:
            if len(picklist_instructions) == 1:
                return filename

            fname, fext = os.path.splitext(filename)
            return fname + str(index + 1) + fext

        for i, instruction in enumerate(picklist_instructions):
            picklist = picklist_generators[platform](
                instruction, self.instruction_to_plate_count[instruction]
            )

            with open(picklist_filename(i), "w") as picklist_file:
                picklist_file.write(picklist)

    def add_instruction(self, instruction: Instruction):
        """Add a single Instruction to this protocol's output.

        Instructions are generated by Steps. Each Step calls this to add the
        step's output to the Protocol for accumulation.

        Args:
            instruction: A single instruction to add to the protocol
        """

        self.instructions.append(instruction)
        self.instruction_to_plate_count[instruction] = self.plate_count

        # add any pippete-able Layout
        if instruction.transfers:
            dest_containers = {t.dest for t in instruction.transfers}
            if all(not isinstance(c, Fridge) for c in dest_containers):
                self.plate_count += len(
                    Layout.from_instruction(
                        instruction, separate_reagents=self.separate_reagents
                    )
                )

    def _check_output(self):
        """Verify that the Protocol has steps and that they have been run.

        If there are no steps, or if there is no list of output Records
        in self.containers, there is nothing to write to a FASTA file.
        """
        if not self.containers:
            self.run()

            if not self.output:
                raise RuntimeError(
                    "Failed to create assemblies after executing all steps."
                )

    # TODO: adjust filename creation
    def _filename(self) -> str:
        """Gather a filename from this protocol's name.

        Convert this protocol's name to a name that's safe for a file.
        The below is from a Github Gist:
            https://gist.github.com/wassname/1393c4a57cfcbf03641dbc31886123b8
        """

        filename = self.name
        whitelist = "-_.() %s%s" % (string.ascii_letters, string.digits)
        filename = filename.replace(" ", "_")

        # keep only valid ascii chars
        cleaned_filename = (
            unicodedata.normalize("NFKD", filename).encode("ASCII", "ignore").decode()
        )

        # keep only whitelisted chars
        return "".join(c for c in cleaned_filename if c in whitelist)

    # Might as well leave the built-ins
    def __str__(self) -> str:
        """Return a human readable summary of the protocol.

        The python built in function str calls this. Should output a summary like:

        ```txt
        Combinatorial MoClo
            how: combinatorial
            design: [[15 x SeqRecord] [10 x SeqRecord] [2 x SeqRecord]]
            steps: []
        ```
        """

        name = self.name
        how = f"\thow: {self.design.__name__}"
        design = f"\tdesign: {str(self.design)}"

        return "\n".join([name, how, design])

    def __iter__(self) -> Iterable[Step]:
        """Iterate over the steps in this protocol."""

        return iter(self.steps)

    def __len__(self) -> int:
        """Return the number of steps in this protocol.

        Returns:
            the number of steps
        """

        return len(self.steps)


class Clone(Protocol):
    def __init__(self, design):
        self.design = design
        super().__init__()


# Taken from DNAbot app
CLIP_OUT_PATH = '1_clip.ot2.py'
MAGBEAD_OUT_PATH = '2_purification.ot2.py'
F_ASSEMBLY_OUT_PATH = '3_assembly.ot2.py'
TRANS_SPOT_OUT_PATH = '4_transformation.ot2.py'
basic_steps = [CLIP_OUT_PATH, MAGBEAD_OUT_PATH, F_ASSEMBLY_OUT_PATH, TRANS_SPOT_OUT_PATH]


class Basic(Clone):
    """ Specific protocol run of BASIC Assembly as described
    in DNAbot
    """

    def __init__(self):
        super().__init__()

        self.steps = [CLIP_OUT_PATH, MAGBEAD_OUT_PATH, F_ASSEMBLY_OUT_PATH, TRANS_SPOT_OUT_PATH]
        inputs = [input_construct_path, output_sources_paths]
        subprotocols = [Subprotocol(str(step), step, inputs) for step in self.steps]        

    def run(self) -> "Protocol":

        for step in subprotocols:
            generate_ot_script(self, assay, template_script)
            raise NotImplementedError

    def _generate_clips_dict(self):
        pass


class Subprotocol(Protocol):
    def __init__(self, name, out_pathname, inputs):
        super().__init__()

        self.name = name
        self.out_pathname = out_pathname
        self.args = kwargs

        print('Processing input csv files...')
        constructs_list = self._generate_constructs_list(inputs[0])
        clips_df = self.generate_clips_df(constructs_list)
        print(output_sources_paths)
        sources_dict = self.generate_sources_dict(output_sources_paths)

        # calculate OT2 script variables
        print('Calculating OT-2 variables...')
        clips_dict = self.generate_clips_dict(clips_df, sources_dict)
        magbead_sample_number = clips_df['number'].sum()
        final_assembly_dict = self.generate_final_assembly_dict(constructs_list,
                                                        clips_df)
        final_assembly_tipracks = self.calculate_final_assembly_tipracks(
            final_assembly_dict)
        spotting_tuples = self.generate_spotting_tuples(constructs_list,
                                                SPOTTING_VOLS_DICT)


        if self.name == str(basic_steps[0]):
            clips_dict=clips_dict
        if self.name == str(basic_steps[1]):
            sample_number=magbead_sample_number
            # THIS IS FOR THE ETHANOL TROUGH WELL IN STEP 2
            ethanol_well=ethanol_well_for_stage_2
        if self.name == str(basic_steps[2]):
            final_assembly_dict=final_assembly_dict
            tiprack_num=final_assembly_tipracks
        if self.name == str(basic_steps[3]):
            spotting_tuples=spotting_tuples
            #Deep well plate for Soc media during
            #soc_well="A{}".format(dnabotinst.soc_column))
            #previously the information about the location of the
            soc_well="A1"z


    def generate_constructs_list(path):
        """Generates a list of dataframes corresponding to each construct. Each 
        dataframe lists components of the CLIP reactions required.
        """

        def process_construct(construct):
            """Processes an individual construct into a dataframe of CLIP reactions
            outlining prefix linkers, parts and suffix linkers.
            """

            def interogate_linker(linker):
                """Interogates linker to determine if the suffix linker is a UTR
                linker.
                """
                if len(linker) >= 4:
                    if linker[:3] == 'UTR':
                        return linker[:4] + '-S'
                else:
                    return linker + "-S"

            clips_info = {'prefixes': [], 'parts': [],
                        'suffixes': []}
            for i, sequence in enumerate(construct):
                if i % 2 != 0:
                    clips_info['parts'].append(sequence)
                    clips_info['prefixes'].append(
                        construct[i - 1] + '-P')
                    if i == len(construct) - 1:
                        suffix_linker = interogate_linker(construct[0])
                        clips_info['suffixes'].append(suffix_linker)
                    else:
                        suffix_linker = interogate_linker(construct[i + 1])
                        clips_info['suffixes'].append(suffix_linker)
            return pd.DataFrame.from_dict(clips_info)

        constructs_list = []
        with open(path, 'r') as csvfile:
            csv_reader = csv.reader(csvfile)
            for index, construct in enumerate(csv_reader):
                if index != 0:  # Checks if row is header.
                    construct = list(filter(None, construct))
                    if not construct[1:]:
                        break
                    else:
                        constructs_list.append(process_construct(construct[1:]))

        # Errors
        if len(constructs_list) > MAX_CONSTRUCTS:
            raise ValueError(
                'Number of constructs exceeds maximum. Reduce construct number in construct.csv.')
        else:
            return constructs_list
        
    def generate_clips_df(constructs_list):
        """Generates a dataframe containing information about all the unique CLIP 
        reactions required to synthesise the constructs in constructs_list.
        """
        merged_construct_dfs = pd.concat(constructs_list, ignore_index=True)
        unique_clips_df = merged_construct_dfs.drop_duplicates()
        unique_clips_df = unique_clips_df.reset_index(drop=True)
        clips_df = unique_clips_df.copy()

        # Error
        if len(unique_clips_df.index) > MAX_CLIPS:
            raise ValueError(
                'Number of CLIP reactions exceeds 48. Reduce number of constructs in construct.csv.')

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
                    mag_wells.append(final_well(x + 1 + 48))
                clips_df.at[index, 'mag_well'] = tuple(mag_wells)
            else:
                mag_wells = []
                for x in range(number):
                    well_count = clips_df.loc[
                        :index - 1, 'number'].sum() + x + 1 + 48
                    mag_wells.append(final_well(well_count))
                clips_df.at[index, 'mag_well'] = tuple(mag_wells)
        return clips_df


    def generate_sources_dict(paths):
        """Imports csvs files containing a series of parts/linkers with 
        corresponding information into a dictionary where the key corresponds with
        part/linker and the value contains a tuple of corresponding information.
        Args:
            paths (list): list of strings each corresponding to a path for a 
                        sources csv file. 
        """
        sources_dict = {}
        for deck_index, path in enumerate(paths):
            print(path)
            with open(path, 'r') as csvfile:
                csv_reader = csv.reader(csvfile)
                for index, source in enumerate(csv_reader):
                    if index != 0:
                        csv_values = source[1:]
                        csv_values.append(SOURCE_DECK_POS[deck_index])
                        sources_dict[str(source[0])] = tuple(csv_values)
        return sources_dict


    def generate_clips_dict(clips_df, sources_dict):
        """Using clips_df and sources_dict, returns a clips_dict which acts as the 
        sole variable for the opentrons script "clip.ot2.py".
        """
        max_part_vol = CLIP_VOL - (T4_BUFF_VOL + BSAI_VOL + T4_LIG_VOL
                                + CLIP_MAST_WATER + 2)
        clips_dict = {'prefixes_wells': [], 'prefixes_plates': [],
                    'suffixes_wells': [], 'suffixes_plates': [],
                    'parts_wells': [], 'parts_plates': [], 'parts_vols': [],
                    'water_vols': []}

        # Generate clips_dict from args
        try:
            for _, clip_info in clips_df.iterrows():
                prefix_linker = clip_info['prefixes']
                clips_dict['prefixes_wells'].append([sources_dict[prefix_linker][0]]
                                                    * clip_info['number'])
                clips_dict['prefixes_plates'].append(
                    [handle_2_columns(sources_dict[prefix_linker])[2]] * clip_info['number'])
                suffix_linker = clip_info['suffixes']
                clips_dict['suffixes_wells'].append([sources_dict[suffix_linker][0]]
                                                    * clip_info['number'])
                clips_dict['suffixes_plates'].append(
                    [handle_2_columns(sources_dict[suffix_linker])[2]] * clip_info['number'])
                part = clip_info['parts']
                clips_dict['parts_wells'].append([sources_dict[part][0]]
                                                * clip_info['number'])
                clips_dict['parts_plates'].append([handle_2_columns(sources_dict[part])[2]]
                                                * clip_info['number'])
                if not sources_dict[part][1]:
                    clips_dict['parts_vols'].append([DEFAULT_PART_VOL] *
                                                    clip_info['number'])
                    clips_dict['water_vols'].append([max_part_vol - DEFAULT_PART_VOL]
                                                    * clip_info['number'])
                else:
                    part_vol = round(
                        PART_PER_CLIP / float(sources_dict[part][1]), 1)
                    if part_vol < MIN_VOL:
                        part_vol = MIN_VOL
                    elif part_vol > max_part_vol:
                        part_vol = max_part_vol
                    water_vol = max_part_vol - part_vol
                    clips_dict['parts_vols'].append(
                        [part_vol] * clip_info['number'])
                    clips_dict['water_vols'].append(
                        [water_vol] * clip_info['number'])
        except KeyError:
            sys.exit('likely part/linker not listed in sources.csv')
        for key, value in clips_dict.items():
            clips_dict[key] = [item for sublist in value for item in sublist]
        return clips_dict


    def generate_final_assembly_dict(constructs_list, clips_df):
        """Using constructs_list and clips_df, returns a dictionary of final
        assemblies with keys defining destination plate well positions and values
        indicating which clip reaction wells are used.
        """
        final_assembly_dict = {}
        clips_count = np.zeros(len(clips_df.index))
        for construct_index, construct_df in enumerate(constructs_list):
            construct_well_list = []
            for _, clip in construct_df.iterrows():
                clip_info = clips_df[(clips_df['prefixes'] == clip['prefixes']) &
                                    (clips_df['parts'] == clip['parts']) &
                                    (clips_df['suffixes'] == clip['suffixes'])]
                clip_wells = clip_info.at[clip_info.index[0], 'mag_well']
                clip_num = int(clip_info.index[0])
                clip_well = clip_wells[int(clips_count[clip_num] //
                                        FINAL_ASSEMBLIES_PER_CLIP)]
                clips_count[clip_num] = clips_count[clip_num] + 1
                construct_well_list.append(clip_well)
            final_assembly_dict[final_well(
                construct_index + 1)] = construct_well_list
        return final_assembly_dict


    def calculate_final_assembly_tipracks(final_assembly_dict):
        """Calculates the number of final assembly tipracks required ensuring
        no more than MAX_FINAL_ASSEMBLY_TIPRACKS are used.
        """
        final_assembly_lens = []
        for values in final_assembly_dict.values():
            final_assembly_lens.append(len(values))
        master_mix_tips = len(list(set(final_assembly_lens)))
        total_tips = master_mix_tips + sum(final_assembly_lens)
        final_assembly_tipracks = total_tips // 96 + (
            1 if total_tips % 96 > 0 else 0)
        if final_assembly_tipracks > MAX_FINAL_ASSEMBLY_TIPRACKS:
            raise ValueError(
                'Final assembly tiprack number exceeds number of slots. Reduce number of constructs in constructs.csv')
        else:
            return final_assembly_tipracks


    def generate_spotting_tuples(constructs_list, spotting_vols_dict):
        """Using constructs_list, generates a spotting tuple
        (Refer to 'transformation_spotting_template.py') for every column of 
        constructs, assuming the 1st construct is located in well A1 and wells
        increase linearly. Target wells locations are equivalent to construct well
        locations and spotting volumes are defined by spotting_vols_dict.
        Args:
            spotting_vols_dict (dict): Part number defined by keys, spottting
                volumes defined by corresponding value.
        """
        # Calculate wells and volumes
        wells = [final_well(x + 1) for x in range(len(constructs_list))]
        vols = [SPOTTING_VOLS_DICT[len(construct_df.index)]
                for construct_df in constructs_list]

        # Package spotting tuples
        spotting_tuple_num = len(constructs_list)//8 + (1
                                                        if len(constructs_list) % 8 > 0 else 0)
        spotting_tuples = []
        for x in range(spotting_tuple_num):
            if x == spotting_tuple_num - 1:
                tuple_wells = tuple(wells[8*x:])
                tuple_vols = tuple(vols[8*x:])
            else:
                tuple_wells = tuple(wells[8*x:8*x + 8])
                tuple_vols = tuple(vols[8*x:8*x + 8])
            spotting_tuples.append((tuple_wells, tuple_wells, tuple_vols))
        return spotting_tuples


    def generate_ot2_script(ot2_script_path, template_path, **kwargs):
        """Generates an ot2 script named 'ot2_script_path', where kwargs are 
        written as global variables at the top of the script. For each kwarg, the 
        keyword defines the variable name while the value defines the name of the 
        variable. The remainder of template file is subsequently written below.        
        """
        print("output location of ot2_script_path:{}".format(ot2_script_path))
        print(os.path.realpath(ot2_script_path))
        this_object_output_path = os.path.realpath(ot2_script_path)

        #current_path = os.getcwd()
        #remove_example = os.path.split("my_examples")
        #writing_path = os.path.join(current_path, "output")
        #output_path = os.path.join(writing_path, ot2_script_path)
        
        with open(ot2_script_path, 'w') as wf:
            with open(template_path, 'r') as rf:
                for index, line in enumerate(rf):
                    if line[:3] == 'def':
                        function_start = index
                        break
                    else:
                        wf.write(line)
                for key, value in kwargs.items():
                    wf.write('{}='.format(key))
                    if type(value) == dict:
                        wf.write(json.dumps(value))
                    elif type(value) == str:
                        wf.write("'{}'".format(value))
                    else:
                        wf.write(str(value))
                    wf.write('\n')
                wf.write('\n')
            with open(template_path, 'r') as rf:
                for index, line in enumerate(rf):
                    if index >= function_start - 1:
                        wf.write(line)
        return this_object_output_path





