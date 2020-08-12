# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-22 22:25:11
#  * @modify date 2020-07-22 22:25:11
#  * @desc [description]
#  */


""" originally inside Clip_Reaction class. Gets non-linker variants 
list from list of variants """
def _get_parts(construct: List[Variant]) -> List[List[Variant]]:
    """ Return list of non-linker Variants, keep module hierarchy through
    nesting, so you get parts = [[M1_part, M1_part], [M2_part]] """
    parts: List[List[Variant]] = [[]]
    same_module_parts = []
    # Build up nested list
    for variant in construct:
        if not variant.is_linker():
            if not same_module_parts:  # first variant
                same_module_parts.append(variant)
            else:
                if variant.module_id == part[-1].module_id:
                    same_module_parts.append(variant)
                else:
                    parts.append(part)
                    part = []
        else:
            continue
    return(parts)





""" OG Subprotocol class before Gabby updated the template scripts
and Clip class was created """


class Subprotocol(Protocol):
    """ A subprotocol is equivalent to one run on a liquid handler
    without human input necessary. The actual steps of the protocol are
    created here. """
    def __init__(self, name, out_pathname, construct, parameters):
        self.name = name
        self.out_pathname = out_pathname

        self.input_construct_path = construct.input_construct_path
        self.output_sources_paths = construct.output_sources_paths

        self.parameters = parameters
        self.kwargs = self._get_basic_kwargs()
        # super().__init__()

        self.template_script = self._get_template()

    def _get_basic_kwargs(self):
        """ Get the keyword arguments unique to each 
        BASIC subprotocol (clip, purification...) needed
        for the creation of opentrons scripts. Completely based
        on DNAbot's opentrons script generator reqs and inputs.

        Returns: 
            Keyword arguments unique to the name of the subprotocol
        """

        print('Processing input csv files...')
        constructs_list = self.generate_constructs_list(
            self.input_construct_path)
        clips_df = self.generate_clips_df(constructs_list)
        sources_dict = self.generate_sources_dict(self.output_sources_paths)

        # calculate OT2 script variables
        print('Calculating OT-2 variables...')
        clips_dict = self.generate_clips_dict(clips_df, sources_dict)
        magbead_sample_number = clips_df['number'].sum()
        final_assembly_dict = self.generate_final_assembly_dict(
            constructs_list, clips_df)
        final_assembly_tipracks = self.calculate_final_assembly_tipracks(
            final_assembly_dict)
        spotting_tuples = self.generate_spotting_tuples(constructs_list,
                                                        self.parameters['SPOTTING_VOLS_DICT'])

        # TODO: make human clear which basic step corresponds to each case
        if self.name == str(basic_steps[0]):
            return {'clips_dict': clips_dict}
        if self.name == str(basic_steps[1]):
            return {'sample_number': magbead_sample_number,
                    # THIS IS FOR THE ETHANOL TROUGH WELL IN STEP 2
                    'ethanol_well': self.parameters['ethanol_well_for_stage_2']}
        if self.name == str(basic_steps[2]):
            return {'final_assembly_dict': final_assembly_dict,
                    'tiprack_num': final_assembly_tipracks}
        if self.name == str(basic_steps[3]):
            return {'spotting_tuples': spotting_tuples,
                    # Deep well plate for Soc media during
                    # soc_well="A{}".format(dnabotinst.soc_column))
                    # previously the information about the location of the
                    'soc_well': "A1"}
        else:
            assert self.name in basic_steps, f"The subprotocol {self.name} is not one of the possible protocols {basic_steps}"
            return 0

    def _get_template(self):
        if self.name == str(basic_steps[0]):
            return 'assembly_template.py'
        if self.name == str(basic_steps[1]):
            return 'clip_template.py'
        if self.name == str(basic_steps[2]):
            return 'purification_template.py'
        if self.name == str(basic_steps[3]):
            return 'transformation_template.py'
        else:
            assert self.name in basic_steps, f"Template script = '', subprotocol name {self.name} does not match any of {basic_steps}"
            return ''

    def generate_constructs_list(self, path):
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
        
    def generate_clips_df(self, constructs_list):
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

    def __str__(self):
        return self.name

    def generate_ot2_script(ot2_script_path, template_path, **kwargs):
        """Generates an ot2 script named 'ot2_script_path', where kwargs are 
        written as global variables at the top of the script. For each kwarg, the 
        keyword defines the variable name while the value defines the name of the 
        variable. The remainder of template file is subsequently written below.        
        """
        print("output location of ot2_script_path:{}".format(ot2_script_path))
        print(os.path.realpath(ot2_script_path))
        this_object_output_path = os.path.realpath(ot2_script_path)

        # current_path = os.getcwd()
        # remove_example = os.path.split("my_examples")
        # writing_path = os.path.join(current_path, "output")
        # output_path = os.path.join(writing_path, ot2_script_path)
        
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

