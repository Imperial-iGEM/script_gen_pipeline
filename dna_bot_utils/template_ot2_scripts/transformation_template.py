from opentrons import protocol_api
from opentrons import legacy_api
import numpy as np

metadata = {'apiLevel': '2.2',
            'protocolName': 'Transformation Template v2',
            'author': 'Gabrielle Johnston',
            'description': 'DNABot updated transformation template'}


def run(protocol: protocol_api.ProtocolContext):
    def generate_transformation_wells(spotting_tuples):
        """Evaluates spotting_tuples and returns transformation wells.

        Args:
        spotting_tuples (list): Sets of spotting reactions are given 
            in the form: ((source wells), (target wells), (spotting volumes)). 
            Each unique transformation well is resuspended once prior to spotting.

        """
        wells = []
        for spotting_tuple in spotting_tuples:
            for source_well in spotting_tuple[0]:
                wells.append(source_well)
        transformation_wells = [well for i, well in enumerate(
            wells) if wells.index(well) == i]
        return transformation_wells

        
    def tiprack_slots(spotting_tuples, max_spot_vol=5):
        """Calculates p10 and p300 tiprack slots required.

        Args:
        spotting_tuples (list): Sets of spotting reactions are given 
            in the form: ((source wells), (target wells), (spotting volumes)). 
            Each unique transformation well is resuspended once prior to spotting.
        max_spot_vol (float): Maximum volume that is spotted per spot reaction.

        """
        # Reactions' number
        transformation_reactions = len(
            generate_transformation_wells(spotting_tuples))
        spotting_reactions = 0
        for spotting_tuple in spotting_tuples:
            spots = np.array(spotting_tuple[2])/max_spot_vol
            np.ceil(spots)
            spotting_reactions = spotting_reactions + int(np.sum(spots))

        # p10 tiprack slots
        p10_tips = transformation_reactions + spotting_reactions
        p10_tiprack_slots = p10_tips // 96 + 1 if p10_tips % 96 > 0 else p10_tips / 96

        # p300 tiprack slots
        p300_tips = transformation_reactions + spotting_reactions
        p300_tiprack_slots = p300_tips // 96 + \
            1 if p300_tips % 96 > 0 else p300_tips / 96
        return int(p10_tiprack_slots), int(p300_tiprack_slots)

    
    def transformation_setup(transformation_wells):
        """Sets up transformation reactions

        Args:
        transformation_wells (list). 

        """

        # Constants
        TEMP = 4  # Incubation temperature.
        ASSEMBLY_VOL = 5  # Volume of final assembly added to competent cells.
        MIX_SETTINGS = (4, 5)  # Mix after setting during final assembly transfers.
        INCUBATION_TIME = 20  # Cells and final assembly incubation time.

        # Set temperature deck to 4 °C and load competent cells
        tempdeck.set_temperature(TEMP)
        protocol.pause()
        protocol.comment('Load competent cells, uncap and resume run')

        assembly_plate_transformation_wells = [assembly_plate.wells_by_name()[well]
        for well in transformation_wells]
        transformation_plate_transformation_wells = [transformation_plate.wells_by_name()[well]
        for well in transformation_wells]

        # Transfer final assemblies
        p10_pipette.transfer(ASSEMBLY_VOL,
                            assembly_plate_transformation_wells,
                            transformation_plate_transformation_wells,
                            new_tip='always',
                            mix_after=(MIX_SETTINGS))
        # Incubate for 20 minutes and remove competent cells for heat shock
        protocol.delay(minutes=INCUBATION_TIME)
        #legacy_api.api.robot.pause()
        protocol.pause()
        protocol.comment('Remove transformation reactions, conduct heatshock and replace.')
        #legacy_api.api.robot.comment(
            #'Remove transformation reactions, conduct heatshock and replace.')

    def phase_switch(comment='Remove final assembly plate. Introduce agar tray and deep well plate containing SOC media. Resume run.'):
        """Function pauses run enabling addition/removal of labware.

        Args:
        comment (str): string to be displayed during run following pause.

        """
        #legacy_api.api.robot.pause()
        protocol.pause()
        #legacy_api.api.robot.comment(comment)
        protocol.comment(comment)


    def outgrowth(
            cols,
            soc_well):
        """Outgrows transformed cells.

        Args:
        cols (list of str): list of cols in transformation plate containing samples.
        soc_well (str): Well containing SOC media in relevant plate.

        """

        # Constants
        SOC_VOL = 125
        SOC_MIX_SETTINGS = (4, 50)
        TEMP = 37
        OUTGROWTH_TIME = 60
        SOC_ASPIRATION_RATE = 25
        P300_DEFAULT_ASPIRATION_RATE = 150

        # Define wells
        transformation_cols = transformation_plate.rows()[0][cols]
        soc = soc_plate.wells_by_name()[soc_well]

        # Add SOC to transformed cells
        p300_pipette.flow_rate.aspirate = SOC_ASPIRATION_RATE
        p300_pipette.transfer(SOC_VOL, soc, transformation_cols,
                            new_tip='always', mix_after=SOC_MIX_SETTINGS)
        p300_pipette.flow_rate.aspirate = P300_DEFAULT_ASPIRATION_RATE

        # Incubate for 1 hour at 37 °C
        tempdeck.set_temperature(TEMP)
        protocol.delay(minutes=OUTGROWTH_TIME)
        tempdeck.deactivate()


    def spotting_cols(spotting_tuples):
        """Evaluates spotting_tuples and returns unique cols (str) 
        associated with each spotting_tuple's source wells.

        Args:
        spotting_tuples (list): Sets of spotting reactions are given 
            in the form: ((source wells), (target wells), (spotting volumes)). 
            Each unique transformation well is resuspended once prior to spotting.

        """
        cols_list = []
        for spotting_tuple in spotting_tuples:
            source_wells_cols = [source_well[1:]
                                for source_well in spotting_tuple[0]]
            unique_cols = [col for i, col in enumerate(
                source_wells_cols) if source_wells_cols.index(col) == i]
            cols_list.append(unique_cols)
        return cols_list


    def spot_transformations(
            spotting_tuples,
            dead_vol=2,
            spotting_dispense_rate=0.025,
            stabbing_depth=2,
            max_spot_vol=5):
        """Spots transformation reactions.

        Args:
        spotting_tuples (list): Sets of spotting reactions are given 
            in the form: ((source wells), (target wells), (spotting volumes)). 
            Each unique source well is resuspended once prior to spotting.
        dead_vol (float): Dead volume aspirated during spotting.
        spotting_dispense_rate (float): Rate p10_pipette dispenses at during spotting.
        stabbing_depth (float): Depth p10_pipette moves into agar during spotting.
        max_spot_vol (float): Maximum volume that is spotted per spot reaction. 

        """

        def spot(
                source,
                target,
                spot_vol):
            """Spots an individual reaction using the p10 pipette.

            Args:
            source (str): Well containing the transformation reaction to be spotted.
            target (str): Well transformation reaction is to be spotted to.
            spot_vol (float): Volume of transformation reaction to be spotted (uL).  

            """
            # Constants
            DEFAULT_HEAD_SPEED = {'x': 400, 'y': 400,
                                'z': 125, 'a': 125}
            SPOT_HEAD_SPEED = {'x': 400, 'y': 400, 'z': 125,
                            'a': 125 // 4}
            DISPENSING_HEIGHT = 5
            SAFE_HEIGHT = 15  # height avoids collision with agar tray.

            # Spot
            p10_pipette.pick_up_tip()
            p10_pipette.aspirate(spot_vol + dead_vol, source)
            p10_pipette.move_to(target.top(SAFE_HEIGHT))
            p10_pipette.move_to(target.top(DISPENSING_HEIGHT))
            p10_pipette.dispense(volume=spot_vol, rate=spotting_dispense_rate)
            for key in SPOT_HEAD_SPEED.keys():
                protocol.max_speeds[key] = SPOT_HEAD_SPEED[key]
            p10_pipette.move_to(target.top(-1 * stabbing_depth))
            for key in DEFAULT_HEAD_SPEED.keys():
                protocol.max_speeds[key] = DEFAULT_HEAD_SPEED[key]
            p10_pipette.move_to(target.top(SAFE_HEIGHT))

            # Dispose of dead volume and tip
            p10_pipette.dispense(dead_vol, spotting_waste)
            p10_pipette.blow_out()
            p10_pipette.drop_tip()

        def spot_tuple(spotting_tuple):
            """Spots all reactions defined by the spotting tuple. Requires the function spot.

                Args:
                spotting_tuple (tuple): Spotting reactions given in the form: (source wells), (target wells), (spotting volumes).

            """
            source_wells = spotting_tuple[0]
            target_wells = spotting_tuple[1]
            spot_vols = list(spotting_tuple[2])
            while max(spot_vols) > 0:
                for index, spot_vol in enumerate(spot_vols):
                    if spot_vol == 0:
                        pass
                    else:
                        vol = spot_vol if spot_vol <= max_spot_vol else max_spot_vol
                        spot(transformation_plate.wells_by_name()[source_wells[index]],
                            agar_plate.wells_by_name()[target_wells[index]], vol)
                        spot_vols[index] = spot_vols[index] - vol

        # Constants
        TRANSFORMATION_MIX_SETTINGS = [4, 50]

        # Spot transformation reactions
        for spotting_tuple in spotting_tuples:
            source_wells_cols = [source_well[1:]
                                for source_well in spotting_tuple[0]]
            unique_cols = [col for i, col in enumerate(
                source_wells_cols) if source_wells_cols.index(col) == i]
            for col in unique_cols:
                transformation_plate_mix_well = transformation_plate.rows()[0][int(col)-1]
                p300_pipette.pick_up_tip()
                p300_pipette.mix(TRANSFORMATION_MIX_SETTINGS[0],
                                TRANSFORMATION_MIX_SETTINGS[1],
                                transformation_plate_mix_well)
                p300_pipette.drop_tip()
            spot_tuple(spotting_tuple)

    # Run protocol

    # Constants
    CANDIDATE_P10_SLOTS = ['9', '2', '5']
    CANDIDATE_P300_SLOTS = ['3', '6']
    P10_TIPRACK_TYPE = 'opentrons_96_tiprack_10ul'
    P300_TIPRACK_TYPE = 'opentrons_96_tiprack_300ul'
    P10_MOUNT = 'right'
    P300_MOUNT = 'left'
    ASSEMBLY_PLATE_TYPE = 'biorad_96_wellplate_200ul_pcr'
    ASSEMBLY_PLATE_SLOT = '8'
    TEMPDECK_SLOT = '10'
    TRANSFORMATION_PLATE_TYPE = 'biorad_96_wellplate_200ul_pcr'
    SOC_PLATE_TYPE = 'usascientific_96_wellplate_2.4ml_deep'
    SOC_PLATE_SLOT = '7'
    TUBE_RACK_TYPE = 'opentrons_24_tuberack_nest_1.5ml_snapcap'
    TUBE_RACK_SLOT = '11'
    SPOTTING_WASTE_WELL = 'A1'
    AGAR_PLATE_TYPE = 'axygen_1_reservoir_90ml'
    AGAR_PLATE_SLOT = '1'

    # Tiprack slots
    p10_p300_tiprack_slots = tiprack_slots(spotting_tuples)
    p10_slots = CANDIDATE_P10_SLOTS[
        :p10_p300_tiprack_slots[0]]
    p300_slots = CANDIDATE_P300_SLOTS[
        :p10_p300_tiprack_slots[1]]

    # Define labware
    p10_tipracks = [protocol.load_labware(P10_TIPRACK_TYPE, slot)
                    for slot in p10_slots]
    p300_tipracks = [protocol.load_labware(P300_TIPRACK_TYPE, slot)
                    for slot in p300_slots]
    p10_pipette = protocol.load_instrument('p10_single',
        P10_MOUNT, tip_racks=p10_tipracks)
    p300_pipette = protocol.load_instrument('p300_multi',
        P300_MOUNT, tip_racks=p300_tipracks)

    assembly_plate = protocol.load_labware(ASSEMBLY_PLATE_TYPE, ASSEMBLY_PLATE_SLOT)
    tempdeck = protocol.load_module('temperature module', TEMPDECK_SLOT)
    #tempdeck = modules.load('tempdeck', TEMPDECK_SLOT)
    #transformation_plate = protocol.load_labware(TRANSFORMATION_PLATE_TYPE,
                                       # TEMPDECK_SLOT, share=True)
    transformation_plate = tempdeck.load_labware(TRANSFORMATION_PLATE_TYPE)
    soc_plate = protocol.load_labware(SOC_PLATE_TYPE, SOC_PLATE_SLOT)
    tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_SLOT)
    spotting_waste = tube_rack.wells_by_name()[SPOTTING_WASTE_WELL]
    agar_plate = protocol.load_labware(AGAR_PLATE_TYPE, AGAR_PLATE_SLOT)

    # Register agar_plate for calibration
    p10_pipette.transfer(1, agar_plate.wells_by_name()[
        'A1'], agar_plate.wells_by_name()['H12'], trash=False)
    #p10_pipette.pick_up_tip(p10_tipracks[0][0])

    # Run functions
    transformation_setup(generate_transformation_wells(spotting_tuples))
    phase_switch()
    spotting_tuples_cols = [col for cols in spotting_cols(
        spotting_tuples) for col in cols]
    unique_cols = [col for i, col in enumerate(
        spotting_tuples_cols) if spotting_tuples_cols.index(col) == i]
    outgrowth(unique_cols, soc_well=soc_well)
    spot_transformations(spotting_tuples)
