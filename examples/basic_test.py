# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 15:11:47
#  * @modify date 2020-07-17 15:11:47
#  * @desc [description]
#  */


import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,os.path.join(parentdir, 'designs')) 

# import mymodule
print('currentdir',currentdir)
print('parentdir',parentdir)

print(__package__ is None)

# from script_gen_pipeline.designs.construct import Construct
from ..designs.construct import Construct
# from .protocol.protocol import Basic

if __name__ == "__main__":

    sbol_path_name = ""
    # sbol_input = sbol_path_name

    input_construct_path = '/Users/oliviagallup/Desktop/Kode/iGEM_2020/script_gen_pipeline/dna_bot_utils/examples/storch_et_al_cons.csv'
    output_sources_paths = '/Users/oliviagallup/Desktop/Kode/iGEM_2020/script_gen_pipeline/dna_bot_utils/examples/part_plate_2_230419.csv'

    design = Construct(input_construct_path, output_sources_paths)
    protocol = Basic(design)

    for i in range(len(protocol.scripts)):
        protocol.generate_ot_script(protocol.subprotocols[i].name,
                                    protocol.subprotocols[i].template_script,
                                    protocol.subprotocols[i].kwargs)
