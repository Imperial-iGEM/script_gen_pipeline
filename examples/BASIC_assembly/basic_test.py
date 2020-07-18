# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 15:11:47
#  * @modify date 2020-07-17 15:11:47
#  * @desc [description]
#  */

from ...input.construct import Construct
from ...protocol.protocol import Basic

sbol_path_name = ""
sbol_input = sbol_path_name

design = Construct(sbol_input)
protocol = Basic(design)