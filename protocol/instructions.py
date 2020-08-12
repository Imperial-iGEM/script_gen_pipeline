# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 18:28:20
#  * @modify date 2020-07-17 18:28:20
#  * @desc [description]
#  */

from typing import List
from uuid import uuid4


class Transfer:
    """Transfer contents from one container to another.

    Args:
        src: the source Container
        dest: the destination Container
        volume: the volume to transfer in microliters
    """

    def __init__(self, src, dest, volume: float = 10.0):
        self.src = src
        self.dest = dest
        self.volume = volume

    def split(self, max_volume: float, multiple_of: float) -> List["Transfer"]:
        """Split a Transfer into multiple other Transfers based on max volume.

        This was necessitated by the Labcyte Echo that only transfers up to 10 uL
        and with transfers that are a multiple of 2.5 nL.

        Args;
            max_volume: the max volume of a single transfer in uL
            multiple_of: each transfer has to be a multiple of this (in uL)

        Returns:
            List[Transfer] -- list of transfers meeting restraints
        """

        transfer_count = math.ceil(self.volume / max_volume)
        volume_per_transfer = self.volume / transfer_count
        volume_per_transfer = round(volume_per_transfer / multiple_of) * multiple_of

        split_transfers: List["Transfer"] = []
        for _ in range(transfer_count):
            split_transfers.append(Transfer(self.src, self.dest, volume_per_transfer))
        return split_transfers

    def __hash__(self):
        return hash(self.src) + hash(self.dest) + hash(self.volume)


class Temperature:
    """A temperature instruction, has a temperature and time component.

    Args:
        temp: temperature in degrees Celcius
        time: length of time in seconds
    """

    def __init__(self, temp: float = 27.0, time: float = 10.0):
        self.temp = temp
        self.time = time


class Instruction:
    """A single instruction set generated by a single Step.
    Serves to track the history and interpretation of protocols.

    As stated in synbio, Instructions are generated by Steps. Each Step 
    calls this to add the step's output to the Protocol for accumulation.
    """

    def __init__(
        self,
        name: str = "",
        transfers: List[Transfer] = None,
        temps: List[Temperature] = None,
        instructions: List[str] = None,
    ):
        self.id = uuid4()
        self.name = name
        self.transfers = transfers
        self.temps = temps
        self.instructions = instructions or []

    def to_txt(self) -> str:
        txt = "not implemented"
        return txt


def instr_to_txt(name: str, instructions: List[Instruction]) -> str:
    """Return a text representation of the instructions for a protocol.

    ```txt
    Combinatorial MoClo:
    1. Mix the Assembly mix
        1.1. Create 200 µL 'assembly-mix' from 1:1 T4 ligase buffer (10X) and NEB Golden Gate Assembly Mix
    2. Setup the Setup Plate as specified
    ```

    Args:
        name: the protocol's name
        instructions: the list of protocol instructions

    Returns:
        a string representation of the protocol
    """

    txt = name + ":\n" if name else ""
    i = 1
    for instruction in instructions:
        instruction_txt = instruction.to_txt(i)
        if instruction_txt:
            txt += instruction_txt
            i += 1

    return txt
