# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-19 16:26:15
#  * @modify date 2020-07-19 16:26:15
#  * @desc [description]
#  */

"""Steps: lab processes to assemble a design."""

from collections import defaultdict
from typing import Sequence, List, Dict, Optional, Callable

from script_gen_pipeline.labware.containers import content_id, Content, Container, Fridge
from script_gen_pipeline.protocol.instructions import Instruction, Temperature, Transfer
from script_gen_pipeline.protocol.biochem_utils import Reagent, Species


class Step:
    """A single Step of a Protocol. Taken
    Straight out of synbio."""

    def __init__(self):
        self.transfers: List[Transfer] = []
        self.temps: List[Temperature] = []
        self.instructions: List[str] = []
        self.inputs: Dict[str, Tuple[Content, float]]

    def __call__(self, protocol: "Protocol"):
        """Execute a step on a protocol, mutating containers and instructions."""
        protocol.add_instruction

        raise NotImplementedError


class Setup(Step):
    """Create a list of transfers to make a list of setup containers

    All transfers for a Setup step come from a Container (Fridge, )

    Attributes:
        target: target container list/ordering (default: {None})
        dest: the type of target container (default: {first target container})
        name: the name of this step in the protocol
        instructions: extra instructions to add to this step,
    """

    def __init__(
        self,
        target: Sequence[Container],
        dest: Container = None,
        name: str = "",
        instructions: List[str] = None,
    ):
        super().__init__()

        if not target:
            raise ValueError

        self.target = target
        if isinstance(target[0], list):
            print("[Setup] Target og a List[Container] and now List[List[Container]]")
        self.dest = dest
        self.name = name
        self.instructions = instructions if instructions else []

    def __call__(self, protocol: "Protocol"):
        """Create setup containers with enough contents to fill the target containers."""

        # accumulate the list of transfer volumes necessary
        id_to_content: Dict[str, Content] = {}
        id_to_volumes: Dict[str, List[float]] = defaultdict(list)
        # for container in self.target:
        for construct in self.target:
            for container in construct:
                for i, content in enumerate(container):
                    cid = content_id(content)
                    id_to_content[cid] = content
                    id_to_volumes[cid].append(container.volumes[i])

        # create the setup containers and transfers based on volume/count
        container = self.dest or (self.target[0][0] if isinstance(self.target[0], list) else self.target[0])
        volume_max = container.volume_max
        setup: List[Container] = []
        for cid, content in id_to_content.items():
            volumes = id_to_volumes[cid]

            new_container = container.create(
                content, volumes=[volumes.pop() + container.volume_dead]
            )
            while volumes:
                volume = volumes.pop()

                if new_container.volume() + volume > volume_max:
                    setup.append(new_container)
                    new_container = container.create(content, volumes=[volume])

                new_container.volumes[0] += volume
            setup.append(new_container)

        # make a transfer to fill each setup container
        transfers = [Transfer(src=Fridge(c), dest=c, volume=c.volume()) for c in setup]
        dest_name = type(setup[0]).__name__
        instructions = [
            f"Setup each {dest_name} with volumes (uL) specified"
        ] + self.instructions

        instruction = Instruction(
            name=self.name, transfers=transfers, instructions=instructions
        )

        protocol.add_instruction(instruction)
        protocol.containers = sorted(setup)

        return protocol


class Pipette(Step):
    """Pipette the input containers to match a 'target' list of containers

    Create transfers to map the input containers to a new list of output containers (target)

    Attributes:
        target: target container list/ordering (default: {None})
        name: the name of this step in the protocol
        instructions: extra instructions to add to this step
    """

    def __init__(
        self,
        target: Sequence[Container],
        name: str = "",
        instructions: List[str] = None,
    ):
        super().__init__()

        self.target = target
        self.name = name
        self.instructions = instructions if instructions else []


class Move(Step):
    """Move a fixed volume from each container to another container of same or new type

    Attributes:
        volume: the amount in millileters to move to the new container
        type: the type of new container to move the contents to
    """

    def __init__(self, volume: float, dest: Container = None, name: str = ""):
        super().__init__()

        self.name = name
        self.volume = volume
        self.dest = dest


class Add(Step):
    """Add contents to the existing containers.

    Attributes:
        add: the src container or content to add
        volume: the volume of the new content to add
        name: the name of this step in the protocol
        instructions: extra instructions to add to this step
    """

    def __init__(
        self,
        add: Content,
        volume: float,
        name: str = "",
        instructions: List[str] = None,
    ):
        super().__init__()

        assert add, "Must select Content to add to each container in Add Step"

        if isinstance(add, Container):
            self.add = add
        else:
            self.add = Fridge(add)
        self.volume = volume
        self.name = name
        self.instructions = instructions if instructions else []


class ThermoCycle(Step):
    """Thermo cycle for PCR, digestion/ligation, etc.

    Make a list of temperature instructions with the temperature
    to set it at, the length of time and the number of cycles.

    Example:
        An example of a ThermoCycle for PCR:
        >>> ThermoCycle(cycles=30, temps=[
        ...    Temperature(temp=97, time=5 * 60),  # denature
        ...    Temperature(temp=55, time=30),  # annealing
        ...    Temperature(temp=72, time=60),  # extension
        ... ])

    Attributes:
        temps: list of temperatures in a gradient
        name: the name of this step in the protocol
        cycles: the number of thermo cycles (default: {1})
        mutate:
            a function to mutate the contents of a container after thermo cycling.
            Used to anneal digested/ligated fragments or amplify DNA with primers
        instructions: list of additional instructions to add
        extension: last extension after other thermocycle steps.
            Example is a final 5 minute extension at the end of PCR that's common
    """

    def __init__(
        self,
        temps: List[Temperature],
        name: str = "",
        cycles: int = 1,
        mutate: Optional[Callable[[Container], Container]] = None,
        extension: Temperature = None,
        instructions: List[str] = None,
    ):
        super().__init__()

        self.name = name
        self.temps = temps
        self.cycles = cycles
        self.mutate = mutate
        self.extension = extension
        self.instructions = instructions if instructions else []


class Incubate(Step):
    """Incubate contents for some time at a set temperature.

    Incubate the contents of containers and in a fridge or some other incubator.

    Attributes:
        name: the name of this step in the protocol
        temp: the target the container should be in (default: {None})
        mutate:
            a function to mutate the contents of a container after thermo cycling.
            Used to anneal digested/ligated fragments or amplify DNA with primers
    """

    def __init__(
        self,
        temp: Temperature,
        name: str = "Incubate",
        mutate: Optional[Callable[[Container], Container]] = None,
    ):
        super().__init__()

        self.name = name
        self.temps = [temp]
        self.mutate = mutate


HeatShock: List[Step] = [
    Move(volume=3.0),
    Add(add=Species("E coli"), volume=10.0),
    ThermoCycle(name="Heat shock", temps=[Temperature(temp=42, time=30)]),
    Add(add=Reagent("SOC"), volume=150.0),
    Incubate(temp=Temperature(temp=37, time=3600)),
]
"""A composite HeatShock step for getting DNA into E coli.

>>> [
...    Move(volume=3.0),
...    Add(add=Species("E coli"), volume=10.0),
...    ThermoCycle(name="Heat shock", temps=[Temperature(temp=42, time=30)]),
...    Add(add=Reagent("SOC"), volume=150.0),
...    Incubate(temp=Temperature(temp=37, time=3600)),
... ]
"""
