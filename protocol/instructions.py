# /**
#  * @author Olivia Gallupová
#  * @email olivia.gallupova@gmail.com
#  * @create date 2020-07-17 18:28:20
#  * @modify date 2020-07-17 18:28:20
#  * @desc [description]
#  */


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