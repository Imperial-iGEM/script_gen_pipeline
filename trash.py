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