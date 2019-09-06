"""
Name: operator

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Functions:
(1) Determine whether two bases are equal or intersect;
(2) Fuse two bases.
"""

from inherent import *


def is_same_base(base_1, base_2):
    """
    Equal judgement of two bases.
    Here, equality can mean that two characters of base are equal,
    or a degenerate base contains another base.

    :param base_1: one base used for comparison.
    :param base_2: another base used for comparison.

    :return: the result of judgement.
    """
    if base_1 is None or base_2 is None:
        return False

    # calculate same bases.
    if base_1 == base_2:
        return True

    # calculate inclusion.
    base_1_values = list(base_mapping[0].get(base_1))
    base_2_values = list(base_mapping[0].get(base_2))

    # judge whether there is an intersection between the two base values.
    if len(list(set(base_1_values).intersection(set(base_2_values)))) > 0:
        return True

    return False


def fused_bases(base_1, base_2):
    """
    Get the intersecting base from two base.

    :param base_1: one base used for fusion.
    :param base_2: another base used for fusion.

    :return: the base fused by two input base.
    """
    if base_1 == base_2:
        return base_1

    base_1_values = list(base_mapping[0].get(base_1))
    base_2_values = list(base_mapping[0].get(base_2))

    fused_values = sorted(list(set(base_1_values).intersection(set(base_2_values))))

    if len(fused_values) > 0:
        fused_base = base_mapping[1].get("".join(fused_values))
    else:
        fused_base = None

    return fused_base
