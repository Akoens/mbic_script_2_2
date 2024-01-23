"""
    A Set of miscellaneous utilities.
"""

# Builtin.
from itertools import zip_longest
from typing import Iterable, Generator


def unpack_2di_to_list(some_iterable: Iterable[Iterable[any]]) -> list[any]:
    """
    Unpacks an iterable of iterables into a single list containing all the elements.
    :param some_iterable: The list of tuples to unpack.
    :return: An unpacked list with all the values of the some_iterable.
    """
    unpacked_list = [element for tup in some_iterable for element in tup]
    return unpacked_list


def sort_dict_to_list(some_dict: dict, sort_by: bool) -> list[tuple]:
    """
    Sort a dict by converting it to a list of key, value pairs and
    sorting that list by key or value.
    :param some_dict: A dict to sort.
    :param sort_by: A bool for sorting by key or by value.
    :return: A sorted list of key value pairs.
    """
    # Get all the key, value pairs.
    occurrences = list(some_dict.items())

    # Sort the key, values of dict
    sorted_key_value_list = sorted(occurrences, key=lambda x: x[int(sort_by)], reverse=True)
    return sorted_key_value_list


def find_matches(first_iter: Iterable, second_iter: Iterable) -> Generator:
    """
    Find the index of positional matches between two iterables.
    :param first_iter: First iterable to match.
    :param second_iter: Second iterable to match.
    :return: A Generator for finding the matches.
    """
    matches = (i for i, mach_pair in enumerate(zip_longest(first_iter, second_iter, " ")) if
               mach_pair[0] == mach_pair[1])
    return matches


def calc_prominence(elements: str or list[int or float or str],
                    element: str or int or float) -> float:
    """
    Calculates the percentage an element prominence in a string or list.
    :param elements: A str or list of any elements to compare.
    :param element: A value to be calculated it prominence of in elements.
    :return: The prominence percentage of the element in the list.
    """

    # Get total length.
    total_count_elements = len(elements)

    # Count the elements in the list.
    count = elements.count(element)

    # Calculate the prominence of an element in the list in percentages.
    prominence_percentage = round(count / total_count_elements * 100, 2)

    return prominence_percentage


def split_str(some_string: str, length: int, append: str = "") -> list[str]:
    """
    Splits a string into a list of segments by length.
    :param some_string: A string to be split into segments.
    :param length: The length of the segments.
    :param append: A string to append to every split segment.
    """
    split_segments = [some_string[i:i + length] + append for i in
                      range(0, len(some_string), length)]
    return split_segments
