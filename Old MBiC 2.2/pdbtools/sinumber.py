"""
    A Collection of Functions for the Scientific notation of Large numbers.
    Created on Wed 1-11-2023 08:42
    Author: A F Koens
    Version: 0.0.1
"""

import math

# https://officespace.zendesk.com/hc/en-us/articles/115000593531-Scientific-Notation-Large-Numbers-Guide
CLASS1_SI_SUFFIX = {3: "K", 6: "M", 9: "B", 12: "T", 15: "Qa", 18: "Qi", 21: "Sx", 24: "Sp",
                    27: "Oc", 30: "N", 33: "Dc"}


def count_digits(number: float) -> int:
    """
    Count the amount of digits in a decimal number.
    :param number: Decimal number to be counted.
    :return: The amount of digits.
    """
    # Check if the number is 0
    if not number:
        return 1

    # Count the number of digits using a log10
    number_of_digits = int(math.log10(abs(number))) + 1

    return number_of_digits


def shift_decimal_left(number: float, n_places: int) -> float:
    """
    Shift the decimal places of a float n_places to the left.
    :param number: A float number to shift places of.
    :param n_places: An int number of places to shift.
    :return: The shifted float number.
    """

    # Shift the decimal by dividing the number by 10 pow [number of places to shift]
    decimal_shifted_number = number / (10 ** n_places)

    return decimal_shifted_number


def right_remove_digit(number: int, n_places: int) -> int:
    """
    Removes digits from the start of a number.
    :param number: A number to remove the digits of.
    :param n_places: The amount of positions to remove.
    :return: The remaining number.
    """

    # Shift the decimal by floor dividing the number by 10 pow [number of places to shift]
    # This shifts the decimal left and removes the float
    remaining_digits = number // (10 ** n_places)

    return remaining_digits


def abbreviate_float(number: float, precision: int = 0) -> str:
    """
    Abbreviate a float to a SI prefixed number.
    :param number: Number to abbreviate.
    :param precision: Precession to abbreviate to.
    :return: Human-readable string of float number.
    """

    #
    digits = count_digits(number)
    closest_si_index = digits - (digits % 3)
    remainder = shift_decimal_left(number, closest_si_index - 1)

    # Get the corresponding Scientific letter for big numbers
    si_suffix_char = CLASS1_SI_SUFFIX[closest_si_index]

    # Format the string with decimal precision
    formatted_string = f"{remainder:.{precision}f}{si_suffix_char}"

    return formatted_string
