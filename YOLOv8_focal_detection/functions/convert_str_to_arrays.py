import numpy as np
import re

def string_to_numpy_array(array_str):
    '''
    Input:
    - array_str: A string representation of a numerical array.

    Action:
    - Extracts all numerical values from the input string using regular expressions.
    - Converts these string representations of numbers into floats.
    - Constructs and returns these numbers as a numpy array.

    Output:
    - Returns a numpy array of floats derived from the string.
    '''

    # Find all numbers in the string
    numbers = re.findall(r'[\d.]+', array_str)
    # Convert strings to float and then to a numpy array
    array = np.array(numbers, dtype=float)
    return array

def string_to_tuples_array(array_str):
    '''
    Input:
    - array_str: A string representation of a numerical array, possibly including Python list or tuple notation.

    Action:
    - Cleans the string by removing text and special characters that are not part of the numerical values.
    - Extracts all numerical values from the cleaned string.
    - Groups every two numbers into a tuple, assuming they represent coordinates or 2D points.
    - Returns these tuples as a list of tuples.

    Output:
    - Returns a list of tuples, each containing two float numbers, derived from the cleaned string.
    '''    
    # Remove 'array' and unnecessary characters
    clean_str = array_str.replace('array', '').replace('[', '').replace(']', '').replace('(', '').replace(')', '')
    # Extract all numbers
    numbers = re.findall(r'[\d.]+', clean_str)
    # Convert every two numbers into a tuple (assuming 2D points) and return as a list of tuples
    tuples_list = [(float(numbers[i]), float(numbers[i+1])) for i in range(0, len(numbers), 2)]
    return tuples_list