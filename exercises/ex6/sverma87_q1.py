#!/usr/bin/env python3
import sys

def createTriangle(symbol: str, length: int) -> str:
    """ Prints a triangle of text with the specified symbol """
    
    triangle = ""
    apex = 0 # the the tallest part of the triangle 
    apex = int(length) / 2
        
    # print the first half of the string
    symbol_print = 1
    for line in range(int(apex)):
        for num_of_symbols in range(symbol_print):
            triangle += symbol
        # this condition ensures that a new line won't be added for the last iteration of the loop
        if line != apex - 1:
            triangle += "\n"
        symbol_print += 1
    
    # get the resversed triangle
    reversed_triangle = triangle[::-1]
    
    # if the length is odd, at the middle apex once 
    apex_odd = apex + 1
    if int(length) % 2 != 0:
        for num_of_symbols in range(int(apex_odd)):
            triangle += symbol
    
    # if the length is odd, don't add a new line
    if int(length) % 2 != 0:
        triangle += reversed_triangle
    else:
        triangle += "\n" + reversed_triangle
    
    return triangle
    
# The only line of code not in a function is calling a function
print(createTriangle(sys.argv[1], sys.argv[2]))