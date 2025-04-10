#!/usr/bin/env python3
import sys

def find_pairs(input: str) -> str:
    """ This function will check if a string of parenthesis is paired or not paired """
    return_string = ""
    
    # we'll use a stack data type
    # whenever we encounter a (, we add to the stack
    # when we don't, we pop from the stack for a matching parenthesis
    # if we can't pop from a stack, that means there's no matching parenthesis and it's not paired
    stack = []
    
    for character in input:
        # add to stack
        if character == "(":
            stack.append(character)
        # encountered a closing parenthesis
        else:
            # no matching
            if len(stack) == 0:
                return_string = "NOT PAIRED"
                return return_string
            stack.pop()
    
    # check stack one more time see if there's anything left in the stack
    if len(stack) == 0:
        return_string = "PAIRED"
    else:
        return_string = "NOT PAIRED"

    return return_string
    

print(find_pairs(sys.argv[1]))