#! /usr/bin/env -S python3 

import argparse

parser = argparse.ArgumentParser(description="script to compute section properties of square tubes")

# Add arguments
#parser.add_argument("-h", "--help", action="store_true", help="help message")
parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose mode")
parser.add_argument("-b", "--outer", type=float, default=2, help="outer dimension")
parser.add_argument("-t", "--wall", type=float, default=0.1, help="wall thickness")
parser.add_argument("-s", "--strength", type=float, default=24, help="yield stress ")

# Parse arguments
args = parser.parse_args()

#if args.help:
#    print(" python  square_tube.py  -b 2.0  -t 0.1  -s 24 ")
if args.verbose:
    print(f"Outer dimension: {args.outer}")
    print(f"Wall thickness : {args.wall}")
    print(f"Yield stress: {args.wall}")

b  = args.outer
t  = args.wall
Fy = args.strength

A = (b**2) - (b-2*t)**2        # cross section area
I = (b**4)/12 - (b-2*t)**4/12  # cross section second moment of area
Z = (b**3)/4  - (b-2*t)**3/4   # plastic section modulus
Mp = Z*Fy

print(f'      I          A          Z  ')
print(f' {I:10.4f} {A:10.4f}  {Z:10.4f} ')
