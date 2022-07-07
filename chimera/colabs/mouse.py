#!/usr/bin/env python3
# this script gives you the ccordinates of your mouse
print("""1) start the script
 2) you have 10 seconds to place the mouse in the right place
 3) wait and do NOT mouve the mouse untill the coordinates are printed 
 you may need to install the pynput module
""")

from pynput.mouse import Button, Controller
from time import sleep

mouse = Controller()

# wait 10 seconds
sleep(10)
# Read and print the pointer position
print('The current pointer position is {0}'.format(
    mouse.position))

