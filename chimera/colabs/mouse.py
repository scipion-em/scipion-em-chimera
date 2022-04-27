from pynput.mouse import Button, Controller
from time import sleep

mouse = Controller()

# Read pointer position
sleep(10)
print('The current pointer position is {0}'.format(
    mouse.position))

