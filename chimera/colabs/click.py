# Is there any way to programmatically prevent Google Colab from disconnecting on a timeout?
# Google Colab notebooks have an idle timeout of 90 minutes and absolute timeout of 12 hours.
# This means, if user does not interact with his Google Colab notebook for more than 90 minutes,
#  its instance is automatically terminated. Also, maximum lifetime of a Colab instance is 12 hours.
# 
# Run this code in your Desktop, Then point mouse arrow over (colabs left panel - file section) 
# directory structure on any directory this code will keep clicking on directory on every 30 seconds
# so it will expand and shrink every 30 seconds so your session will not get expired Important 
# - you have to run this code in your pc
from pynput.mouse import Controller,Button
import time

mouse = Controller()
TITLE = (975, 10)
CLICK = (975, 121) 
while True:
    oldPosition = mouse.position
    print("oldPosition", oldPosition)
    mouse.position = TITLE
    mouse.click(Button.left,1)
    mouse.position = CLICK
    mouse.click(Button.left,1)
    #mouse.position = CLICK
    #mouse.click(Button.left,1)
    mouse.position = (oldPosition)
    mouse.click(Button.left,1)
    print('clicked')

    time.sleep(30)
