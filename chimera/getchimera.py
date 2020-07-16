import sys

import webview
import os

TGZ = 'ChimeraX-1.0.tar.gz'


def getChimeraX():

    def on_loaded():
        window = webview.windows[0]
        loadedUrl = window.get_current_url()
        if loadedUrl == "https://www.cgl.ucsf.edu/chimerax/cgi-bin/secure/chimerax-get.py":
            anchors = window.get_elements("a")
            downloadAnchor = anchors[0]
            downloadUrl = downloadAnchor['href']
            window.destroy()
            if "chimerax-get.py" in downloadUrl:
                os.system('wget "%s" -nv -c -O %s' %(downloadUrl, TGZ))
            else:
                print("Licence declined. Installation can't continue.")

    if os.path.exists(TGZ):
        print(TGZ, " found. Skipping download.")
        sys.exit()

    # Create a standard webview window
    window = webview.create_window('Chimera license agreement',
                                   'https://www.cgl.ucsf.edu/chimerax/cgi-bin/secure/chimerax-get.py?file=linux/ChimeraX-1.0.tar.gz',
                                   on_top=True)
    window.loaded += on_loaded
    webview.start(gui='qt')

if __name__ == '__main__':
    try:
        getChimeraX()
    except Exception as e:
        print("There was an error trying to launch Chimera license agreement page (%s)" % str(e))
        print('Please, follow this link (https://github.com/scipion-em/scipion-em-chimera) for instructions to connect Scipion with ChimeraX' )