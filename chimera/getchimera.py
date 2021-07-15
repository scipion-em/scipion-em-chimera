import sys

import os
import tkinter as tk
from tk_html_widgets import HTMLLabel
from urllib import request, parse
import re

DOMAIN = "https://www.cgl.ucsf.edu"
LICENSE_URL = DOMAIN + "/chimerax/cgi-bin/secure/chimerax-get.py"

def getTgz(version):

    return 'ChimeraX-%s.tar.gz' %version

def getFile(version):
    # Examples:
    #    https://www.rbvi.ucsf.edu/chimerax/cgi-bin/secure/chimerax-get.py?file=1.1/linux/ChimeraX-1.1.tar.gz
    #    https://www.rbvi.ucsf.edu/chimerax/cgi-bin/secure/chimerax-get.py?file=1.2/linux/ChimeraX-1.2.5.tar.gz

    file = ".".join(version.split(".")[:2])
    return '%s/linux/' % file + getTgz(version)

def getLicenseURLwithFile(version):
    return LICENSE_URL + '?file=' + getFile(version)

#http://www.rbvi.ucsf.edu/chimerax/cgi-bin/secure/chimerax-get.py?file=1.1/linux/ChimeraX-1.1.tar.gz
def getChimeraX(version):

    print("Getting %s version" % version)
    TGZ = getTgz(version)

    def licenceAccepted():
        """ Called when licence has been accepted"""

        root.destroy()
        data = {"file": getFile(version), "choice": "Accept"}
        data = parse.urlencode(data).encode()
        req = request.Request(LICENSE_URL, data)
        resp = request.urlopen(req).read().decode('utf-8')
        # Get the download url
        match = re.search('url=(.*)"', resp)
        if not match:
            raise Exception("Response does not contain expected pattern: %s" % resp)

        downloadUrl = match.group(1)
        downloadUrl = DOMAIN + downloadUrl
        print("Downloading chimera from " + downloadUrl)
        os.system('wget "%s"  --progress=bar -c -O %s' % (downloadUrl, TGZ))
        sys.exit()

    if os.path.exists(TGZ):
        print(TGZ, " found. Skipping download.")
        sys.exit()

    licenseHTML = request.urlopen(getLicenseURLwithFile(version)).read().decode()

    licenseHTML = re.sub("<style.*/style>", "", licenseHTML, flags=re.DOTALL)
    licenseHTML = re.sub("<script.*/script>", "", licenseHTML, flags=re.DOTALL)

    root = tk.Tk()
    root.title("ChimeraX Binary Installation: License Page")
    frame = root
    frame.rowconfigure(0, weight=1)
    html_label = HTMLLabel(frame, html=licenseHTML)
    html_label.grid(row=0, column=0, sticky="news")
    html_label.fit_height()

    acceptButton = tk.Button(frame, text="Accept", command=licenceAccepted)
    acceptButton.grid(row=1, column=0, sticky="s")

    root.mainloop()

    print("License not accepted.")
    sys.exit(1)


if __name__ == '__main__':
    try:
        if len(sys.argv) != 2:
            print("must pass version to download:  1.1 or 1.2,...")
            sys.exit(2)
        version =  sys.argv[1]
        getChimeraX(version)
    except Exception as e:
        print("There was an error trying to launch Chimera license agreement page (%s)" % str(e))
        print('Please, follow this link (https://github.com/scipion-em/scipion-em-chimera) for instructions to connect Scipion with ChimeraX' )