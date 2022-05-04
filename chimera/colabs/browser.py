
class createColabScript():
    """Create script file used to run colabfold"""
    def __init__(self, 
                scriptFilePointer, 
                extraPath, 
                url,
                injectJavaScriptList,
                transferFn,
                resultsFile):
        print("createColabScript", resultsFile)
        colabCommand = f'''
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtWebEngineWidgets import *
import sys
import os

class WebPage(QWebEnginePage):
    """override chooseFiles method so we can inject the filename"""
    def __init__(self, parent=None, transferFn=None):
        super(WebPage, self).__init__(parent)
        self.transferFn = transferFn

    def chooseFiles(self, mode, oldfiles, mimetypes):
        print("self.transferFn", self.transferFn)
        return [self.transferFn]

# The QMainWindow class provides a main application window
class MainWindow(QMainWindow):
    """ Create a minimal web browser which run alphafold
    using a colab book and retrieves the results
    """
    CHIMERA = 0
    PHENIX = 1
    TEST = 2

    def __init__(self, extraPath, url,
                 injectJavaScriptList, transferFn, resultsFile):
        super(MainWindow,self).__init__()
        self.extraPath = extraPath
        self.injectJavaScriptList = injectJavaScriptList
        self.transferFn=transferFn
        self.resultsFile=resultsFile
        ###
        # Phase 1 create browser and connect
        ###
        # The QWebEngineView class provides a widget 
        # that is used to view and edit web documents.
        self.browser = QWebEngineView()
        # URL to open in the browser
        self.browser = QWebEngineView()
        # inject my version of WebPage
        self.p = WebPage(self.browser, 
                        transferFn=transferFn)
        self.browser.setPage(self.p)  
        self.browser.setUrl(QUrl(url))

        ###
        # Phase 2 java script injection with data
        # and start execution
        ###
        # execute after page is loaded
        self.browser.page().loadFinished.connect(self.executeJavaScript)

        self.setCentralWidget(self.browser)
        # show browser
        self.show()

        ###
        # Phase 3 retrieve data
        ###

        ## In this minimal web browser downloads
        # are not automatic This signal is emitted whenever a 
        # download has been triggered. The download argument holds 
        # the state of the download. The download has to be explicitly accepted
        # with QWebEngineDownloadItem::accept() or it will be cancelled.
        print("call to on_downloadRequested")
        self.p.profile().downloadRequested.connect(
            self.on_downloadRequested
        )

    def executeJavaScript(self):
        """ Execute javascript command in self.command"""
        # p.runJavaScript("throw document.readyState ")  # For DEBUGING
        for javaCommand in self.injectJavaScriptList:
            self.p.runJavaScript(javaCommand)
        
    def on_downloadRequested(self, download):
        """ Handle downloads.
        In this minimal web browser downloads
        are not automatic This signal is emitted whenever a 
        download has been triggered. The download argument holds 
        the state of the download. The download has to be explicitly accepted
        with QWebEngineDownloadItem::accept() or it will be cancelled.
        """

        try:
            download.setPath(self.resultsFile)
            download.accept()
            download.finished.connect(self.foo)
            ##download.finished.connect(lambda:self.downloadfinished(outPutFile.rsplit('/')[-1]))
        except Exception as e:
            print("ERROR", e)

        # download.finished.connect(lambda: print("download finished"))

    def downloadfinished(self,filename):
        print(filename+' Downloaded!')

    def foo(self):
        """ Close colab"""
        print("finished")
        self.hide();
        self.close();


# construct  QApplication 
app = QApplication(sys.argv)
# 
# create QWidget
window = MainWindow(extraPath='{extraPath}',
                    url='{url}',
                    injectJavaScriptList={injectJavaScriptList},
                    transferFn='{transferFn}',
                    resultsFile='{resultsFile}',
                    )
app.exec_()
'''
        scriptFilePointer.write(colabCommand)


#function KeepClicking(){
#   console.log("Clicking");
#   document.querySelector("colab-toolbar-button#connect").click()
#}setInterval(KeepClicking,60000)        