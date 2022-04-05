
class createColabScript():
    """Create script file used to run colabfold"""
    def __init__(self, 
                scriptFilePointer, 
                extraPath, 
                url,
                injectJavaScriptList):
        print("createColabScript::init")
        colabCommand = f'''
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtWebEngineWidgets import *
import sys

# The QMainWindow class provides a main application window
class MainWindow(QMainWindow):
    """ Create a minimal web browser which run alphafold
    using a colab book and retrieves the results
    """
    CHIMERA = 0
    PHENIX = 1

    def __init__(self, extraPath, url, injectJavaScriptList):
        super(MainWindow,self).__init__()
        self.extraPath = extraPath
        self.injectJavaScriptList = injectJavaScriptList
        self.counter = 0
        ###
        # Phase 1 create browser and connect
        ###
        # The QWebEngineView class provides a widget 
        # that is used to view and edit web documents.
        self.browser = QWebEngineView()
        # URL to open in the browser
        self.browser.setUrl(QUrl(url))
        
        ###
        # Phase 2 java script injection with data
        # and start execution
        ###
        # execute after page is loaded
        for javaCommand in injectJavaScriptList:
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
        p = self.browser.page()
        p.profile().downloadRequested.connect(
            self.on_downloadRequested
        )

    def executeJavaScript(self):
        """ Execute javascript command in self.command"""
        p = self.browser.page()
        # p.runJavaScript("throw document.readyState ")  # For DEBUGING
        p.runJavaScript(self.injectJavaScriptList[self.counter])
        self.counter += 1
        
    def on_downloadRequested(self, download):
        """ Handle downloads.
        In this minimal web browser downloads
        are not automatic This signal is emitted whenever a 
        download has been triggered. The download argument holds 
        the state of the download. The download has to be explicitly accepted
        with QWebEngineDownloadItem::accept() or it will be cancelled.
        """
        download.setPath(self.extraPath)
        download.accept()
        download.finished.connect(self.foo)
        # download.finished.connect(lambda: print("download finished"))

    def foo(self):
        """ Prepare files for Scipion"""
        print("finished")
        self.hide();
        self.close();


# construct  QApplication 
app = QApplication(sys.argv)
# 
# create QWidget
window = MainWindow(extraPath='{extraPath}',
                    url='{url}',
                    injectJavaScriptList={injectJavaScriptList}
                    )
app.exec_()
'''
        scriptFilePointer.write(colabCommand)