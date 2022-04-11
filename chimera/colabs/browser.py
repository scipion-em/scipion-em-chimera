
class createColabScript():
    """Create script file used to run colabfold"""
    def __init__(self, 
                scriptFilePointer, 
                extraPath, 
                url,
                injectJavaScriptList,
                uploadfile=None):

        colabCommand = f'''
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtWebEngineWidgets import *
import sys
import os

class WebPage(QWebEnginePage):
    """override chooseFiles method so we can inject the filename"""
    def __init__(self, parent=None, uploadFile=None):
        super(WebPage, self).__init__(parent)
        self.uploadFile = uploadFile

    def chooseFiles(self, mode, oldfiles, mimetypes):
        #TODO: change returned file
        print('uploadFile', self.uploadFile)
        return ['/home/roberto/ScipionUserData/projects/TestPhenixProtSearchFit/Runs/000850_ProtImportPdb/extra/5ni1.cif']


# The QMainWindow class provides a main application window
class MainWindow(QMainWindow):
    """ Create a minimal web browser which run alphafold
    using a colab book and retrieves the results
    """
    CHIMERA = 0
    PHENIX = 1
    TEST = 2

    def __init__(self, extraPath, url,
                 injectJavaScriptList, uploadFile):
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
        self.browser = QWebEngineView()
        # inject my version of WebPage
        self.p = WebPage(self.browser, uploadFile=uploadFile)
        self.browser.setPage(self.p)  
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
        print("call to on_downloadRequested")
        self.p.profile().downloadRequested.connect(
            self.on_downloadRequested
        )

    def executeJavaScript(self):
        """ Execute javascript command in self.command"""
        # p.runJavaScript("throw document.readyState ")  # For DEBUGING
        self.p.runJavaScript(self.injectJavaScriptList[self.counter])
        self.counter += 1
        
    def on_downloadRequested(self, download):
        """ Handle downloads.
        In this minimal web browser downloads
        are not automatic This signal is emitted whenever a 
        download has been triggered. The download argument holds 
        the state of the download. The download has to be explicitly accepted
        with QWebEngineDownloadItem::accept() or it will be cancelled.
        """
        old_path = download.url().path()  # download.path()
        print("old_path", old_path)
        #path, _ = QFileDialog.getSaveFileName(
        #    self, "Save File", old_path, "*.zip"
        #)
        try:
            extraPath = os.path.join(self.extraPath, "results.zip")
            print("begin on_downloadRequested", extraPath)
            download.setPath(extraPath)
            # download.setPath(extraPath)
            download.accept()
            ##download.finished.connect(self.foo)
            download.finished.connect(lambda:self.downloadfinished(extraPath.rsplit('/')[-1]))

            from chimerax.core.commands import run # execute chimera cli command
            print("on_downloadRequested_exit")
            run(session, 'exit') # exit script
        except Exception as e:
            print("ERROR", e)

        # download.finished.connect(lambda: print("download finished"))

    def downloadfinished(self,filename):
        print(filename+' Downloaded!')

    def foo(self):
        """ Prepare files for Scipion"""
        print("finished")
        #self.hide();
        #self.close();


# construct  QApplication 
app = QApplication(sys.argv)
# 
# create QWidget
window = MainWindow(extraPath='{extraPath}',
                    url='{url}',
                    injectJavaScriptList={injectJavaScriptList},
                    uploadFile='{uploadfile}',
                    )
app.exec_()
'''
        scriptFilePointer.write(colabCommand)