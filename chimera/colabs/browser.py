
class createColabScript():
    """Create script file used to run colabfold"""
    def __init__(self, sequence, scriptFilePointer, extraPath, flavour, url):
        
        colabCommand = f'''
        from PyQt5.QtCore import *
        from PyQt5.QtWidgets import *
        from PyQt5.QtGui import *
        from PyQt5.QtWebEngineWidgets import *

        # The QMainWindow class provides a main application window
        class MainWindow(QMainWindow):
            """ Create a minimal web browser which run alphafold
            using a colab book and retrieves the results
            """
            CHIMERA = 0
            PHENIX = 1

            def __init__(self, sequence, extraPath, flavour, url):
                #TODO; pass colab book URL
                super(MainWindow,self).__init__(*args, **kwargs)
                self.sequence = sequence
                self.extraPath = extraPath
                # The QWebEngineView class provides a widget 
                # that is used to view and edit web documents.
                self.browser = QWebEngineView()
                # URL to open in the browser
                self.browser.setUrl(QUrl(url))
                # execute after page is loaded
                self.browser.page().loadFinished.connect(self._set_colab_sequence)

                self.setCentralWidget(self.browser)
                # execute after page is loaded
                if flavour == self.CHIMERA:
                    self.browser.page().loadFinished.connect(self._run_colab) 
                # show browser
                self.show()

                ## In this minimal web browser downloads
                # are not automatic This signal is emitted whenever a 
                # download has been triggered. The download argument holds 
                # the state of the download. The download has to be explicitly accepted
                # with QWebEngineDownloadItem::accept() or it will be cancelled.
                p = self.browser.page()
                p.profile().downloadRequested.connect(
                    self.on_downloadRequested
                )

            #def update_title(self): 
            #    title = self.browser.page().title() 
            #    self.setWindowTitle("Initializacion finished") 
        
            def _set_colab_sequence(self, sequence):
                """ Upload sequece to colab."""
                p = self.browser.page()
                #p.runJavaScript("throw document.readyState ")
                command = ('document.querySelector("paper-input").setAttribute("value", "%s");' % self.sequence)  + 'document.querySelector("paper-input").dispatchEvent(new Event("change"));'
                set_seq_javascript = (command)
                p.runJavaScript(set_seq_javascript)

                # document.querySelectorAll("paper-input")[1].setAttribute("value", "kkgg"); for second entry


            def _run_colab(self):
                """ press start button. TODO not sure if we can use this for phenix""" 
                p = self.browser.page()
                p.runJavaScript('document.querySelector("colab-run-button").click()')
                
            #@QtCore.pyqtSlot("QWebEngineDownloadItem*")
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
        app = QApplication([])
        # create QWidget
        window = MainWindow(sequence={sequence},
                            extraPath={extraPath},
                            url={url},
                            flavour={flavour}  # Flavour.CHIMERA
                            )
        app.exec_()
        '''
        scriptFilePointer.write(colabCommand)