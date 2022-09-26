try:
    from PyQt6 import QtCore, QtWidgets, QtWebEngineWidgets
except Exception:
    from PyQt5 import QtCore, QtWidgets, QtWebEngineWidgets



class Widget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(Widget, self).__init__(parent)

        self.view = QtWebEngineWidgets.QWebEngineView()
        self.view.page().profile().downloadRequested.connect(
            self.on_downloadRequested
        )
        url = "https://colab.research.google.com/github/scipion-em/scipion-em-chimera/blob/devel/chimera/colabs/test_colab.ipynb"
        self.view.load(QtCore.QUrl(url))
        hbox = QtWidgets.QHBoxLayout(self)
        hbox.addWidget(self.view)

    @QtCore.pyqtSlot("QWebEngineDownloadItem*")
    def on_downloadRequested(self, download):
        #old_path = download.url().path()  # download.path()
        #suffix = QtCore.QFileInfo(old_path).suffix()
        #path, _ = QtWidgets.QFileDialog.getSaveFileName(
        #    self, "Save File", old_path, "*." + suffix
        #)

        download.setPath("kk.zip")
        download.accept()


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    w = Widget()
    w.show()
    sys.exit(app.exec_())
