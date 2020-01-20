try:  # python 2
    import Tkinter as tk
    import tkFont
    import ttk
except ImportError:  # Python 3
    import tkinter as tk
    import tkinter.font as tkFont
    import tkinter.ttk as ttk

import collections
import json

textFont1 = ("Arial", 10, "bold italic")
textFont2 = ("Arial", 16, "bold")
textFont3 = ("Arial", 8, "bold")

kkk = None


class LabelWidget(tk.Entry):
    def __init__(self, master, x, y, text):
        self.text = tk.StringVar()
        self.text.set(text)
        tk.Entry.__init__(self, master=master)
        self.config(relief="ridge", font=textFont1,
                    bg="#ffffff000", fg="#000000fff",
                    readonlybackground="#ffffff000",
                    justify='center', width=8,
                    textvariable=self.text,
                    state="readonly")
        self.grid(column=x, row=y)


class EntryWidget(tk.Entry):
    def __init__(self, master, x, y):
        tk.Entry.__init__(self, master=master)
        self.value = tk.StringVar()
        self.config(textvariable=self.value, width=8,
                    relief="ridge", font=textFont1,
                    bg="#ddddddddd", fg="#000000000",
                    justify='center')
        self.grid(column=x, row=y)
        self.value.set("")


class EntryGrid:
    """ Dialog box with Entry widgets arranged in columns and rows."""

    def __init__(self, colList, rowList, form, formAttribute, title="Entry Grid"):
        self.cols = colList[:]
        self.colList = colList[:]
        self.colList.insert(0, "")
        self.rowList = rowList
        self.win = tk.Toplevel()
        self.form = form
        self.formAttribute = formAttribute

        self.win.wm_title(title)

        self.mainFrame = tk.Frame(self.win)
        self.mainFrame.config(padx='3.0m', pady='3.0m')
        self.mainFrame.grid()
        self.make_header()

        self.gridDict = collections.OrderedDict()
        for i in range(1, len(self.colList)):
            for j in range(len(self.rowList)):
                w = EntryWidget(self.mainFrame, i, j + 1)
                self.gridDict[(i - 1, j)] = w.value

                def handler(event, col=i - 1, row=j):
                    return self.__entryhandler(col, row)

                w.bind(sequence="<FocusOut>", func=handler)

        self.okButton = tk.Button(self.mainFrame, text="OK",
                                  command=self.quit).grid(row=j + 2, column=1, pady=4)
        self.data = collections.OrderedDict()

    def __entryhandler(self, col, row):
        s = self.gridDict[(col, row)].get()
        if s.upper().strip() == "EXIT":
            self.win.destroy()
        elif s.upper().strip() == "TEST":
            self.test()
        elif s.upper().strip() == "ADENO":
            self.ico()

    def ico(self):
        """ enter right values for icosahedron """
        labelDict = collections.OrderedDict()

        labelDict['A'] = 'h1'
        labelDict['B'] = 'h1'
        labelDict['C'] = 'h1'
        labelDict['D'] = 'h2'
        labelDict['E'] = 'h2'
        labelDict['F'] = 'h2'
        labelDict['G'] = 'h3'
        labelDict['H'] = 'h3'
        labelDict['I'] = 'h3'
        labelDict['J'] = 'h4'
        labelDict['K'] = 'h4'
        labelDict['L'] = 'h4'
        labelDict['M'] = 'p'
        labelDict['N'] = 'iiia'
        labelDict['O'] = 'viiiO'
        labelDict['P'] = 'viiiP'
        labelDict['Q'] = 'lh3psdeudo'
        labelDict['R'] = 'lh3psdeudo'
        labelDict['S'] = 'lh3psdeudo'
        labelDict['T'] = 'lh3sym'

        counter = 0
        for row in self.rowList:
            # print counter, row, labelDict[row]
            self.set(0, counter, labelDict[row])
            counter += 1

    def test(self):
        """ enter right values for icosahedron """
        for i in range(len(self.cols)):
            for j in range(len(self.rowList)):
                self.set(i, j, "")
                self.win.update_idletasks()
                self.set(i, j, i + 1 + j)
                self.win.update_idletasks()

    def __str__(self):
        return "data = {}".format(self.data)

    def quit(self):
        for k, v in zip(self.rowList, self.gridDict.values()):
            # for k, v in zip(self.rowList, list(self.gridDict.values())):
            self.data[k] = v.get()
        value = json.dumps(self.data)
        self.form.setVar(self.formAttribute, value)
        self.win.destroy()

    def make_header(self):
        self.hdrDict = {}
        for i, label in enumerate(self.colList):
            def handler(event, col=i, row=0, text=label):
                return self.__headerhandler(col, row, text)

            w = LabelWidget(self.mainFrame, i, 0, label)
            self.hdrDict[(i, 0)] = w
            w.bind(sequence="<KeyRelease>", func=handler)

        for i, label in enumerate(self.rowList):
            def handler(event, col=0, row=i + 1, text=label):
                return self.__headerhandler(col, row, text)

            w = LabelWidget(self.mainFrame, 0, i + 1, label)
            self.hdrDict[(0, i + 1)] = w
            w.bind(sequence="<KeyRelease>", func=handler)

    def __headerhandler(self, col, row, text):
        """ has no effect when Entry state=readonly """
        self.hdrDict[(col, row)].text.set(text)

    def get(self, x, y):
        return self.gridDict[(x, y)].get()

    def set(self, x, y, v):
        self.gridDict[(x, y)].set(v)
        return v
