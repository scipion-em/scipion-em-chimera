# -*- coding: utf-8 -*-
import os
import sqlite3
import collections
import sys
import csv

# dictionary that relates chains with structures
labelDict = collections.OrderedDict()

firstValue = labelDict['A'] = 'h1'
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
labelDict['O'] = 'viiia'
labelDict['P'] = 'viiib'
labelDict['Q'] = 'lh3psdeudo'
labelDict['R'] = 'lh3a'
labelDict['S'] = 'lh3a'
labelDict['T'] = 'lh3sym'

# database file name
sqliteFN = "overlaps.sqlite"
tableName = "contacts"
pdbFileName = '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/cootOut0006_real_space_refined_renamed.cif'
# drop table and recrearte it 
commandDropTable = """DROP table IF EXISTS {}"""

commandCreateTable = """
CREATE TABLE {}(
     id integer primary key autoincrement,
     modelId_1  char(8),
     protId_1   char(8),
     chainId_1  char(8),
     aaName_1   char(3),
     aaNumber_1 int,
     atomId_1   char(8),
     modelId_2  char(8),
     protId_2   char(8),
     chainId_2  char(8),
     aaName_2   char(3),
     aaNumber_2 int,
     atomId_2   char(8),
     overlap float,
     distance float
     );"""

commandDropView = """DROP VIEW IF EXISTS {};"""
commandGroupAAView = """
CREATE VIEW {} AS
  SELECT modelId_1, 
         protId_1, 
         chainId_1, 
         aaName_1, 
         aaNumber_1, 
         atomId_1,
         modelId_2, 
         protId_2, 
         chainId_2, 
         aaName_2, 
         aaNumber_2,
         atomId_2,
         count(*), 
         ROUND(avg(overlap), 2), 
         ROUND(avg(distance),2) 
FROM contacts 
WHERE protId_1 = '{}'
GROUP BY modelId_1, 
         protId_1, 
         chainId_1, 
         aaName_1, 
         aaNumber_1, 
         atomId_1,
         modelId_2, 
         protId_2, 
         chainId_2, 
         aaName_2, 
         aaNumber_2, 
         atomId_2
ORDER BY protId_1, aaNumber_1, atomId_1, modelId_2, protId_2, 
                   aaNumber_2, atomId_2;

"""

commandEliminateDuplicates = """CREATE VIEW {} AS 
SELECT *
FROM {}

EXCEPT

SELECT cb.*
FROM {} ca, {} cb
WHERE ca.protId_1    = cb.protId_2
  AND cb.protId_1    = ca.protId_2
  AND ca.chainId_1   = cb.chainId_2
  AND cb.chainId_1   = ca.chainId_2
  AND ca.aaNumber_1  = cb.aaNumber_2
  AND cb.aaNumber_1  = ca.aaNumber_2
  AND ca.atomId_1  = cb.atomId_2
  AND cb.atomId_1  = ca.atomId_2
  AND ca.modelId_1   < cb.modelId_2


ORDER BY              protId_1, chainId_1, aaNumber_1, atomId_1, 
         modelId_2, protId_2, chainId_2, aaNumber_2, atomId_2;
"""
commandReport = """
SELECT * 
FROM {}
{}
"""


def connectDB(sqliteFN, tableName):
    conn = sqlite3.connect(sqliteFN)
    c = conn.cursor()
    c.execute(commandDropTable.format(tableName))
    c.execute(commandCreateTable.format(tableName))
    return c, conn


def createChimeraScript(labelDict, pdbFileName, firstValue):
    outFiles = []
    f = open("/tmp/chimera.cmd", "w")
    f.write("open {}\n".format(pdbFileName))
    f.write("sym #0 group i,222r contact 3\n")  # apply symmetry
    # rainbow model
    # findclash  #0:.Q,.R,.S test other savefile /home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/lh3a_test_v3.over  overlap -0.4 hbond 0.0 namingStyle simple
    protId = firstValue
    chains = ""
    comma = ''
    for k, v in labelDict.items():
        if protId == v:
            chains += "{}.{}".format(comma, k)
            comma = ','
            outFileBase = v
        else:
            outFile = "/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/{}.over".format(outFileBase)
            outFiles.append(outFile)
            f.write(
                """echo {}\nfindclash  #0:{} test other savefile {} overlap -0.4 hbond 0.0 namingStyle simple\n""".format(
                    chains, chains, outFile))
            protId = v
            chains = ".{}".format(k)
            outFileBase = v
    outFile = "/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/{}.over".format(outFileBase)
    outFiles.append(outFile)
    f.write("""echo {}\nfindclash  #0:{} test other savefile {} overlap -0.4 hbond 0.0 namingStyle simple\n""".format(
        chains, chains, outFile))
    f.close()
    os.system("/home/roberto/Scipion/scipion_plugin/software/em/chimera-1.13.1/bin/chimera --nogui /tmp/chimera.cmd")
    return outFiles


def parseFiles(outFiles, c):
    d = {}
    AND = " "
    WHERE = []
    for inFile in outFiles:
        print(inFile)
        counter = 0
        for line in open(inFile):
            if counter < 8:
                # print ("skip line", line
                counter += 1
            else:
                info = line.split()  # ['#2', 'TYR', '851.B', 'CE1', '#0.7', 'PRO', '78.N', 'CA', '3.059', '0.581']
                d['modelId_1'] = "'" + info[0] + "'"  # '#2'
                d['aaName_1'] = "'" + info[1] + "'"  # "'TYR'"
                info2 = info[2].split(".")  # ['851', 'B']
                d['aaNumber_1'] = info2[0]  # 851
                d['chainId_1'] = "'" + info2[1] + "'"  # B
                d['atomId_1'] = "'" + info[3] + "'"  # CE1
                d['protId_1'] = "'" + labelDict[info2[1]] + "'"
                d['modelId_2'] = "'" + info[4] + "'"  # '#0.7'
                d['aaName_2'] = "'" + info[5] + "'"  # PRO
                info2 = info[6].split(".")
                d['aaNumber_2'] = info2[0]  # 78
                d['chainId_2'] = "'" + info2[1] + "'"  # N
                d['protId_2'] = "'" + labelDict[info2[1]] + "'"
                d['atomId_2'] = "'" + info[7] + "'"  # CA
                d['overlap'] = info[8]  # 3.059
                d['distance'] = info[9]  # 0.5

                command = "INSERT INTO contacts "
                keys = "("
                values = " ("
                for key, value in d.items():
                    keys += key + ", "
                values += str(value) + ", "
                keys = keys[:-2] + ")"
                values = values[:-2] + ")"

                command += keys + " VALUES " + values
                ##print command
                c.execute(command)


def createReport(unique_value, c):
    WHERE = ""  # protId_2<>'h1'
    AND = "WHERE "
    for protein in unique_value:
        print("processing protein: ", protein)
        # get contact by prot_ID
        c.execute(commandDropView.format("view_" + protein))
        c.execute(commandGroupAAView.format("view_" + protein, protein))
        # eliminate duplications, most contacts are two times
        c.execute(commandDropView.format("view_ND_" + protein))
        c.execute(commandEliminateDuplicates.format("view_ND_" + protein,
                                                    "view_" + protein,
                                                    "view_" + protein,
                                                    "view_" + protein))
        # eliminate duplications, for h2 do not show h1 contacts
        outfile = "/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/{}.csv".format(protein)
        # c.execute('.header on')
        # c.execute('.mode csv')
        # c.execute('.once {}'.format(outfile))
        if sys.version_info < (3,):
            f = open(outfile, 'wb')
        else:
            f = open(outfile, 'w', newline="")

        data = c.execute(commandReport.format("view_ND_" + protein, WHERE))
        writer = csv.writer(f)  # ,delimiter=';')
        # write the rest
        writer.writerows(data)
        f.close()
        print(commandReport.format("view_ND_" + protein, WHERE))
        WHERE += "   {} protId_2<>'{}'\n ".format(AND, protein)
        AND = 'AND '


# execute chimera find clash
outFiles = createChimeraScript(labelDict, pdbFileName, firstValue)
# outFiles = ['/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/h1.over', '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/h2.over', '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/h3.over', '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/h4.over', '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/p.over', '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/iiia.over', '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/viiia.over', '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/viiib.over', '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/lh3a.over', '/home/roberto/latex/Documents/2019/paper_ladv2/OVERLAPS/lh3b.over']

# connect to database, delete table and recreate it
c, conn = connectDB(sqliteFN, tableName)

# parse all files created by chimera
parseFiles(outFiles, c)

# create report by prot_ID
# ordered list
seen = set()
unique_value = []
for x in labelDict.values():
    # for x in list(labelDict.values()):
    if x not in seen:
        unique_value.append(x)
        seen.add(x)
createReport(unique_value, c)

conn.commit()
conn.close()

# ---- anallyze one by one, put dict as ordered dict produce csv
