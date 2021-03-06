#NOTES:
#1) sys.argv[1] must be GRASS lines vector
#2) Script must be run in GRASS mapset with EPSG: 4326 (WGS84) for Maxent formatting

import sys
import os
import grass.script as gs
from grass.pygrass.vector import Vector
from grass.pygrass import *
import sqlite3
import csv
import argparse

#Parse Arguments
parser = argparse.ArgumentParser(description='Convert vector lines to Maxent formatted points csv file')
parser.add_argument('damsIn', type=str, help='Input vector lines file')
parser.add_argument('pointsPerline', type=int, help='Number of points per line')
parser.add_argument('csvOutputdir', type=str, help='Output directory for Maxent formatted csv file')
#could add a v.clip call on a 4th argument: study area, which extracts dams in study area

args = parser.parse_args()

damsIn = args.damsIn
pointsPerline = args.pointsPerline
csvOutputdir = args.csvOutputdir


damLines = Vector(damsIn)
damLines.open(mode='r')
damLink = damLines.dblinks[0]
damAtts = damLink.table()
cursor = damAtts.execute()
lineCats = []
for row in cursor.fetchall():
    lineCats.append(row[0])

    

totalPoints = int(pointsPerline*len(lineCats))
percentIncrement = 100/int(pointsPerline)
rulesFile = open(csvOutputdir + '/rules.txt', mode='w')

pointID = 0
for lineCat in lineCats:
    percent = 0
    for point in range(0, int(pointsPerline)):
        ruleRow = 'P{0}{1}{0}{2}{0}{3}{4}'.format(' ', pointID, lineCat, percent,'%\n')
        rulesFile.write(ruleRow)
        percent = percent + percentIncrement
        pointID = pointID + 1



rulesFile.close()


#CREATE POINTS FROM LINE FEATURES

rules = os.path.realpath(rulesFile.name)
damPoints = damsIn + '_points'
gs.run_command('v.segment', input=damsIn, output=damPoints, rules=rules, verbose='True', overwrite='True')

#add attribute table and Maxent columns
gs.run_command('v.db.addtable', map=damPoints)
gs.run_command('v.to.db', map=damPoints, type='point', option='coor', columns='Long,Lat')

points = Vector(damPoints)
points.open(mode='r')
pointsLink = points.dblinks[0]
pointsAtts = pointsLink.table()
cursor_2 = pointsAtts.execute()


#OUTPUT TO CSV WITH MAXENT FORMATTING
csvOutputdir = csvOutputdir
csvOutputFile = os.path.join(csvOutputdir, damPoints + '_maxent_input.csv')

gs.run_command('db.out.ogr', input=damPoints, output=csvOutputFile, format='CSV', overwrite='True') 














