# Calculate area weighted monthly temperatures from map data
# Usage:
#  python temp.py my.dat > my.temp


import sys, math, numpy


# read a month of map data
def read_map( lines ):
  w = lines[0].split()
  month,year = sorted( [int(w[0]),int(w[1])] )
  date = year+month/12.0-1.0/24.0
  smap = [[numpy.nan for i in range(72)] for j in range(36)]
  for j in range(len(smap)):
    w = lines[j+1].split()
    for i in range(len(smap[j])):
      if not '.' in w[i]:
        t = 0.01*float(w[i])
      else:
        t = float( w[i] )
      if t > -99.0:
        smap[j][i] = t
  smap.reverse()
  return year, month, smap


# write a month of map data
def write_map( year, month, smap ):
  tmap = reversed( smap )
  lines = ["%4d %2d\n"%(year,month)]
  for row in tmap:
    s = ""
    for val in row:
      if not numpy.isnan(val):
        s += "%7.3f "%(val)
      else:
        s += "-99.9 "
    lines.append( s[:-1] + "\n" )
  return lines


# area of a latitude band by index
def areas( grid ):
  area = grid*[0.0]
  for i in range(grid):
    area[i] = ( ( math.sin(math.radians(180.0*(i+1)/grid-90.0)) -
                  math.sin(math.radians(180.0*(i  )/grid-90.0)) ) /
                math.sin(math.radians(180.0/grid)) )
  return area


def meanSH( tmap1 ):
  # average over the cells
  area = areas( len(tmap1) )
  s0, s1 = 0.0, 0.0
  for lati in range(0,len(tmap1)/2):
    for lngi in range(len(tmap1[lati])):
      if not numpy.isnan(tmap1[lati][lngi]):
        a = area[lati]
        s0 += a
        s1 += a*tmap1[lati][lngi]
  return s1/s0

def meanNH( tmap1 ):
  # average over the cells
  area = areas( len(tmap1) )
  s0, s1 = 0.0, 0.0
  for lati in range(len(tmap1)/2,len(tmap1)):
    for lngi in range(len(tmap1[lati])):
      if not numpy.isnan(tmap1[lati][lngi]):
        a = area[lati]
        s0 += a
        s1 += a*tmap1[lati][lngi]
  return s1/s0


# MAIN PROGRAM
# default values
datafile1 = None
if len(sys.argv) > 1:
  datafile1 = sys.argv[1]

# read data
nmonths = 9999
f = open( datafile1 )
lines1 = f.readlines()
f.close()
nmonths = min(nmonths,len(lines1)/37)

# calculate maps
maps = []
for m in range(nmonths):
  # read land data
  year1,month1,tmap1 = read_map( lines1[37*m:37*m+37] )
  n,s = meanNH(tmap1), meanSH(tmap1)
  print year1+month1/12.0-1.0/24.0, 0.5*(n+s), n, s

