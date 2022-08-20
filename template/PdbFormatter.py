#!/usr/bin/env python
def splitta(something):
   txt = open(something)
   lines = [line for line in txt.readlines()]
   row = [val.split() for val in lines]
   return row

def fun(listaz):
    L=[]
    for val in listaz:
        L.append(Fun(val))
    return L
  
def Fun(lista):
    L = []
    for k in range(10):
        try:
            L.append(lista[k])
        except:
            L.append('')
    return L

def HTtoH1H2(lists):
    giove=fun(lists)
    L=[]
    for n in range(len(giove)):
      if giove[n][0] == 'HETATM':
        if giove[n][2] == 'OT':
            giove[n+1][2]='H1'
            giove[n+2][2]='H2'
            L.append(giove[n])
        else:
            L.append(giove[n])
    return L
          


def writer(file1,lista):
    f = open(file1,'w')
    for lines in lista:
        v0,v1,v2,v3,v4,v5,v6,v7,v8,v9 = Fun(lines)
#        print v0,v1,v2,v3,v4,v5,v6,v7,v8,v9
        f.write('%-6s %4s  %2s  %3s %5s %11s %7s %7s %1s %1s \n' % (v0,v1,v2,v3,v4,v5,v6,v7,v8,v9))
#        f.write(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9))
#        f.write('{0:6s} {1:>6s} {2:>4s} {3:>5s} {4:>5s} {5:>12s} {6:>7s} {7:>7s} {8:>7s} {9:>4s} \n '.format(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9))
    f.close()

def printwell(lists):
    for val in lists:
       print val

lol=HTtoH1H2(splitta('H2O.pdb_2'))
writer('file.out',lol)
