#!/bin/python
import commands

d=commands.getoutput("grep rate\"\=\" "+sys.argv[-1]).split('\n')
# f=open("mdot.dat")
# d = f.readlines()
d = array([float(entry.split("rate")[1][1:].split('he')[0][:-1]) for entry in d])
# f.close()
print "<dotm>_t =",average(d)
