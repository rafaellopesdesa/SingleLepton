#!/bin/env python

import os, sys

fileList = open('list.txt')
fileInfo = [item for sublist in [map(lambda x:'%s%s' % (info[1],x), os.listdir(info[1])) for info in [line.split() for line in fileList.readlines() if (len(line) > 1 and line.split()[0] == sys.argv[1])]] for item in sublist]
fileList.close()

newFile = open('file.list', 'w')
newFile.write(fileInfo[int(sys.argv[2])])
newFile.close()
