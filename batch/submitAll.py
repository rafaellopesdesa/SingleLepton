#!/bin/env python

import os, re, sys, ROOT, time, math, filecmp, shutil
def main():
    if len(sys.argv) != 2:
        print 'Usage: '
        print '   submitAll.py sample'
        sys.exit(0)
    fileList = open('list.txt')
    fileInfo = [item for sublist in [map(lambda x:'%s%s' % (info[1],x), os.listdir(info[1])) for info in [line.split() for line in fileList.readlines() if (len(line) > 1 and line.split()[0] == sys.argv[1])]] for item in sublist]
    fileList.close()
    print 'Doing sample %s in %d jobs' % (sys.argv[1], len(fileInfo))
    writeCondor('submit.%s.condor' % sys.argv[1], sys.argv[1], len(fileInfo))
    
def writeCondor(fileName, sampleName, numJobs):
    configFile = open(fileName, 'w')
    configFile.write('universe=grid\n')
    configFile.write('Grid_Resource=condor cmssubmit-r1.t2.ucsd.edu glidein-collector.t2.ucsd.edu\n')
    configFile.write('+DESIRED_Sites="T2_US_UCSD"\n')
    configFile.write('Should_Transfer_Files = YES\n')
#    configFile.write('WhenToTransferOutput = ON_EXIT\n')
#    configFile.write('Transfer_Input_Files = ../SingleLepton.tar.gz\n')
    configFile.write('Output = logs/condorLog.%s.$(Process).stdout\n' % sampleName)
    configFile.write('Error = logs/condorLog.%s.$(Process).stderr\n' % sampleName)
    configFile.write('Log = logs/condorLog.%s.$(Process).log\n' % sampleName)
    configFile.write('+Owner = undefined\n')
    configFile.write('notify_user = rclsa@fnal.gov\n')
    configFile.write('x509userproxy = /tmp/x509up_u31150\n')
    configFile.write('Executable = run.sh\n')
    configFile.write('Transfer_executable = True\n')
    configFile.write('getenv = True\n')
    configFile.write('arguments = $(Process) %s\n' % sampleName)
    configFile.write('Queue %d\n' % numJobs)
    configFile.close()

if __name__ == "__main__":
    main()
