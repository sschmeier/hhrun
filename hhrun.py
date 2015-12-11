#!/usr/bin/env python
"""
NAME: hhrun.py
=========

DESCRIPTION
===========
## Information
A wrapper for ``hhsearch`` (http://wwwuser.gwdg.de/~compbiol/data/hhsuite/releases/). Run ``hhsearch`` on each sequence entry in a fasta-file using multiple distributed concurrent processes.

## Databases
Databases for ``hhsearch`` can be downloaded at http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/.

INSTALLATION
============

1. Download the ``hhsuite`` at http://wwwuser.gwdg.de/~compbiol/data/hhsuite/releases/.
2. Download databases at http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/.
3. Download or clone ``hhrun``.
4. Run ``hhrun.py``.

USAGE
=====

VERSION HISTORY
===============

0.1.1   2015/12/10    Minor documentation changes.
0.1.0   2015/12/10    Initial version.

LICENCE
=======

2015, copyright Sebastian Schmeier (s.schmeier@gmail.com), http://sschmeier.com
"""
__version__='0.1.1'
__date__='2015/12/10'
__email__='s.schmeier@gmail.com'
__author__='Sebastian Schmeier'
import sys, os, os.path, argparse, csv, collections, glob
import gzip, bz2, zipfile
import time, tempfile, subprocess, re, uuid

# When piping stdout into head python raises an exception
# Ignore SIG_PIPE and don't throw exceptions on it...
# (http://docs.python.org/library/signal.html)
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

class Seqrecord:
    """ One fastA record """
    def __init__(self):
        self.description = None
        self.seq = None

class FastaParser:
    """
    Used if BioPython is not installed.
    A relative quick hack.
    The BioPython implementation is better and should be used if possible.
    Nevertheless, it is working.
    """
    def __init__(self):
        """"""
        self.fileobj = None
        
    def parse(self, file_obj, type="fasta"):
        if type!="fasta":
            sys.stderr.write('Without BioPython only fasta-format allowed. EXIT.\n')
            sys.exit()
        self.fileobj = file_obj
        oReg = re.compile('>(.+)')
        aLines = [s.strip() for s in self.fileobj.readlines()]
        while len(aLines)>0:
            sLine = aLines.pop(0)
            oRes = oReg.search(sLine)
            if oRes:
                obj = Seqrecord()
                obj.description = oRes.group(1)
            else:
                obj.seq = sLine
                yield obj
        
def parse_cmdline():
    ## parse cmd-line -----------------------------------------------------------
    sDescription = 'Read a fasta-file and run for each entry hhsearch in parallel depending on the submited processes to use. It will look in the database folder for all database files and run hhsearch with all of them.' 
    sVersion='version %s, date %s' %(__version__,__date__)
    sEpilog = 'Copyright %s (%s)' %(__author__, __email__)

    oParser = argparse.ArgumentParser(description=sDescription,
                                      #version=sVersion,  # does not work in Python 3
                                      epilog=sEpilog)
    oParser.add_argument('sFile',
                         metavar='FILE',
                         help='Peptide sequences in fastA-format. [if set to "-" or "stdin" reads from standard in]')
    oParser.add_argument('--hh',
                         dest='sHH',
                         default='~/temp/hhsuite/bin',
                         metavar='PATH',
                         help='Path to the directory containing "hhsearch" executable. "hhsearch" needs to reside in the hhsuite directory, as it uses internally links to libraries in that directory and we are not setting any global variables (Download at: http://wwwuser.gwdg.de/~compbiol/data/hhsuite/releases/). [default = ~/temp/hhsuite/bin]')
    oParser.add_argument('--db',
                         dest='sDB',
                         default='~/temp/hhsuite/db',
                         metavar='PATH',
                         help='Path to directory containing databases for hhsearch (Download at: http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/). [default = ~/temp/hhsuite/db]')
    oParser.add_argument('-o', '--out',
                         metavar='STRING',
                         dest='sOut',
                         default=None,
                         help='Out-file. [default: "stdout"]')
    
    group1 = oParser.add_argument_group('HHsuite', 'optional arguments:')
    group1.add_argument('-n', '--number',
                         metavar='INT',
                         dest='iTop',
                         default=10,
                         type=int,
                         help='Number of top matches to return. [default: 10]')

    group2 = oParser.add_argument_group('Multithreading', 'optional arguments:')
  
    group2.add_argument('-p', '--processes',
                         metavar='INT',
                         type=int,
                         dest='iP',
                         default=1,
                         help='Number of sub-processes (workers) to use. It is only logical to not give more processes than cpus/cores are available. [default: 1]')
    group2.add_argument('-t', '--time',
                         action='store_true',
                         dest='bTIME',
                         default=False,
                         help='Time the runtime and print to stderr. [default: False]')
    group2.add_argument('--noprogress',
                         action='store_true',
                         dest='bNoProgress',
                         default=False,
                         help='Turn the progress-bar off. A progress-bar will force a "chunksize" of 1 for the threading. This might slow things down for very large job numbers, but allows for a realistic progress-bar. [default: Show progress-bar]')
    
    oArgs = oParser.parse_args()
    return oArgs, oParser

def load_file(s):
    """ LOADING FILES """
    if s in ['-', 'stdin']:
        oF = sys.stdin
    elif s.split('.')[-1] == 'gz':
        oF = gzip.open(s)
    elif s.split('.')[-1] == 'bz2':
        oF = bz2.BZFile(s)
    elif s.split('.')[-1] == 'zip':
        oF = zipfile.Zipfile(s)
    else:
        oF = open(s)
    return oF

def my_func(args):
    """
    THIS IS THE ACCTUAL WORKFUNCTION THAT HAS TO BE EXECUTED MULTPLE TIMES.
    This funion will be distributed to the cores requested.
    # do stuff
    res = ...
    return (args, res)
    """
    # Do stuff and create result
    # EXAMPLE: Here we add up arg1 and arg2 and wait a bit.
    sJob = str(args[0])
    sFasta = args[1]
    sDB = args[2]  # path to db
    sHH = args[3] # path to hhsearch
    sUUID = args[4]  # unique id
    
    #=======================================
    # Where to put intermidate tempfiles
    # can be changed to e.g. /temp
    sTempfileDir = os.getcwd()
    #=======================================

    #oFfasta = tempfile.NamedTemporaryFile()
    # should use tempfile here
    oFfasta = open(os.path.join(sTempfileDir,'job_%s_%s.in'%(sUUID,sJob)), 'w')
    oFfasta.write(sFasta)
    oFfasta.close()
    
    sCmd = "%s -i %s -d %s -o %s"%(os.path.join(sHH,'hhsearch'),
                                        oFfasta.name,
                                        sDB,
                                        os.path.join(sTempfileDir,'job_%s_%s.out'%(sUUID, sJob)))

    # standard out to /dev/null
    oFNULL = open(os.devnull, 'w')
    # hhsearch prints info to stdout, get rid of it => print to oFNULL
    subprocess.call(sCmd, shell=True, stdout=oFNULL)
    oFres = open(os.path.join(sTempfileDir, 'job_%s_%s.out'%(sUUID,sJob)))
    sHHres = oFres.read()
    oFres.close()
    # delete tempfiles manually, can be deleted if real tempfiles are used
    # can be turned off to keep the intermediate results
    os.remove(oFres.name)
    os.remove(oFfasta.name)
    # work on the result
    oRe = re.compile('>(.+?)\n(.+Sum_probs=.+?)\n')
    aRes = oRe.findall(sHHres)    
    return (args, aRes)

def main():
    oArgs, oParser = parse_cmdline()

    # Load input fastA-file
    oF = load_file(oArgs.sFile)

    if not oArgs.sOut:
        oFout = sys.stdout
    elif oArgs.sOut in ['-', 'stdout']:
        oFout = sys.stdout
    elif oArgs.sOut.split('.')[-1] == 'gz':
        oFout = gzip.open(oArgs.sOut, 'wb')
    else:
        oFout = open(oArgs.sOut, 'w')
 
    # ------------------------------------------------------
    #  THREADING
    # ------------------------------------------------------
    from timeit import default_timer as timer
    from multiprocessing import Pool
    
    # get number of subprocesses to use
    iNofProcesses = oArgs.iP
    if iNofProcesses<1: oParser.error('-p has to be > 0: EXIT.')
       
    #
    # FILL ARRAY WITH PARAMETER SETS TO PROCESS
    #
    # this array contains all jobs that have to be run
    aJobs = []
    # Clean the path variables up
    sHH = os.path.abspath(os.path.expanduser(oArgs.sHH.strip()))
    sDB = os.path.abspath(os.path.expanduser(oArgs.sDB.strip()))

    # create db string here, based on avaialbe databases, break if none available
    aDBs = glob.glob(sDB+'/*_hhm_db')
    if len(aDBs)<1:
        oParser.errer('At least one database needs to be available. EXIT.')
    sDB = "'%s'" %(' '.join(aDBs))

    # Non standard library for reading fasta-files
    try:
        from Bio import SeqIO as Parser  # BioPython import for reading seq files
    except ImportError:
        Parser = FastaParser()
        #oParser.error('BioPython needs to be installed. Exit.')

    # create a unique id
    sUUID = uuid.uuid4()
    iJob = 1
    # e.g. read tasks from file, here one fasta-entry one task
    for record in Parser.parse(oF, "fasta") :
        aJobs.append((iJob,
                      '>%s\n%s'%(record.description,record.seq),
                      sDB,
                      sHH,
                      sUUID))
        iJob+=1
    oF.close()

    # For timing
    fStart_time = timer() ## very crude

    iNumJobs = len(aJobs)
    sys.stderr.write('#JOBS TO RUN: %i | #CONCURRENT PROCESSES TO USE: %i\n' %(iNumJobs, iNofProcesses))
    
    # create pool of workers ---------------------
    pool = Pool(processes=iNofProcesses)

    if not oArgs.bNoProgress:
        #====================================================================
        # "chunksize" usually only makes a noticable performance
        # difference for very large iterables
        # Here, I set it to one to get the progress bar working nicly
        # Otherwise it will not give me the correct number of processes left
        # but chunksize number.
        chunksize = 1
        #====================================================================
        aResults = pool.map_async(my_func, aJobs, chunksize=chunksize)
    else: 
        aResults = pool.map_async(my_func, aJobs)
        
    # No more work to add to pool
    pool.close()

    # Number of total jobs
    iNumJobs = len(aJobs)

    if not oArgs.bNoProgress:
        # Progress bar
        #==============================
        # This can be changed to make progressbar bigger or smaller
        iProgressBarLength = 50
        #==============================
        while not aResults.ready():
            iNumNotDone = aResults._number_left
            iNumDone = iNumJobs-iNumNotDone
            iBarDone = iNumDone*iProgressBarLength/iNumJobs
            sBar = ('=' * iBarDone).ljust(iProgressBarLength)
            iPercent = int(iNumDone*100/iNumJobs)
            sys.stderr.write("[%s] %s%%\r" % (sBar, str(iPercent).rjust(2)))
            sys.stderr.flush()
            time.sleep(0.1)  # wait a bit: here we test all .1 secs
        # Finish the progress bar
        sBar = '=' * iProgressBarLength
        sys.stderr.write("[%s] 100%%\r"%(sBar))
        
    aResults = aResults.get()
    # --------------------------------------------
    
    fEnd_time = timer()
    if oArgs.bTIME: sys.stderr.write('\nRUNTIME(s): %.4f | AVG/JOB: %.4f\n' %(fEnd_time - fStart_time, (fEnd_time - fStart_time)/iNumJobs))
    
    # Do stuff with the results
    for res in aResults:
        args = res[0]
        aHHres = res[1]
        oFout.write('%s\n--\n'%(args[1]))
        for tRes in aHHres[0:oArgs.iTop]:
            oFout.write('%s\n%s\n'%(tRes[0], tRes[1]))
        oFout.write('#\n')

    oFout.close()
   
    return
        
if __name__ == '__main__':
    sys.exit(main())


