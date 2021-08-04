import os
import csv
import gzip
import pysam
from Bio import SeqIO
import re
import numpy as np
import pandas as pd
import yaml
from collections import OrderedDict
import logging
import subprocess
import itertools
import time
from multiprocessing import Pool
import shutil

import sys
csv.field_size_limit(sys.maxsize)

logger = logging.getLogger(__name__)

def callShell(cmd, noException=False, noPipefail=False):
    defaults = ['set', '-o', 'nounset;', 'set', '-o', 'errexit;']
    if not noPipefail:
        defaults.extend(['set', '-o', 'pipefail;'])

    command = ' '.join(defaults + cmd)
    call(command, shell=True, noException=noException)

def call(cmd, shell=False, executable=None, noException=False,
         stdin=None, stdout=None, stderr=None):
    if type(cmd) is list:
        logger.info('Running command: ' + ' '.join(cmd))
    else:
        logger.info('Running command: %s' % cmd)

    if shell is True and executable is None:
        executable = '/bin/bash'

    retc = subprocess.call(cmd, shell=shell, executable=executable,
                           stdin=stdin, stdout=stdout, stderr=stderr)
    if not noException and retc != 0:
        raise RuntimeError('Failed to run command: %s' % cmd)

    if retc != 0:
        if type(cmd) is list:
            logger.info('Command: %s FAILED with exit code %s' %
                        (' '.join(cmd), retc))
        else:
            logger.info('Command: %s FAILED with exit code %s' % (cmd, retc))

    return retc

def runEnvCmd(condaEnv, cmd):
    # run bash cmd in a conda env, no need to deactivate
    fullcmd = ['set', '+eu', '&&',
               '.', '$(conda info --base)/etc/profile.d/conda.sh', '&&',
               'conda', 'activate', condaEnv, '&&',
               'set', '-eu', '&&']
    fullcmd.extend(cmd)
    call(' '.join(fullcmd), shell=True)

def _unpackConfigPath(params, rootDir):
    for key in params:
        value = params[key]
        if type(value) is dict:
            _unpackConfigPath(value, rootDir)
        else:
            if re.match('^\/', str(value)): # with absolute path
                pass
            elif re.match('.+\/', str(value)): # need to expand
                params[key] = os.path.join(rootDir, value)
            else:
                pass

def getConfigSection(configPath, rootDir, section):
    tmpParams = OrderedDict()
    with open(configPath) as configFile:
        tmpParams = yaml.load(configFile, yaml.Loader)
    _unpackConfigPath(tmpParams, rootDir)
    paramSection = tmpParams[section]
    return paramSection

def emptyDumpster(list2delete):
    """a list can contain files and folders"""
    for item in list2delete:
        if os.path.isfile(item):
            os.remove(item)
        elif os.path.isdir(item):
            shutil.rmtree(item)
        else:
            logger.info('Unable to delete: %s', item)

def linkFile(a, b):
    """ force to link a file to b file"""
    if os.path.isfile(b):
        os.remove(b)
    os.link(a, b)

def runParlFuncs(numJobs, applyFunc, applyArgs=None, initFunc=None,
                callbackFunc=None, numProcesses=1):
    """
    Uses multiprocessing to launch multiple instances of applyFunc which is defined at top level
    """

    pool = Pool(numProcesses, initFunc)
    try:
        results = []
        for i in range(numJobs):
            args = None
            if applyArgs:
                args = applyArgs[i]
            result = pool.apply_async(applyFunc, [args], callback=callbackFunc)
            results.append(result)

        while len(results) > 0:
            temp = results
            time.sleep(2)
            for r in temp:
                if r.ready():
                    r.get()
                    results.remove(r)

        pool.close()
        pool.join()

    except KeyboardInterrupt:
        pool.terminate()
        pool.join()
        raise

def runParlFuncsReturn(numJobs, applyFunc, applyArgs=None, initFunc=None,
                callbackFunc=None, numProcesses=1):
    """
    Uses multiprocessing to launch multiple instances of applyFunc which is defined at top level
    return run values in a list
    """
    pool = Pool(numProcesses, initFunc)
    try:
        results = []
        returnList = []
        for i in range(numJobs):
            args = None
            if applyArgs:
                args = applyArgs[i]
            result = pool.apply_async(applyFunc, [args], callback=callbackFunc)
            results.append(result)

        while len(results) > 0:
            temp = results
            time.sleep(2)
            for r in temp:
                if r.ready():
                    rval = r.get()
                    returnList.append(rval)
                    results.remove(r)

        pool.close()
        pool.join()

    except KeyboardInterrupt:
        pool.terminate()
        pool.join()
        raise

    return(returnList)

def runParlCmds(cmds, outList=list(), errList=list(), numProcesses=None):
    # input cmd should be a simple cmd without pipe
    def waitProcesses(processes):
        for p in processes:
            if p.wait() != 0:
                raise RuntimeError(
                    'Child process exited with error code %s' %
                    p.returncode)

    if len(outList) > 0 and len(outList) != len(cmds):
        raise RuntimeError('len(outList) != len(cmds)')
    if len(errList) > 0 and len(errList) != len(cmds):
        raise RuntimeError('len(errList) != len(cmds)')

    if numProcesses is None:
        numProcesses = len(cmds)

    processes = []
    filesToClose = []
    for cmd, out, err in itertools.zip_longest(cmds, outList, errList,
                                                fillvalue='-'):
        stdout = None
        stderr = None
        if isinstance(out, str) and out != '-':
            stdout = open(out, 'w')
            filesToClose.append(stdout)
        if isinstance(err, str) and err != '-':
            stderr = open(err, 'w')
            filesToClose.append(stderr)

        logger.info('Running CMD ...: %s' % ' '.join(cmd))
        logger.info('stdout: %s' % out)
        logger.info('stderr: %s' % err)
        p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, shell=False, bufsize=-1)
        processes.append(p)

        # Loop until one (or more) precess(es) terminate
        while len(processes) == numProcesses:
            for p in processes:
                p.poll()
                if p.returncode is None:
                    continue
                if p.returncode != 0:
                    raise RuntimeError(
                        'Child process exited with error code %s' %
                        p.returncode)
                processes.remove(p)
            time.sleep(0.1)

    # All processes were launched, now just wait for them to terminate
    waitProcesses(processes)
    # close all files
    [file.close() for file in filesToClose]
    

def file2dict(filename, keyname=None, delimiter='\t'):
    """function for making a dictionary with feature dictionaries
    from a text file with a header
    """
    dname = OrderedDict()
    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=delimiter)
        headers = reader.fieldnames
        l = 0
        for row in reader:
            l += 1
            if keyname is None:
                dname[l] = row
            else:
                dname[row[keyname]] = row
    return [dname, headers]

def tfile2dict(filename, fieldnames, keyname=None, delimiter='\t'):
    """function for making a dictionary
    from a text file without a header
    """
    dname = OrderedDict()
    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile, fieldnames=fieldnames, delimiter=delimiter)
        l = 0
        for row in reader:
            l += 1
            if keyname is None:
                dname[l] = row
            else:
                dname[row[keyname]] = row
    return dname

def dict2tsv(dictD, tsvFile, dumpKey="Id"):
    dictFeatures = list(next(iter(dictD.values())).keys())
    outfh = open(tsvFile, 'w')
    if dumpKey=="False":
        outfh.write('\t'.join([feature for feature in dictFeatures]) + '\n')
        for i in dictD:
            outfh.write('\t'.join(str(dictD[i][feature]) for feature in dictFeatures) + '\n')
    else:
        outfh.write('\t'.join([dumpKey, '\t'.join([feature for feature in dictFeatures])]) + '\n')
        for i in dictD:
            outfh.write('\t'.join([str(i), '\t'.join(str(dictD[i][feature]) for feature in dictFeatures)]) + '\n')
    outfh.close()

def delDictFeatures(dictD, fList):
    dictD2 = dictD
    for i in dictD:
        for f in fList:
            if f in dictD2[i]:
                del dictD2[i][f]
    return dictD2

def fillDict(dictD, fill=None):
    # some dict has missing key and val pairs:
    # get all keys
    dictD2 = dictD
    fList = list(next(iter(dictD.values())).keys())
    for i in dictD:
        for f in fList:
            if f not in dictD[i]:
                if fill is None:
                    dictD2[i][f] = np.nan
                else:
                    dictD2[i][f] = fill 
    return dictD2 
    
def simpleD2tsv(dictD, tsvFile):
    dictFeatures = list(dictD.keys())
    outfh = open(tsvFile, 'w')
    outfh.write('\t'.join(dictFeatures) + '\n')
    outfh.write('\t'.join([str(dictD[feature]) for feature in dictFeatures]) + '\n')
    outfh.close()

def list2tsv(inList, tsvFile):
    with open(tsvFile, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for i in inList:
            writer.writerows([[i]])

def list2freqCount(inList):
    # report aggregated counts
    listName = '|'.join(list((pd.Series(inList).value_counts()).index))
    listCount = '|'.join([str(i) for i in list(pd.Series(inList).value_counts())])
    freqList = ':'.join([listName, listCount])
    return(freqList)

def splitBed(workBed, splitBedNum, outDir):
    # split work bed file into chuncks while keeping same chromosome in the same file
    workBedD = tfile2dict(filename=workBed, fieldnames=["chr", "start", "end"])
    splitline = int(round(len(workBedD)/int(splitBedNum)))
    def unique(sequence): # get uniq sequence without reorder
        seen = set()
        return [x for x in sequence if not (x in seen or seen.add(x))]
    x = [workBedD[j]["chr"] for j in [i for i in range(0, len(workBedD), splitline)][1:]+[len(workBedD)]]
    lastChrs = unique(x) # last chr for each chunck
    x0 = [workBedD[line]["chr"] for line in workBedD]
    allChrs = unique(x0)
    chunkfileD = OrderedDict() # make ordered dict
    startChr = allChrs[0]
    f = 0
    for lastChr in lastChrs:
        f += 1
        chunkBedFile = os.path.join(outDir, "chunk"+str(f)+"_work.bed")
        chunkfileD[chunkBedFile] = OrderedDict()
        chunkfileD[chunkBedFile]["chrs"] = allChrs[allChrs.index(startChr):allChrs.index(lastChr)+1]
        if allChrs.index(lastChr)+1 < len(allChrs):
            startChr = allChrs[allChrs.index(lastChr)+1]

    # open file for write
    chunkBedFileList = []
    for chunkBedFile in chunkfileD:
        chunkBedFileList.append(chunkBedFile)
        outfh = open(chunkBedFile, 'w')
        for line in workBedD:
            if workBedD[line]["chr"] in chunkfileD[chunkBedFile]["chrs"]:
                outfh.write("\t".join([workBedD[line]["chr"], workBedD[line]["start"], workBedD[line]["end"]]) + "\n")
        outfh.close()
    return chunkBedFileList


def getReadlenfromFastq(fastqfile):
    # find read length from the first read
    if fastqfile.endswith('.gz'):
        fastq = gzip.open(fastqfile, 'r')
    else:
        fastq = open(fastqfile, 'r')
    i = 0
    for line in fastq:
        i += 1
        if i == 2:
            readlen = len(line.strip())
            break
    return readlen

def getMinReadlenfromFastq(fastqfile, n=100):
    # find min read lengh of first 100 reads
    if fastqfile.endswith('.gz'):
        fastq = gzip.open(fastqfile, 'rt')
    else:
        fastq = open(fastqfile, 'rt')
    i = 0
    minReadlen = 100000
    for line in fastq:
        i += 1
        if re.search('^\+$|^@', line):
            pass
        else:
            readlen = len(line.strip())
            if readlen < minReadlen:
                minReadlen = readlen
        if i > n*4:
            break
    return minReadlen

def getReadcountfromFastq(fastqfile):
    if fastqfile.endswith('.gz'):
        fastq = gzip.open(fastqfile, 'rt')
    else:
        fastq = open(fastqfile, 'rt')
    readcount = 0
    for line in fastq:
        line = line.strip()
        if re.match("^\+$", line):
            readcount += 1
    return readcount        

def clipFqEnds(fqIn, fqOut, clip5, clip3):
    if fqIn.endswith('.gz'):
        fastq = gzip.open(fqIn, 'r')
    else:
        fastq = open(fqIn, 'r')
    with open(fqOut, 'w') as outH:
        c = 0
        for line in fastq:
            line = line.strip()
            c += 1
            if c == 2: # sequence
                if clip3 == 0:
                    line = line[clip5:len(line)]
                else:
                    line = line[clip5:-clip3]
                line = line.upper()    
            elif c == 4:
                if clip3 == 0:
                    line = line[clip5:len(line)]
                else:
                    line = line[clip5:-clip3]
                c = 0
            else:
                pass
            outH.write(line + '\n')    

def getReadlenfromBam(bamfile):
    samfile = pysam.AlignmentFile(bamfile, "rb")
    for read in samfile:
        readLen = len(read.get_forward_sequence())
        break
        """
        readlen = 0
        warning = 0
        for s in read.cigartuples:
            # 0:match, 1:insertion, 4:softclip
            if s[0] in (0, 1, 4):
                readlen += s[1]
            if s[0] in (5, 6, 7, 8):
                warning = 1
        if warning == 0:
            break
        """
    return readLen

def ab12fqW(ab1File, fqFile, newId):
    records = SeqIO.parse(ab1File, "abi")
    for record in records:
        record.id = newId
        count = SeqIO.write(record, fqFile, "fastq")
        logger.info("Converted %i records" % count)

def getOneSeqfromFa(faFile, oneFaId, outFile):
    records = SeqIO.parse(faFile, "fasta")
    oneSeq = None
    for record in records:
        if record.id == oneFaId:
            oneSeq = str(record.seq)
            SeqIO.write(record, outFile, "fasta")
    return oneSeq        

def getSeqsfromFa(faFile, faList, outFile):
    records = SeqIO.parse(faFile, "fasta")
    foundList = []
    with open(outFile, 'w') as outH:
        for record in records:
            if record.id in faList:
                SeqIO.write(record, outH, "fasta")
                foundList.append(record.id) 
    missingSeqs = list(set(faList) - set(foundList))
    if len(missingSeqs) > 0:
        logger.info("unfound seqs:"+ "|".join(missingSeqs))

def refGb2faW(refId):
     gbFile = refId + '.gb'
     dummy = SeqIO.read(gbFile, "gb")
     faFile = refId + '.fa'
     SeqIO.write(dummy, faFile, "fasta")

def resetQual(fqIn, fqOut, reset = 0):
    with open(fqIn, 'r') as inH, open(fqOut, 'w') as outH:
        c = 0
        for line in inH:
            c += 1
            if c == 2: # sequence
                iter = re.finditer("N|n", line)
                indices = [m.start(0) for m in iter] # 0 based
            elif c == 4: # qscorea
                line = list(line)
                for i in indices:
                    line[i] = chr(reset+33) # assign 0 to N/n 
                line = ''.join(line)
            outH.write(line)   

def readExl(exlFile):
    #exlFile = "/mnt/r/Sequences/WorkSpace/Projects/Deliverables/PACT133_PT0461-M485-UCLA/PT0461-M485_for_oligo_ordering_with_Y84C_BlpI_framework_columnbased.xlsx"
    xls = pd.ExcelFile(exlFile)
    sheet_to_df_map = OrderedDict()
    for sheet_name in xls.sheet_names:
        if re.search("_P\d_Sense", sheet_name) or re.search("_P\d_S\.", sheet_name):
            sheet_to_df_map[sheet_name] = xls.parse(sheet_name)
    df = pd.concat(sheet_to_df_map.values())
    df = df.reset_index(drop=True) # reset index
    # convert df into dict
    outD = df.to_dict('index')
    return outD

def runPileup(args):
    # a function run in a subprocess
    # args stored in a dict
    logger.info('Running ...')
    samtoolsExe = args["samtoolsExe"]
    hg19Fa = args["hg19Fa"]
    workBed = args["workBed"]
    bam = args["bam"]
    pileup = args["pileup"]
    cmd = [samtoolsExe, "mpileup", "-d", "5000", "-xBQ0", "-f", hg19Fa, "-l", workBed, bam, "|",
           "bgzip", "-c", ">", pileup, ";",
           "tabix -s 1 -b 2 -e 2", pileup]
    mytools.callShell(cmd)
    logger.info('Finished.')
