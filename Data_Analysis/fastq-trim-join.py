import os
import glob

import gzip
import shutil


fastqDir ='Rawreads/' 

fastxDir = '00-postfastx/'
trimmomaticDir = '00-posttrimmomatic/'
fastqjoinDir = '00-postfastqjoin/'

if not os.path.exists(fastxDir):
    os.makedirs(fastxDir)

if not os.path.exists(trimmomaticDir):
	os.makedirs(trimmomaticDir)


if not os.path.exists(fastqjoinDir):
	os.makedirs(fastqjoinDir)

def getIDs(myFileList):
    id = []
    for f in myFileList:
        id.append(f.split('_L001')[0].split(fastqDir[1:])[-1])
    return(id)

def runFASTX(fwd,rev):
    fid = fwd[:fwd.rfind('fastq')]
    os.system('fastx_trimmer -z -f 5 -i %s -o %s' % (fwd,fid+'cropped.fastq.gz'))
    rid = rev[:rev.rfind('fastq')]
    os.system('fastx_trimmer -z -f 5 -l 250 -i %s -o %s' % (rev,rid+'cropped.fastq.gz'))
    print(fid,rid)


def runTrimmomatic(fwd,rev):

    fid = fwd[fwd.find('/')+1:fwd.rfind('.fastq')]
    rid = rev[rev.find('/')+1:rev.rfind('.fastq')]#[:rev.rfind('.fastq.gz')]
    pairedFWD =  trimmomaticDir + fid + '.paired.fastq.gz'
    unpairedFWD = trimmomaticDir + fid + '.unpaired.fastq.gz'
    pairedREV = trimmomaticDir + rid + '.paired.fastq.gz'
    unpairedREV = trimmomaticDir + rid +  '.unpaired.fastq.gz'
    cmd = 'trimmomatic PE %s %s %s %s %s %s SLIDINGWINDOW:100:28' % (fwd, rev,pairedFWD, unpairedFWD, pairedREV, unpairedREV)
    os.system(cmd)

def joinFastq(pairedFWD,pairedREV,h):
    output = fastqjoinDir + h + '.cropped.trimmed.joined.fastq'
    cmd = 'fastq-join -p 3 -m 20 %s %s -o %s' % (pairedFWD, pairedREV, output)
    os.system(cmd)


flist = glob.glob(fastqDir+'*R1*fastq.gz')
print(len(flist),fastqDir+'*R1*fastq.gz')

idList = getIDs(flist)

for myID in idList:
    print(myID)
    with gzip.open('%s%s_L001_R1_001.fastq.gz'%(fastqDir,myID), 'rb') as f_in:
        with open('%s%s_L001_R1_001.fastq'%(fastqDir,myID), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    with gzip.open('%s%s_L001_R2_001.fastq.gz'%(fastqDir,myID), 'rb') as f_in:
        with open('%s%s_L001_R2_001.fastq'%(fastqDir,myID), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f_in,f_out)
    runFASTX('%s%s_L001_R1_001.fastq'%(fastqDir,myID),'%s%s_L001_R2_001.fastq'%(fastqDir,myID))
    runTrimmomatic('%s%s_L001_R1_001.cropped.fastq.gz'%(fastqDir,myID),'%s%s_L001_R2_001.cropped.fastq.gz'%(fastqDir,myID))
    joinFastq('%s%s_L001_R1_001.cropped.paired.fastq.gz'%(trimmomaticDir,myID), '%s%s_L001_R2_001.cropped.paired.fastq.gz'% (trimmomaticDir,myID), myID)

