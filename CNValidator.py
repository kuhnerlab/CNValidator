# CNValidator.py
# Copyright 2018 Mary K. Kuhner

#Permission is hereby granted, free of charge, to any person obtaining a
#copy of this software and associated documentation files (the "Software"),
#to deal in the Software without restriction, including without limitation
#the rights to use, copy, modify, merge, publish, distribute, sublicense,
#and/or sell copies of the Software, and to permit persons to whom the
#Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included
#in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
#OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
#ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#OTHER DEALINGS IN THE SOFTWARE.


min_useable_baf = 0.4    # if a position in the patient normal sample
                         # has a BAF less than this score, we will
                         # discard it

max_useable_baf = 0.65   # if a position in the patient normal sample
                         # has a BAF greater than this score, we will
                         # discard it

# BAF scores with values contained in "badbafvals" will be ignored
# when determing whether a segment "matches" another
badbafvals = ["NA","?"]

# segments which contain fewer than "min_bafcount" positions with valid
# BAF scores will be discarded
min_bafcount = 10

# at least this fraction of positions must "match" for two segments to be
# considered "matching"
min_matching = 0.95

# These chromosomes should not be used (generally because they are not
# diploid).  If your data, for example, calls the Y chromosome "chr24",
# you can add that name to this list.  (We do not exclude chr24 in the
# default version because in non-human data this could be an autosome.)

invalid_chromosomes = ["X","Y","chrX","chrY"]

# These strings give the column headers; if your input file has non-standard
# column headers, you can insert yours into the ends of these lists, but
# this will FAIL MISERABLY if you use a name that the program expects for, say,
# patient to indicate your sample.

patheaders = ["patient"]
sampleheaders = ["sample","biopsy"]
chromheaders = ["chr","chrom"]
startheaders = ["startpos","segstart"]
endheaders = ["endpos","segend"]
aheaders = ["intA"]
bheaders = ["intB"]

#######################################################################

def local_open(filename,mode):
  import gzip
  if filename.endswith("gz"):
    return gzip.open(filename,mode)

  return open(filename,mode)

def noquotes(mystring):
  return mystring.strip('"')

def readinputfile(infilename):
  rsegfiledict = {}
  rbaflist = []
  rnormfilename = None
  rtags = ["#SEGMENTATION","#NORMAL_BAF","#SAMPLE_BAF"]
  allrfiles = local_open(infilename,"r").readlines()
  if allrfiles[0].rstrip() not in rtags:
    print "File",infilename,"is not a proper input for validator.py"
    print "The file must begin with one of the three legal headers:"
    print "   #SEGMENTATION"
    print "   #NORMAL_BAF"
    print "   #SAMPLE_BAF"
    print "There should be a SEGMENTATION header for each segmentation"
    print "method being considered, but only one of each of the BAF headers"
    exit()
  import itertools
  iterator = allrfiles.__iter__()
  for line in iterator:
    line = line.rstrip()
    if line == "":  continue    # skip blank lines
    if line in rtags:
      doing = line
      if doing == rtags[0]:
        line = iterator.next()    
        line = line.rstrip()
        rsegfiledict[line[1:]] = []
        rsegref = rsegfiledict[line[1:]]
    else:
      if doing == rtags[0]:     # SEGMENTATION
        rsegref.append(line)
      if doing == rtags[1]:     # NORMAL_BAF
        rnormfilename = line
      if doing == rtags[2]:     # SAMPLE_BAF
        rbaflist.append(line)

  if len(rsegfiledict) == 0:
    print "No segmentation files found in",infilename
    exit()
  for method in rsegfiledict:
    if len(rsegfiledict) == 0:
      print "Method",method,"has no segmentation files in",infilename
      exit()
  if len(rbaflist) == 0:
    print "No sample BAF files found in",infilename
    exit()
  if rnormfilename is None:
    print "No normal BAF file found in",infilename
    exit()

  return rsegfiledict,rbaflist,rnormfilename

def makekeyforsegments(ksample,kchrom):
  return ksample + "_" + kchrom

def makekeyforbafs(kchrom,kstart,kend):
  return kchrom + str(kstart) + "_" + str(kend)

def insegment(newseg,chrom,origsegments):
  newstart = newseg[0]
  newend = newseg[1]
  for method in origsegments:
    for samplechromkey in origsegments[method]:
      mysample, mychrom = samplechromkey.split("_")
      if mychrom != chrom:
        continue
      samplesegs = origsegments[method][samplechromkey]
      for seg in samplesegs:
        segstart = seg[0]
        segend = seg[1]
        if segstart >= newend:
          break
        if segend <= newstart:
          continue
        if newstart >= segstart and newend <= segend:
          return True

  return False

def isstart(pos,chrom,origsegments):
  for method in origsegments:
    for samplechromkey in origsegments[method]:
      if ('_'+chrom) not in samplechromkey:
        continue
      samplesegs = origsegments[method][samplechromkey]
      for seg in samplesegs:
        segstart = seg[0]
        if pos == segstart:
          return True
  return False

def isend(pos,chrom,origsegments):
  for method in origsegments:
    for samplechromkey in origsegments[method]:
      if ('_'+chrom) not in samplechromkey:
        continue
      samplesegs = origsegments[method][samplechromkey]
      for seg in samplesegs:
        segend = seg[1]
        if pos == segend:
          return True
  return False

def findallelespecificcn(segs,segboundaries):
  searchmin = segboundaries[0]
  searchmax = segboundaries[1]
  for seg in segs:
    segend = seg[1]
    if segend <= searchmin:
      continue
    if segend >= searchmax:
      return [seg[2],seg[3]]

  return None

# this function assumes that the "filtering" BAFs are from the
# patient normal, and therefore there can only be one value!
def filterbadpositions(fbafs,fbmin,fbmax):
  fbads = {}
  for chrom in fbafs:
    chrombads = []
    for pospair in fbafs[chrom]:
      fpos = pospair[0]
      bafval = pospair[1].values()[0]
      if bafval < fbmin or bafval > fbmax:
        chrombads.append(fpos)
    fbads[chrom] = sorted(chrombads)
  return fbads

#adapted/copied from https://docs.python.org/2/library/bisect.html
# the "index" function
import bisect
def inbadpos(filterlist,searchpos):
  # locate the leftmost value exactly equal to searchpos
  ind = bisect.bisect_left(filterlist,searchpos)
  if ind != len(filterlist) and filterlist[ind] == searchpos:
    return True

  return False

def countmarkers(bafs):
  total = 0
  for chrom in bafs:
    for pospair in bafs[chrom]:
      total += 1
  return total

def readbafs(filename,filterthese):
  bafchrindex = None
  bafsbychrom = {}
  for line in local_open(filename,"r"):
    if line.startswith("##"):  continue
    if line == "":  continue    # skip blank lines
    line = line.rstrip().split("\t")
    line = [noquotes(x) for x in line]
    if bafchrindex is None:
      bafchrindex = line.index("Chr")
      bafposindex = line.index("Position")
      sids = line[3:]
      continue

    bafchr = line[bafchrindex]
    if bafchr in invalid_chromosomes:  # silently discard positions on
      continue                         # non-diploid chromosomes

    bafpos = int(line[bafposindex])

    # silently throw out markers whose name begins with "cnvi"
    marker_name = line[0]             # in good R fashion, the first column
    if marker_name.startswith("cnvi"):# has no header and we assume that it
      continue                        # is the marker's name...we hope so

    if filterthese:
      if inbadpos(filterthese[bafchr],bafpos):
        continue

    bafvals = line[3:]
    bafsbysid = {}
    for sid,bafval in zip(sids,bafvals):
      if bafval in badbafvals:
        bafsbysid[sid] = -1.0
      else:
        bafsbysid[sid] = float(bafval)
    if bafchr not in bafsbychrom:
      bafsbychrom[bafchr] = []
    bafsbychrom[bafchr].append([bafpos,bafsbysid])
    
  return bafsbychrom

def addmorebafs(origbafs,morebafs):
  for chrom in morebafs:
    if chrom not in origbafs:
      origbafs[chrom] = morebafs[chrom]
      continue
    newlist = []
    for newpair in morebafs[chrom]:
      positionnotfound = True
      for origpair in origbafs[chrom]:
        if newpair[0] == origpair[0]:
          positionnotfound = False
          for sid in newpair[1]:
            if sid in origpair[1]:
              if newpair[1][sid] != origpair[1][sid]:
                print "Found two BAF values for position",newpair[0]
                print newpair[1][sid],"and",origpair[1][sid]
                print "in sample",sid
                exit()
            origpair[1][sid] = newpair[1][sid]
      if positionnotfound:
        newlist.append(newpair)
    origbafs[chrom] += newlist

def makesidlist(bafsbychrom):
  sids = set()
  for chrom in bafsbychrom:
    for pospair in bafsbychrom[chrom]:
      for sid in pospair[1]:
        sids.add(sid)
  sidlist = list(sids)
  return sidlist

def sortchroms(chromset):
  numerics = []
  nonnumerics = []
  for chrom in chromset:
    if chrom.isdigit():  # would like to use isnumeric(), but we don't
                         # seem to have it...
      numerics.append(chrom)
    else:
      nonnumerics.append(chrom)
  numerics = sorted(numerics, key = lambda x: int(x))
  nonnumerics.sort()
  newchroms = numerics + nonnumerics
  return newchroms

def findbafsin(allbafs,boundaries,sid):
  fbafvals = []
  startpos = boundaries[0]
  endpos = boundaries[1]
  for pospair in allbafs:
    curpos = pospair[0]
    if curpos < startpos or curpos > endpos:
      continue
    fbafvals.append(pospair[1][sid])
  return fbafvals

def converttohilow(cbafs):
  hilows = []
  for baf in cbafs:
    if baf == -1.0:
      hilows.append("NA")
      continue
    if baf < 0.5:
      hilows.append("low")
      continue
    hilows.append("high")

  return hilows

def ismatching(mbafs1,mbafs2):
  hilow1 = converttohilow(mbafs1)
  hilow2 = converttohilow(mbafs2)
  matches = []
  for hl1,hl2 in zip(hilow1,hilow2):
    if hl1 == "NA" or hl2 == "NA":
      continue
    if hl1 == hl2:
      matches.append(1)
      continue
    matches.append(0)

  if not matches:  # since there were no BAFs in at least one of
    return False   # the inputs, we set "matching" arbitrarily now,
                   # and will filter on not enough BAFs for the
                   # "evalutation" field later

  pmatch = sum(matches)/float(len(matches))
  if pmatch <= 1.0 - min_matching or pmatch >= min_matching:
    return True

  return False

def howmanyuseable(pbafs):
  count = len(pbafs) - sum([n < 0 for n in pbafs])
  return count

def addscoredsegto(testdict,sid,chrom,pos,bafs,matches,scored):
  newkey = sid + "_" + chrom + "_" + str(pos)
  if newkey not in testdict:
    testdict[newkey] = {}
    # length for closed interval
    testdict[newkey]["length"] = scored["segend"] - scored["segstart"] + 1
    testdict[newkey]["matches"] = matches
    if scored["ascn"] is None:
      testdict[newkey]["aascn"] = -1
      testdict[newkey]["bascn"] = -1
    else:
      testdict[newkey]["aascn"] = scored["ascn"][0]
      testdict[newkey]["bascn"] = scored["ascn"][1]
    testdict[newkey]["bafs"] = bafs
  if matches:
    testdict[newkey]["matches"] = matches

def findindex(line,headers):
  found = False
  for item in headers:
    try:
      myindex = line.index(item)
      found = True
      break
    except:
      continue
  if not found:
    print "Could not find column header",headers,"in line",line
  return myindex

#######################################################################

import sys
import gzip

if len(sys.argv) != 2:
  print "USAGE: python CNValidator.py masterfile.txt"
  exit()

infilename = sys.argv[1]
segfiledict,baffilelist,normfilename = readinputfile(infilename)

dummy = {}
print "Reading normal BAFs from",normfilename
normalbafsbychrom = readbafs(normfilename,dummy)
badpos = filterbadpositions(normalbafsbychrom,min_useable_baf,max_useable_baf)
totalignore = 0
for chrom in badpos.iterkeys():
  totalignore += len(badpos[chrom])
totalmarkers = countmarkers(normalbafsbychrom)
#print
#print "A grand total of",str(totalignore),"out of",str(totalmarkers),
#print "markers will be ignored across all chromosomes"
print

print "Reading sample BAFs from",baffilelist[0]
samplebafsbychrom = readbafs(baffilelist[0],badpos)
for file in baffilelist[1:]:
  print "  ",file
  ubafs = readbafs(file,badpos)
  addmorebafs(samplebafsbychrom,ubafs)

sidlist = makesidlist(samplebafsbychrom)
for chrom in samplebafsbychrom:
  for pospair in samplebafsbychrom[chrom]:
    for sid in sidlist:
      if pospair[1][sid]:
        continue
      pospair[1][sid] = -1.0

segmethods = []
for method in segfiledict:
  segmethods.append(method)

check_pid = None
orig_patientfile = ""
allchroms = set()
allsamples = set()
allstarts = {}
allends = {}
origsegments = {}
for method in segmethods:
  segfilelist = segfiledict[method]
  origsegments[method] = {}
  osegsref = origsegments[method]
  for file in segfilelist:
    patindex = None
    check_sid = None
    for line in local_open(file,"r"):
      line = line.rstrip().split("\t")
      if line == "":  continue    # skip blank lines
      line = [noquotes(x) for x in line]
      if line[0] == "#":
        continue
  
      if patindex is None:
        patindex = findindex(line,patheaders)
        sampleindex = findindex(line,sampleheaders)
        chromindex = findindex(line,chromheaders)
        startindex = findindex(line,startheaders)
        endindex = findindex(line,endheaders)
        aindex = findindex(line,aheaders)
        bindex = findindex(line,bheaders)
        continue
  
      pid = line[patindex]
      sid = line[sampleindex]
      chrom = line[chromindex]
      try:
        startp = int(line[startindex])
      except ValueError:
        print "File",file,"has a bad starting position in line"
        print line
        exit()
      try:
        endp = int(line[endindex])
      except ValueError:
        print "File",file,"has a bad ending position in line"
        print line
        exit()
      try:
        aval = int(line[aindex])
      except ValueError:
        if line[aindex] == "NA":
          aval = -1
        else:
          print "File",file,"has a bad a-allele count in line"
          print line
          exit()
      try:
        bval = int(line[bindex])
      except ValueError:
        if line[bindex] == "NA":
          bval = -1
        else:
          print "File",file,"has a bad b-allele count in line"
          print line
  
      if check_pid is None:
        check_pid = pid
        orig_patientfile = file
      else:
        if check_pid != pid:
          print "found patient",check_pid,"in file",orig_patientfile,
          print "and patient",pid,"in file",file
          print "All samples must come from the same patient."
          exit()

      if check_sid is None:
        check_sid = sid
      else:
        if check_sid != sid:
          print "found sample",check_sid,"and sample",sid,"in file",file
          print "Each file must contain only one sample."
          exit()
  
      allchroms.add(chrom)
      allsamples.add(sid)
      newseg = [startp,endp,aval,bval]
      newkey = makekeyforsegments(sid,chrom)
      if newkey in osegsref:
        osegsref[newkey].append(newseg)
      else:
        osegsref[newkey] = [newseg,]
  
      if chrom in allstarts:
        allstarts[chrom].add(startp)
      else:
        allstarts[chrom] = set([startp,])
  
      if chrom in allends:
        allends[chrom].add(endp)
      else:
        allends[chrom] = set([endp,])
print "Processing data from",str(len(allchroms)),"chromosomes for",
print str(len(allsamples)),"samples"
print

# we create the canonical list of mini-segments and a new dictionary of
# lists for allstarts and allends which now contain the starts and ends
# of the minisegments.
#
# allchroms is also now sorted (conceptually) into 2 groups,
# first numerics in ascending order, then all non-numerics in alpha order
# and made into a dictionary of lists.
#
# we do rely on "sorted" returning a list
allchroms = sortchroms(allchroms)
minisegs = {}
newstarts = {}
newends = {}
for chrom in allchroms:
  allpos = list(allstarts[chrom]) + list(allends[chrom])
  allpos = sorted(set(allpos))  # make a unique sorted list
  newstarts[chrom] = set()
  newends[chrom] = set()
  # use indices to move startpos and endpos down the list
  for i1 in xrange(len(allpos) - 1):
    i2 = i1 + 1
    pos1 = allpos[i1]
    pos2 = allpos[i2]
    if insegment([pos1,pos2],chrom,origsegments):
      if isend(pos1,chrom,origsegments):
        pos1 -= 1
      if isstart(pos2,chrom,origsegments):
        pos2 += 1
      newstarts[chrom].add(pos1)
      newends[chrom].add(pos2)
      if chrom in minisegs:
        minisegs[chrom].append([pos1,pos2])
      else:
        minisegs[chrom] = [[pos1,pos2],]
  newstarts[chrom] = sorted(newstarts[chrom])
  newends[chrom] = sorted(newends[chrom])

allstarts = newstarts
allends = newends

scoredsegments = {}
for method in segmethods:
  scoredsegments[method] = {}
  ssegsref = scoredsegments[method]
  for sample in allsamples:
    for chrom in allchroms:
      samplekey = makekeyforsegments(sample,chrom)
      ssegsref[samplekey] = []
      samplesegs = origsegments[method][samplekey]
      chrombafs = samplebafsbychrom[chrom]
      miniseg = minisegs[chrom]
      for mseg in miniseg:
        cnscores = findallelespecificcn(samplesegs,mseg)
        bafvals = findbafsin(chrombafs,mseg,sample)
        newseg = {}
        newseg["segstart"] = mseg[0]
        newseg["segend"] = mseg[1]
        newseg["ascn"] = cnscores
        newseg["bafs"] = bafvals
        ssegsref[samplekey].append(newseg)
      

import itertools
allsamplepairs = []
for pair in itertools.combinations(allsamples,2):
  allsamplepairs.append(pair)

testsegs = {}
for method in segmethods:
  testsegs[method] = {}
  tsegref = testsegs[method]
  for sid1,sid2 in allsamplepairs:
    for chrom in allchroms:
      skey1 = makekeyforsegments(sid1,chrom)
      skey2 = makekeyforsegments(sid2,chrom)
      segs1 = scoredsegments[method][skey1]
      segs2 = scoredsegments[method][skey2]
      for seg1,seg2 in zip(segs1,segs2):
        bafs1 = seg1["bafs"]
        bafs2 = seg2["bafs"]
        if len(bafs1) != len(bafs2):
          print "Samples",sid1,"and",sid2,"have a differing number of BAF"
          print "scores for the positions on chromosome",chrom,"from",
          print str(seg1["segstart"]),"to",str(seg1["segend"]),"."
          print "Sample",sid1,"has",len(bafs1),"positions and sample",sid2,
          print "has",len(bafs2),"positions."
          exit()
        matches = ismatching(bafs1,bafs2)
        addscoredsegto(tsegref,sid1,chrom,seg1["segstart"],bafs1,matches,seg1)
        addscoredsegto(tsegref,sid2,chrom,seg2["segstart"],bafs2,matches,seg2)
  
notenoughbafs = {}
unknownallelecounts = []
outsegs = {}
for method in segmethods:
  outsegs[method] = {}
  osegref = outsegs[method]
  notenoughbafs[method] = {}
  for sample in allsamples:
    if sample not in osegref:
      osegref[sample] = {}
    notenoughbafs[method][sample] = []
    for chrom in allchroms:
      if chrom not in osegref[sample]:
        osegref[sample][chrom] = {}
      skey = makekeyforsegments(sample,chrom)
      seglist = scoredsegments[method][skey]
      for pos in allstarts[chrom]:
        testkey = sample + "_" + chrom + "_" + str(pos)
        testseg = testsegs[method][testkey]
        invalid = False
        if howmanyuseable(testseg["bafs"]) < min_bafcount:
          notenoughbafs[method][sample].append(testkey)
          invalid = True
        balanced = (testseg["aascn"] == testseg["bascn"])
        if testseg["matches"]:
          condition = "Match"
        else:
          condition = "No match"
        osegref[sample][chrom][pos] = {}
        osegref[sample][chrom][pos]["balanced"] = balanced
        osegref[sample][chrom][pos]["condition"] = condition
        osegref[sample][chrom][pos]["length"] = testseg["length"]
        osegref[sample][chrom][pos]["nbafs"] = len(testseg["bafs"])
  
        if testseg["aascn"] == -1 or testseg["bascn"] == -1:
          unknownallelecounts.append(testkey)
          invalid = True
          testcall = "UN/UN"
        else:
          testcall = str(testseg["aascn"]) + "/" + str(testseg["bascn"])
          if testseg["aascn"] > testseg["bascn"]:
            testcall = str(testseg["bascn"]) + "/" + str(testseg["aascn"])
        osegref[sample][chrom][pos]["call"] = testcall
  
        if invalid:
          success = "UN"
        else:
          if not balanced and condition == "Match":
            success = "TP"
          elif not balanced and condition == "No match":
            success = "FP"
          elif condition == "Match": # we must be balanced
            success = "FN"
          else:
            success = "TN"
        osegref[sample][chrom][pos]["evaluation"] = success

#useful debug code not needed in production
#for method in segmethods:
#  for sample in allsamples:
#    print method,"discarded",len(notenoughbafs),"segments in sample",
#    print sample
#    print "\tfor having fewer than",str(min_bafcount),
#    print "positions with good enough BAFs."
#print

# segments which lack a "condition" field with a score of "match"
# for any of the samples, need to have their "evaluation" field changed
# to "UN", for unvalidateable.
unvalidatable = {}
totalsegs = {}
for method in segmethods:
  unvalidatable[method] = 0
  totalsegs[method] = 0
  for chrom in allchroms:
    for pos in allstarts[chrom]:
      foundmatch = False
      for sample in allsamples:
        totalsegs[method] += 1
        oseg = outsegs[method][sample][chrom][pos]
        if oseg["condition"] == "Match":
          foundmatch = True
          break
      if foundmatch:
        continue
      for sample in allsamples:
        outsegs[method][sample][chrom][pos]["evaluation"] = "UN"
      unvalidatable[method] += 1

# useful debug code not needed in production
#for method in segmethods:
#  print method,"found and discarded",str(unvalidatable[method]),
#  print "segment/sample pairs out of",str(totalsegs[method])
#  print "\tfor failing to have any matches"
#print

possible_coverage = {}
achieved_coverage = {}
discarded_coverage = {}
validated_coverage = {}
contradicted_coverage = {}
segcount = {}
segvalidated = {}
segcontradicted = {}
for method in segmethods:
  possible_coverage[method] = {}
  achieved_coverage[method] = {}
  discarded_coverage[method] = {}
  validated_coverage[method] = {}
  contradicted_coverage[method] = {}
  segcount[method] = {}
  segvalidated[method] = {}
  segcontradicted[method] = {}
  for sample in allsamples:
    possible_coverage[method][sample] = {}
    achieved_coverage[method][sample] = {}
    discarded_coverage[method][sample] = {}
    validated_coverage[method][sample] = {}
    contradicted_coverage[method][sample] = {}
    segcount[method][sample] = {}
    segvalidated[method][sample] = {}
    segcontradicted[method][sample] = {}
    for chrom in allchroms:
      possible_coverage[method][sample][chrom] = 0.0
      achieved_coverage[method][sample][chrom] = 0.0
      discarded_coverage[method][sample][chrom] = 0.0
      validated_coverage[method][sample][chrom] = 0.0
      contradicted_coverage[method][sample][chrom] = 0.0
      segcount[method][sample][chrom] = 0
      segvalidated[method][sample][chrom] = 0
      segcontradicted[method][sample][chrom] = 0
      for pos in allstarts[chrom]:
        seg = outsegs[method][sample][chrom][pos]
        length = seg["length"]
        possible_coverage[method][sample][chrom] += length
        segcount[method][sample][chrom] += 1
        if seg["evaluation"] == "UN":
          discarded_coverage[method][sample][chrom] += length
        else:
          achieved_coverage[method][sample][chrom] += length
          if (not seg["balanced"]) and seg["condition"] == "Match":
            validated_coverage[method][sample][chrom] += length
            segvalidated[method][sample][chrom] += 1
            continue
          if seg["balanced"] and seg["condition"] == "No match":
            validated_coverage[method][sample][chrom] += length
            segvalidated[method][sample][chrom] += 1
            continue
          contradicted_coverage[method][sample][chrom] += length
          segcontradicted[method][sample][chrom] += 1

overallname = check_pid+"_overall_output.tsv"
overallout = local_open(overallname,"w")
outline = "Sample\tMethod\tMB_total\tMB_validated\tMB_contradicted\t"
outline += "MB_accuracy\tSegments_total\tSegments_validated\t"
outline += "Segments_contradicted\tSegments_accuracy\n"
overallout.write(outline)
for method in segmethods:
  for sample in allsamples:
    outline = sample + "\t" + method + "\t"
    mbtotal = 0.0
    validated = 0.0
    contradicted = 0.0
    nsegs = 0
    nsegsvalid = 0
    nsegscontradict = 0
    for chrom in allchroms:
      mbtotal += possible_coverage[method][sample][chrom]
      validated += validated_coverage[method][sample][chrom]
      contradicted += contradicted_coverage[method][sample][chrom]
      nsegs += segcount[method][sample][chrom]
      nsegsvalid += segvalidated[method][sample][chrom]
      nsegscontradict += segcontradicted[method][sample][chrom]
    mbtotal = mbtotal/1000000
    outline += ("%.6f" % mbtotal) + "\t"
    validated = validated/1000000
    outline += ("%.6f" % validated) + "\t"
    contradicted = contradicted/1000000
    outline += ("%.6f" % contradicted) + "\t"
    if validated + contradicted != 0.0:
      accuracy = validated / (validated + contradicted)
      outline += ("%.6f" % accuracy) + "\t"
    else:
      outline += "NA\t"
    outline += str(nsegs) + "\t" + str(nsegsvalid) + "\t"
    outline += str(nsegscontradict) + "\t"
    if nsegsvalid + nsegscontradict != 0:
      accuracy = float(nsegsvalid) / float(nsegsvalid + nsegscontradict)
      outline += ("%.2f" % accuracy) + "\n"
    else:
      outline += "NA\n"
    overallout.write(outline)
      
print "Writing overall report to",overallname
overallout.close()

for sample in allsamples:
  detailname = sample+"_detailed_output.tsv"
  detailout = local_open(detailname,"w")
  print "Writing detailed report for sample",sample,"to",detailname
  outline = "# Detailed report for sample "  + sample
  outline += "\n"
  detailout.write(outline)
  outline = "chrom\tstartpos\tendpos\tnBAFs"
  for method in segmethods:
    outline += "\t" + method + "_call\t" + method + "_evaluation"
  outline += "\n"
  detailout.write(outline)
  for chrom in allchroms:
    for pos in allstarts[chrom]:
      for method in segmethods:
        seg = outsegs[method][sample][chrom][pos]
        if method == segmethods[0]:
          outline = chrom + "\t"
          outline += str(pos) + "\t" + str(pos+seg["length"]-1) + "\t"
          outline += str(seg["nbafs"])
        outline += "\t" + seg["call"] + "\t"
        outline += seg["evaluation"]
      outline += "\n"
      detailout.write(outline)
  detailout.close()
