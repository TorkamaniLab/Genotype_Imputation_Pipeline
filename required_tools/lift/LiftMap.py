#!/usr/bin/python
import sys, os

def usage():
    print("%s -m old.map -p old.ped -d old.dat -o new.prefix -c chain.file.chain" % sys.argv[0] )
    print("lift old.map (optionally .ped and .dat) to new.prefix.map (.ped/.dat if possible)")

def die(msg):
    print msg
    sys.exit(2)
    
def myopen(fn):
    import gzip
    try:
        h = gzip.open(fn)
        ln = h.read(2) # read arbitrary bytes so check if @param fn is a gzipped file
    except:
        # cannot read in gzip format
        return open(fn)
    h.close()
    return gzip.open(fn)

def map2bed(fin, fout):
    fo = open(fout, 'w')
    for ln in myopen(fin):
        chrom, rs, mdist, pos = ln.split()
	chrom = 'chr' + chrom
        pos = int(pos)
        fo.write('%s\t%d\t%d\t%s\n' % (chrom, pos-1, pos, rs))
    fo.close()
    return True

# global var:
LIFTED_SET = set()
UNLIFTED_SET = set()
def liftBed(fin, fout, funlifted, chainFile):
    params = dict()
    params['LIFTOVER_BIN'] = 'liftOver'
    params['OLD'] = fin
    params['CHAIN'] = chainFile
    params['NEW'] = fout
    params['UNLIFTED'] = fout + '.unlifted'
    from string import Template
    cmd = Template('$LIFTOVER_BIN $OLD $CHAIN $NEW $UNLIFTED')
    cmd = cmd.substitute(params)
    os.system(cmd)
    #record lifted/unliftd rs
    for ln in myopen(params['UNLIFTED']):
	if len(ln) == 0 or ln[0] == '#':continue
	UNLIFTED_SET.add(ln.strip().split()[-1])
    for ln in myopen(params['NEW']):
	if len(ln) == 0 or ln[0] == '#':continue
	LIFTED_SET.add(ln.strip().split()[-1])
    
    return True

def bed2map(fin, fout):
    fo = open(fout, 'w')
    for ln in myopen(fin):
	chrom, pos0, pos1, rs = ln.split()
	chrom = chrom.replace('chr', '')
	fo.write('%s\t%s\t0.0\t%s\n' % (chrom, rs, pos1))
    fo.close()
    return True

def liftDat(fin, fout):
    fo = open(fout, 'w')
    for ln in myopen(fin):
	if len(ln) == 0 or ln[0] != 'M':
	    fo.write(ln)
	else:
	    t, rs = ln.strip().split()
	    if rs in LIFTED_SET:
		fo.write(ln)
    fo.close()
    return True
    
def liftPed(fin, fout, fOldMap):
    # two ways to do it:
    # 1. write unlifted snp list
    #    use PLINK to do this job using --exclude
    # 2. alternatively, we can write our own method
    # we will use method 2
    marker = [i.strip().split()[1] for i in open(fOldMap)]
    flag = map(lambda x: x not in UNLIFTED_SET, marker)
    # print marker[:10]
    # print flag[:10]
    fo = open(fout, 'w')
    for ln in myopen(fin):
        f = ln.strip().split()
        l = len(f)
        f = f[:6] + [ f[i*2] + ' '+f[i*2 +1] for i in xrange(3, l/2 )]
        fo.write('\t'.join(f[:6]))
        fo.write('\t')
        if len(f[6:]) != len(flag):
            die('Inconsistent length of ped and map files')
        newMarker = [m for i, m in enumerate(f[6:]) if flag[i]]
        fo.write('\t'.join(newMarker))
        fo.write('\n')
        #print marker[:10]      
        #die('test')
    return True

def makesure(result, succ_msg, fail_msg = "ERROR"):
    if result:
        print 'SUCC: ', succ_msg
    else:
        print 'FAIL: ', fail_msg
        sys.exit(2)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-m', dest='mapFile', required = True)
    parser.add_argument('-p', dest='pedFile')
    parser.add_argument('-d', dest='datFile')
    parser.add_argument('-c', dest='chainFile', required = True)
    parser.add_argument('-o', dest='prefix', required = True)
    args = parser.parse_args()
    
    oldBed = args.mapFile + '.bed'
    makesure(map2bed(args.mapFile, oldBed),
             'map->bed succ')

    newBed = args.prefix + '.bed'
    unlifted = args.prefix + '.unlifted'
    makesure(liftBed(oldBed, newBed, unlifted, args.chainFile),
             'liftBed succ')

    newMap = args.prefix + '.map'
    makesure(bed2map(newBed, newMap),
	     'bed->map succ')
    
    if args.datFile:
	newDat = args.prefix + '.dat'
	makesure(liftDat(args.datFile, newDat),
		 'liftDat succ')

    if args.pedFile:
	newPed = args.prefix + '.ped'
	makesure(liftPed(args.pedFile, newPed, args.mapFile),
		 'liftPed succ')
