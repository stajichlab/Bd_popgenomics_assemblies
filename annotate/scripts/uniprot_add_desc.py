#!/usr/bin/env python3
import  urllib.request, csv, sys, warnings,re
import codecs, time
urlbase="https://www.uniprot.org/uniprot/?query={}&sort=score&columns=id,entry%20name,reviewed,organism,genes(PREFERRED),genes,protein%20names&format=tab&limit=1"
if len(sys.argv) < 2:
    print("No file argument")
    sys.exit(2)

file = sys.argv[1]
pat = re.compile(r'^(.+?)\s+\(')
with open(file,'r') as genes:
    dat = csv.reader(genes,delimiter="\t")
    for line in dat:
        time.sleep(5)
        if line[0].startswith("#"):
            continue
        url=urlbase.format(line[0])
#        warnings.warn("url is {}".format(url))
        with urllib.request.urlopen(url) as uniprot:
            uniprotdat = csv.reader(codecs.iterdecode(uniprot,'utf-8'),delimiter="\t")
            header = 0
            for row in uniprotdat:
                if header == 0:
                    header = 1
                    continue
                newdesc = row[6]
                m = pat.match(newdesc)
                if m:
                    newdesc = m.group(1)

                print("\t".join([line[0],newdesc]))
