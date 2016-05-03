# files in input/enrichGene folder created by Yan Kou
# each file reports significant peaks from a ChIP-seq experiment
# peaks are reported as distance to nearest gene

import time
import zipfile
import codecs

folder = 'input/'
date = time.strftime('%Y%m%d', time.localtime())
zippedfilename = 'enrichGene.zip'
mergedfilename = 'dataset_{}_original.txt'.format(date)
maxdistance = 5000
uspecies = set()
profile_count = {}
with zipfile.ZipFile(folder + zippedfilename, mode='r') as zf, codecs.open(folder + mergedfilename, encoding='utf-8', mode='w') as fw:
	unzippedfilenames = zf.namelist()[2:]
	for chipseqprofile in unzippedfilenames:
		specs = chipseqprofile.replace('enrichGene/', '').split('_')
		tfORhm = specs[0].strip()
		cellORtissue = specs[1].strip()
		peaktype = specs[-2].strip()
		if specs[2].strip() == 'hg19':
			species = 'Homo sapiens'
		elif specs[2].strip() == 'mm9':
			species = 'Mus musculus'
		else:
			continue
		uspecies.add(species)
		if (tfORhm,cellORtissue,species) not in profile_count:
			profile_count[(tfORhm,cellORtissue,species)] = 1
		else:
			profile_count[(tfORhm,cellORtissue,species)] += 1
		genedistlist = []
		with zf.open(chipseqprofile, mode='r') as csp:
			for line in csp:
				entries = line.strip().split('\t')
				if len(entries) == 2:
					genesym = entries[0].strip()
					distance = entries[1].strip()
					if abs(float(distance)) < maxdistance:
						genedistlist += [genesym, distance]
		fw.write('\t'.join([tfORhm, cellORtissue, species, str(profile_count[(tfORhm,cellORtissue,species)]), peaktype] + genedistlist) + '\n')
print(uspecies)
