import time
import codecs
import re
import os
import urllib
import gzip

# this dataset is no longer maintained but data are still hosted on ftp site as of 4/4/16

def main():
	folder = 'input/'
	date = time.strftime('%Y%m%d', time.localtime())
	fileurl = 'http://cgap.nci.nih.gov/Pathways/BioCarta_Pathways'
	fileprefix = 'pathwaylist_{0}'.format(date)

	# download pathways
	urllib.urlretrieve(fileurl, folder + '{0}_{1}'.format(fileprefix, 'original.html'))

	# parse pathway name, id, url
	with codecs.open(folder + '{0}_{1}'.format(fileprefix, 'original.html'), encoding='utf-8', mode='r') as fr:
		hrefs = re.findall('<A class=genesrch href=.+</A>', fr.read())
	pathways = {}
	for href in hrefs:
		m = re.search('href="(?P<pathwayurl>\S+)"', href)
		pathwayurl = m.group('pathwayurl')
		pathwayid = pathwayurl.split('/')[-1].strip()
		m = re.search('>(?P<pathwayname>.+)</A>', href)
		pathwayname = m.group('pathwayname').strip()
		if '<IMG SRC=' not in pathwayname:
			pathways[pathwayid] = {'name':pathwayname, 'id':pathwayid, 'url':pathwayurl}
	with codecs.open(folder + '{0}_{1}'.format(fileprefix, 'parsed.txt'), encoding='utf-8', mode='w') as fw:
		fw.write('\t'.join(['pathwayname', 'pathwayid', 'pathwayurl']) + '\n')
		for pathwayid, pathway in pathways.iteritems():
			fw.write('\t'.join([pathway['name'], pathway['id'], pathway['url']]) + '\n')

	# parse old pathway name, id, url
	genelisturltemplate = 'ftp://ftp1.nci.nih.gov/pub/PID/molList/{}.mol.tab.gz'
	with codecs.open(folder + 'pathwaylist_20141217_original.html', encoding='utf-8', mode='r') as fr:
		hrefs = re.finditer('<a href="http://pid.nci.nih.gov/search/pathway_landing.shtml\?pathway_id=(?P<pathwayid>\d+)\S+">(?P<pathwayname>.+)</a>', fr.read())
	oldpathways = {}
	for href in hrefs:
		pathwayid = href.group('pathwayid').strip()
		pathwayname = href.group('pathwayname').strip()
		pathwaydataurl = genelisturltemplate.format(pathwayid)
		oldpathways[pathwayid] = {'name':pathwayname, 'id':pathwayid, 'dataurl':pathwaydataurl}
	
	# download pathway data
	os.mkdir('molfiles')
	for pathwayid, pathway in oldpathways.iteritems():
		urllib.urlretrieve(pathway['dataurl'], folder + 'molfiles/{}.mol.tab.gz'.format(pathwayid))

	# parse pathway data
	for pathwayid in oldpathways:
		with gzip.open(folder + 'molfiles/{}.mol.tab.gz'.format(pathwayid), mode='r') as fr:
			fr.readline()
			genes = re.finditer('\S+ {1}.+ {3}(?P<geneid>\d+)', fr.read().decode('utf-8'))
		geneids = set()
		for gene in genes:
			geneids.add(gene.group('geneid').strip())
		oldpathways[pathwayid]['geneids'] = list(geneids)
	with codecs.open(folder + 'dataset_20141217_original.txt', encoding='utf-8', mode='w') as fw:
		for pathwayid, pathway in oldpathways.iteritems():
			fw.write('\t'.join([pathway['name'], pathway['id'], pathway['dataurl']] + pathway['geneids']) + '\n')

main()
