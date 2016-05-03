import re
import codecs
import json
import requests
import time

def main():
	folder = 'input/'
	date = time.strftime('%Y%m%d', time.localtime())
	baseurl = 'http://humancyc.org'
	classinstancesendpoint = '/HUMAN/class-instances'
	classinstancesparams = {'object':'Pathways'}
	pathwaygenesendpoint = '/HUMAN/pathway-genes'
	pathwaygenesparams = {'object':''}
	species = 'Homo sapiens'

	# get pathways
	pathwaylisthtml = requests.get(baseurl + classinstancesendpoint, params=classinstancesparams).text
	hrefs = re.finditer('HREF="(?P<pathwayendpoint>\S+)">(?P<pathwayname>.+)</A>', pathwaylisthtml)
	pathways = {}
	for href in hrefs:
		pathwayendpoint = href.group('pathwayendpoint').strip()
		pathwayurl = baseurl + pathwayendpoint
		pathwayid = pathwayendpoint.split('=')[-1].strip()
		pathwayname = re.sub('</.>', '', re.sub('<.>', '', href.group('pathwayname'))).strip()
		enrichrname = '_'.join([pathwayname, species, pathwayid])
		pathways[pathwayid] = {'enrichrname':enrichrname, 'name':pathwayname, 'id':pathwayid, 'url':pathwayurl, 'species':species, 'genes':[]}

	# get gene lists
	for pathwayid in pathways:
		pathwaygenesparams['object'] = pathwayid
		genelist = requests.get(baseurl + pathwaygenesendpoint, params=pathwaygenesparams).text.strip().split('\n')
		del genelist[:3]
		for row in genelist:
			entries = row.strip().split('\t')
			genesymbol = entries[2].strip()
			# hcgeneid = entries[0]
			# genes.append({'humancycgeneid':hcgeneid, 'symbol':genesymbol})
			pathways[pathwayid]['genes'].append(genesymbol)
		time.sleep(0.2)

	# write pathway gene lists
	with codecs.open(folder + 'dataset_{}_original.gmt'.format(date), encoding='utf-8', mode='w') as fw:
		for pathwayid, pathway in pathways.iteritems():
			fw.write('\t'.join([pathway['enrichrname'], pathway['name'], pathway['id'], pathway['species'], pathway['url']]+pathway['genes']) + '\n')
	with codecs.open(folder + 'dataset_{}_original.json'.format(date), encoding='utf-8', mode='w') as fw:
		json.dump(pathways, fw, indent=2, ensure_ascii=False)

main()
