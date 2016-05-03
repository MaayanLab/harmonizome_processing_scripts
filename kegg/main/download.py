import requests
import codecs
import json
import time

def main():
	folder = 'input/'
	date = time.strftime('%Y%m%d', time.localtime())
	baseurl = 'http://rest.kegg.jp/'
	species = 'Homo sapiens'
	pathwayurltemplate = 'http://www.genome.jp/kegg/pathway/hsa/{}.html'

	endpoint = 'list/pathway/hsa'
	params = {}
	doc = requests.get(baseurl+endpoint, params=params).text
	lines = doc.strip().split('\n')
	pathways = {}
	for line in lines:
		entries = line.strip().split('\t')
		pathwayid = entries[0].replace('path:', '').strip()
		pathwayname = entries[1].replace(' - Homo sapiens (human)', '').strip()
		pathwayurl = pathwayurltemplate.format(pathwayid).strip()
		enrichrname = '{0}_{1}_{2}'.format(pathwayname, species, pathwayid)
		pathways[pathwayid] = {'name':pathwayname, 'id':pathwayid, 'enrichrname':enrichrname, 'species':species, 'url':pathwayurl, 'geneids':[]}

	endpoint = 'conv/hsa/ncbi-geneid'
	params = {}
	doc = requests.get(baseurl+endpoint, params=params).text
	lines = doc.strip().split('\n')
	genes = {}
	for line in lines:
		entries = line.strip().split('\t')
		hsaid = entries[1].replace('hsa:', '').strip()
		geneid = entries[0].replace('ncbi-geneid:', '').strip()
		genes[hsaid] = {'hsaid':hsaid, 'geneid':geneid}

	endpoint = 'link/hsa/pathway'
	params = {}
	doc = requests.get(baseurl+endpoint, params=params).text
	lines = doc.strip().split('\n')
	for line in lines:
		entries = line.strip().split('\t')
		pathwayid = entries[0].replace('path:', '').strip()
		hsaid = entries[1].replace('hsa:', '').strip()
		pathways[pathwayid]['geneids'].append(genes[hsaid]['geneid'])

	with codecs.open(folder + 'dataset_{}_original.gmt'.format(date), encoding='utf-8', mode='w') as fw:
		for pathwayid, pathway in pathways.iteritems():
			fw.write('\t'.join([pathway['enrichrname'], pathway['name'], pathway['id'], pathway['species'], pathway['url']]+pathway['geneids']) + '\n')

	with codecs.open(folder + 'dataset_{}_original.json'.format(date), encoding='utf-8', mode='w') as fw:
		json.dump(pathways, fw, indent=2, ensure_ascii=False)

main()
