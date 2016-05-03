import downloadutils as du
import time
import codecs

def main():
	folder = 'input/'
	date = time.strftime('%Y%m%d', time.localtime())
	baseurl = 'http://mips.helmholtz-muenchen.de/genre/proj/corum/'
	datasetendpoint = 'allComplexes.csv'
	datasetendpointparams = {}
	datasetprefix = 'dataset_{0}'.format(date)
	featureendpoint = 'complexdetails.html'
	featureendpointparams = {'id':''}
	featureendpointtemplate = 'complexdetails.html?id={}'
	species = 'Homo sapiens'
	description = ''
	featuretype = 'protein complex'
	datasetname = 'CORUM Protein Complexes'

	# download dataset
	du.downloadEndpoint(baseurl, datasetendpoint, datasetendpointparams, folder + '{0}_{1}'.format(datasetprefix, 'original.csv'))

	# parse dataset
	dlm = ';'
	field_header = {'id':'Complex id', 'name':'Complex name', 'organism':'organism', 'geneids':'subunits (Entrez IDs)'}
	field_type = {'id':'string', 'name':'string', 'organism':'string', 'geneids':'string', 'species':'string', 'description':'string', 'url':'string', 'enrichrname':'string'}
	field_case = {'id':'as is', 'name':'as is', 'organism':'lower', 'geneids':'as is', 'species':'as is', 'description':'as is', 'url':'as is', 'enrichrname':'as is'}
	field_dlm = {'id':None, 'name':None, 'organism':None, 'geneids':',', 'species':None, 'description':None, 'url':None, 'enrichrname':None}
	field_rem = {'id':None, 'name':None, 'organism':None, 'geneids':['(', ')'], 'species':None, 'description':None, 'url':None, 'enrichrname':None}
	field_uentry = {'id':set(), 'name':set(), 'geneids':set()}
	field_entry_list = []
	organisms = set()
	with codecs.open(folder + '{0}_{1}'.format(datasetprefix, 'original.csv'), encoding='utf-8', mode='r') as fr:
		field_idx = du.parseHeader(fr.readline(), field_header, dlm)
		for line in fr:
			field_entry = du.parseLine(line, field_idx, field_type, field_case, field_dlm, field_rem, dlm)
			organisms.add(field_entry['organism'])
			if field_entry['organism'] == 'human' and len(field_entry['geneids']) > 0 and field_entry['name'] != '':
				del field_entry['organism']
				field_entry['species'] = species
				field_entry['description'] = description
				field_entry['url'] = baseurl + featureendpointtemplate.format(field_entry['id'])
				field_entry['enrichrname'] = '{0}_{1}_{2}'.format(field_entry['name'], field_entry['species'], field_entry['id'])
				# print(field_entry)
				field_entry_list.append(field_entry)
				field_uentry = du.appendUniqueEntries(field_entry, field_uentry)

	# write parsed dataset
	du.writeFile(folder + '{0}_{1}'.format(datasetprefix, 'parsed.txt'), field_entry_list, alsojson=False)

	# print summary sets
	print('organism')
	print(organisms)
	for field, uentry in field_uentry.iteritems():
		print(field)
		print(uentry)

main()
