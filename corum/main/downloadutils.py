# import re
import codecs
import json
# import requests
import urllib
import time

def cleanString(string, remstrings):
	for remstring in remstrings:
		string = string.replace(remstring, '')
	return string

def appendUniqueEntries(field_entry, field_uentry):
	for field in field_uentry:
		if type(field_entry[field]) is list:
			field_uentry[field] = field_uentry[field].union(field_entry[field])
		else:
			field_uentry[field].add(field_entry[field])
	return field_uentry

def writeFile(filename, field_entry_list, alsojson=True):
	fieldlist = field_entry_list[0].keys()
	fieldlist.sort()
	with codecs.open(filename, encoding='utf-8', mode='w') as fw:
		fw.write('\t'.join(fieldlist) + '\n')
		for field_entry in field_entry_list:
			writeLine(fw, fieldlist, field_entry)
	if alsojson:
		with codecs.open(filename.replace('.txt', '.json'), encoding='utf-8', mode='w') as fw:
			json.dump(field_entry_list, fw, indent=2, ensure_ascii=False, encoding='utf-8')
	
def writeLine(fw, fieldlist, field_entry):
	entries = []
	for field in fieldlist:
		entry = field_entry[field]
		if type(entry) is list:
			if type(entry[0]) is not str:
				entry = '|'.join([str(x) for x in entry])
			else:
				entry = '|'.join(entry)
		elif type(entry) is not str:
			entry = str(entry)
		entries.append(entry)
	fw.write('\t'.join(entries) + '\n')

def parseLine(line, field_idx, field_type, field_case, field_dlm, field_rem, dlm):
	entries = line.strip().split(dlm)
	field_entry = {}
	for field, idx in field_idx.iteritems():
		entry = entries[idx].strip()
		if field_rem[field] is not None:
			entry = cleanString(entry, field_rem[field])
		if field_type[field] == 'float':
			if field_dlm[field] is not None:
				entry = list(set([float(x.strip()) for x in entry.split(field_dlm[field]) if x.strip() != '']))
			else:
				entry = float(entry)
		elif field_type[field] == 'int':
			if field_dlm[field] is not None:
				entry = list(set([int(x.strip()) for x in entry.split(field_dlm[field]) if x.strip() != '']))
			else:
				entry = int(entry)
		else:
			if field_dlm[field] is not None:
				entry = list(set([changeCase(x.strip(), field_case[field]) for x in entry.split(field_dlm[field]) if x.strip() != '']))
			else:
				entry = changeCase(entry, field_case[field])
		field_entry[field] = entry
	return field_entry

def changeCase(string, case):
	if case == 'upper':
		string = string.upper()
	elif case == 'lower':
		string = string.lower()
	return string

def parseHeader(line, field_header, dlm):
	headers = line.strip().split(dlm)
	field_idx = {}
	for field, header in field_header.iteritems():
		field_idx[field] = headers.index(header)
	return field_idx

def downloadEndpoint(baseurl, endpoint, params, filename):
	# urllib.urlretrieve(baseurl + endpoint, filename, data=urllib.urlencode(params)) # post?
	urllib.urlretrieve('{0}{1}?{2}'.format(baseurl, endpoint, urllib.urlencode(params)), filename) # get?
	# with codecs.open(filename, mode='wb', encoding='utf-8') as fw:
	# 	data = requests.get(baseurl + endpoint, params=params, stream=True)
	# 	if not data.ok:
	# 		print('error downloading data')
	# 	else:
	# 		for block in data.iter_content(1024):
	# 			fw.write(block)

def entrezGeneIds2entrezGeneSymbols(entrezgeneids, entrezgeneid_entrezgenesymbol):
	return [{'geneid':x, 'genesymbol':entrezgeneid_entrezgenesymbol[x]} for x in entrezgeneids if x in entrezgeneid_entrezgenesymbol]
