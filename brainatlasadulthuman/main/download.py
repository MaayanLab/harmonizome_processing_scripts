import time
import urllib
import zipfile
import os

date = time.strftime('%Y%m%d', time.localtime())
folder = 'input/'

filename_url = {'normalized_microarray_donor9861.zip':'http://human.brain-map.org/api/v2/well_known_file_download/178238387',
'normalized_microarray_donor10021.zip':'http://human.brain-map.org/api/v2/well_known_file_download/178238373',
'normalized_microarray_donor12876.zip':'http://human.brain-map.org/api/v2/well_known_file_download/178238359',
'normalized_microarray_donor14380.zip','http://human.brain-map.org/api/v2/well_known_file_download/178238316',
'normalized_microarray_donor15496.zip','http://human.brain-map.org/api/v2/well_known_file_download/178238266',
'normalized_microarray_donor15697.zip','http://human.brain-map.org/api/v2/well_known_file_download/178236545'}

for filename, fileurl in filename_url.iteritems():d
	urllib.urlretrieve(fileurl, folder + filename)
	if filename[-4:] == '.zip':
		with zipfile.ZipFile(folder + filename, mode='r') as zf:
			zf.extractall(folder)
		os.remove(folder + filename)
