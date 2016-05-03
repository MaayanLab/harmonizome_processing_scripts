import time
import urllib
import zipfile
import os

date = time.strftime('%Y%m%d', time.localtime())
folder = 'input/'

filename_url = {'NCI60_U133A_20070815.raw.csv.zip':'http://plugins.biogps.org/download/NCI60_U133A_20070815.raw.csv.zip',
'NCI60_sample_info.xls':'http://plugins.biogps.org/download/NCI60_sample_info.xls',
'GPL96-15653.txt':'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&amp;is_datatable=true&amp;acc=GPL96&amp;id=15653&amp;db=GeoDb_blob82'}

for filename, fileurl in filename_url.iteritems():
	urllib.urlretrieve(fileurl, folder + filename)
	if filename[-4:] == '.zip':
		with zipfile.ZipFile(folder + filename, mode='r') as zf:
			zf.extractall(folder)
		os.remove(folder + filename)
