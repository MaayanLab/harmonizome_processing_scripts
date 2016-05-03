import time
import urllib
import zipfile
import os

date = time.strftime('%Y%m%d', time.localtime())
folder = 'input/'

filename_url = {'gnf1m-anntable.zip':'http://plugins.biogps.org/download/gnf1m-anntable.zip',
'GNF1M_geneatlas_20120817.zip':'http://plugins.biogps.org/download/GNF1M_geneatlas_20120817.zip',
'geneatlas_MOE430_20090327.raw.avg.csv.zip':'http://plugins.biogps.org/download/geneatlas_MOE430_20090327.raw.avg.csv.zip',
'GPL1261-14790.txt':'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&amp;is_datatable=true&amp;acc=GPL1261&amp;id=14790&amp;db=GeoDb_blob82'}

for filename, fileurl in filename_url.iteritems():
	urllib.urlretrieve(fileurl, folder + filename)
	if filename[-4:] == '.zip':
		with zipfile.ZipFile(folder + filename, mode='r') as zf:
			zf.extractall(folder)
		os.remove(folder + filename)
