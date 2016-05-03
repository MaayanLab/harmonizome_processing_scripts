import time
import urllib
import zipfile
import os

date = time.strftime('%Y%m%d', time.localtime())
folder = 'input/'

filename_url = {'human_sample_annot.csv':'http://plugins.biogps.org/download/human_sample_annot.csv',
'gnf1h-gcrma.zip':'http://plugins.biogps.org/download/gnf1h-gcrma.zip',
'gnf1h-gcrma-unaveraged.zip':'http://plugins.biogps.org/download/gnf1h-gcrma-unaveraged.zip',
'GPL96-15653.txt':'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&amp;is_datatable=true&amp;acc=GPL96&amp;id=15653&amp;db=GeoDb_blob82',
'gnf1h-anntable.zip':'http://plugins.biogps.org/download/gnf1h-anntable.zip'}

for filename, fileurl in filename_url.iteritems():
	urllib.urlretrieve(fileurl, folder + filename)
	if filename[-4:] == '.zip':
		with zipfile.ZipFile(folder + filename, mode='r') as zf:
			zf.extractall(folder)
		os.remove(folder + filename)
