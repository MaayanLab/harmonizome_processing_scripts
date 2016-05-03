import time
import urllib
import zipfile
import os

folder = 'input/'
date = time.strftime('%Y%m%d', time.localtime())
fileurl = 'http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.135/BIOGRID-ORGANISM-3.4.135.mitab.zip'
zippedfilename = 'BIOGRID-ORGANISM-3.4.135.mitab.zip'
urllib.urlretrieve(fileurl, folder + zippedfilename)
with zipfile.ZipFile(folder + zippedfilename, mode='r') as zf:
	unzippedfilenames = zf.namelist()
	unzippedfilename = [x for x in unzippedfilenames if 'sapiens' in x][0]
	zf.extract(unzippedfilename, folder)
newfilename = 'dataset_{}_original.mitab'.format(date)
os.rename(folder + unzippedfilename, folder + newfilename)
os.remove(folder + zippedfilename)
