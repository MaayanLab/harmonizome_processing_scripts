import time
import urllib
import zipfile
import os

folder = 'input/'
date = time.strftime('%Y%m%d', time.localtime())
fileurl = 'http://www.reactome.org/download/current/ReactomePathways.gmt.zip'
zippedfilename = 'ReactomePathways.gmt.zip'
urllib.urlretrieve(fileurl, folder + zippedfilename)
with zipfile.ZipFile(folder + zippedfilename, mode='r') as zf:
	unzippedfilename = zf.namelist()[0]
	zf.extract(unzippedfilename, folder)
newfilename = 'dataset_{}_original.gmt'.format(date)
os.rename(folder + unzippedfilename, folder + newfilename)
os.remove(folder + zippedfilename)
