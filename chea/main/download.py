import shutil

folder = 'input/'
shutil.copyfile(folder + 'Enrichr_ChEA_20141013.gmt', folder + 'dataset_20141013_original.gmt')

# import urllib
# import time

# folder = 'input/'
# date = time.strftime('%Y%m%d', time.localtime())
# fileurl = 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=ChEA_2015'
# filename = 'datset_{}_original.gmt'.format(date)
# urllib.urlretrieve(fileurl, folder + filename)
