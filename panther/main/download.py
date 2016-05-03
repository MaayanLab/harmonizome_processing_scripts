import time
import urllib

folder = 'input/'
date = time.strftime('%Y%m%d', time.localtime())
fileurl = 'ftp://ftp.pantherdb.org//pathway/current_release/SequenceAssociationPathway3.4.txt'
filename = 'dataset_{}_original.txt'.format(date)
urllib.urlretrieve(fileurl, folder + filename)
