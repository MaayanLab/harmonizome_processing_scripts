import time
import urllib

# This script doesn't actually work.
# Broad Institute requires users to register before downloading the data.
# Registration is free.
# Must sign in and manually download the file.

folder = 'input/'
date = time.strftime('%Y%m%d', time.localtime())
fileurl = 'http://www.broadinstitute.org/achilles/datasets/5/download/Achilles_QC_v2.4.3.rnai.Gs.gct'
filename = 'dataset_{}_original.gct'.format(date)
urllib.urlretrieve(fileurl, folder + filename)
