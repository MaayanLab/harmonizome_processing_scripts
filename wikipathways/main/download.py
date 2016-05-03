import time
import urllib

folder = 'input/'
date = time.strftime('%Y%m%d', time.localtime())
fileurl = 'http://www.pathvisio.org/data/bots/gmt/wikipathways.gmt'
filename = 'dataset_{}_original.gmt'.format(date)
urllib.urlretrieve(fileurl, folder + filename)
