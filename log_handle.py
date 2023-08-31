# log handler using logging

import logging

logfile = None
logging.basicConfig(filename=logfile,
                    format='%(asctime)s - %(levelname)s : %(message)s',
                    level=logging.INFO)

if logfile is None:
    print('Start logging')
else:
    print('Start logging into', logfile)
