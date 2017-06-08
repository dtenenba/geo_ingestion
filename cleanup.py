
"""
Drop database in order to have a clean run.
This should ONLY be run on a development database where
you don't care if data is deleted. It DROPS
the database and REMOVES ALL DATA in it!
"""

import os
import sys

import pymongo

URL = os.getenv('MONGO_WRITE_URL')

if URL is None:
    print("MONGO_WRITE_URL is not set.")
    sys.exit(1)

CLIENT = pymongo.MongoClient(URL)
RESULT = CLIENT.tcga.command("dropDatabase")
print(RESULT)
