
"""
Drop database in order to have a clean run.
"""

import os
import sys

import pymongo

URL = os.getenv('MONGO_WRITE_URL')

if URL is None:
    print("MONGO_WRITE_URL is not set.")
    sys.exit(1)

CLIENT = pymongo.MongoClient(URL)
# FIXME database name will change at some point
RESULT = CLIENT.some_db.command("dropDatabase")
print(RESULT)
