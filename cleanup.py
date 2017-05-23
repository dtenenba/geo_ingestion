
"""
Drop database in order to have a clean run.
"""

import os
import sys

import pymongo

url = os.getenv('MONGO_WRITE_URL')

if url is None:
    print("MONGO_WRITE_URL is not set.")
    sys.exit(1)

client = pymongo.MongoClient(url)
# FIXME database name will change at some point
result = client.some_db.command("dropDatabase")
print(result)
