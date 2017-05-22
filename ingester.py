import os
import sys
import tempfile

import GEOparse
import pymongo

ok = all([x in os.environ for x in ['MONGO_READ_URL', 'MONGO_WRITE_URL']] )
if not ok:
    print("MONGO_READ_URL and MONGO_WRITE_URL environment")
    print("variables must be set.")
    sys.exit(1)


def ingest(accession, disease):
    geodir = tempfile.TemporaryDirectory()
    # import IPython;IPython.embed()
    print("geodir is {}".format(geodir.name))
    gse = GEOparse.get_GEO(geo=accession, destdir=geodir.name)
    clinical_collection = "{}_geo".format(disease)
    metadata = gse.metadata
    write_client = pymongo.MongoClient(os.getenv('MONGO_WRITE_URL'))
    write_db = write_client.some_db # TODO FIXME
    write_coll = write_db[clinical_collection]
    write_coll.insert_one(metadata)

    import IPython;IPython.embed()

    return "something"

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("usage: python ingester.py <ACCESSION> <DISEASE>")
        sys.exit(0)
    result = ingest(sys.argv[1], sys.argv[2])
    print(result)
