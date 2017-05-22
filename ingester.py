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
    clinical_collection_name = "{}_geo".format(disease)
    metadata = gse.metadata
    write_client = pymongo.MongoClient(os.getenv('MONGO_WRITE_URL'))
    write_db = write_client.some_db # TODO FIXME
    clinical_collection = write_db[clinical_collection_name]
    clinical_collection.insert_one(metadata)
    # FIXME this doesn't seem right:
    for key, value in gse.gsms.items():
        clinical_collection.insert_one(value.metadata)
    mol_collection_name = "geo_{}_{}".format(disease, accession)
    mol_collection = write_db[mol_collection_name]

    import IPython;IPython.embed()

    return "something"

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("usage: python ingester.py <ACCESSION> <DISEASE>")
        sys.exit(0)
    result = ingest(sys.argv[1], sys.argv[2])
    print(result)
