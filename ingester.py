import os
import sys
import tempfile
import statistics

import GEOparse
import pymongo
import pandas as pd

ok = all([x in os.environ for x in ['MONGO_READ_URL', 'MONGO_WRITE_URL']] )
if not ok:
    print("MONGO_READ_URL and MONGO_WRITE_URL environment")
    print("variables must be set.")
    sys.exit(1)


def ingest(accession, disease):
    geodir = tempfile.TemporaryDirectory()
    # import IPython;IPython.embed()
    print("geodir is {}".format(geodir.name))
    print("Downloading data set {} from GEO....".format(accession),
          flush=True)
    gse = GEOparse.get_GEO(geo=accession, destdir=geodir.name, silent=True)
    clinical_collection_name = "{}_geo".format(disease)
    metadata = gse.metadata
    write_client = pymongo.MongoClient(os.getenv('MONGO_WRITE_URL'))
    write_db = write_client.some_db # TODO FIXME
    clinical_collection = write_db[clinical_collection_name]
    clinical_collection.insert_one(metadata)
    # FIXME this doesn't seem right:
    print("Inserting metadata into clinical collection (1)...",
          flush=True)
    for key, value in gse.gsms.items():
        clinical_collection.insert_one(value.metadata)
    mol_collection_name = "{}_geo_{}".format(disease, accession)
    mol_collection = write_db[mol_collection_name]
    key0 = list(gse.gsms.keys())[0]
    numrows = len(gse.gsms[key0].table)
    docs = []
    print("Inserting into {}....".format(mol_collection_name), flush=True)
    for idx in range(numrows):
        if idx % 1000 == 0:
            print("Inserted {} of {} rows...      \r".format(idx, numrows),
                  end="", flush=True)
        data = {}
        _min = _max = None
        for samplename, sample in gse.gsms.items():
            tbl = gse.gsms[samplename].table
            probe, value = tbl.iloc[idx, :]
            data[samplename] = value
            _min = value if _min is None else min(_min, value)
            _max = value if _max is None else max(_max, value)
        doc = dict(id=probe, data=data, min=_min, max=_max)
        # mol_collection.insert_one(doc)
        docs.append(doc)
    print("                              \r", flush=True)
    print("Done.", flush=True)
    print("here we go")
    mol_collection.insert_many(docs)
    print("really done")

    # for now, assume there is just one platform (GPL) for each data set (GSE)
    gplkey = list(gse.gpls.keys())[0]
    gpl = gse.gpls[gplkey]

    gpl_collection_name = "geo_{}".format(gplkey)
    if not gpl_collection_name in write_db.collection_names():
        docs = []
        gpl_collection = write_db[gpl_collection_name]
        print("inserting into {}....".format(gpl_collection_name),flush=True)
        numrows = len(gpl.table)
        for idx in range(numrows):
            if idx % 1000 == 0:
                print("Inserted {} of {} rows...      \r".format(idx, numrows),
                      end="", flush=True)
            row = gpl.table.iloc[idx, :].to_dict()
            row['id'] = row['ID']
            del row['ID']
            row['symbol'] = row['Gene Symbol']
            del row['Gene Symbol']
            row['entrez'] = row['ENTREZ_GENE_ID']
            del row['ENTREZ_GENE_ID']
            # gpl_collection.insert_one(row)
            docs.append(row)

    print("                              \r", flush=True)
    print("Done.", flush=True)
    print("here we go")
    gpl_collection.insert_many(docs)
    print("really done")


    # map probe IDs to gene symbol, and reduce duplicates to only the
    # most variable gene:

    # make a new data frame that has all the sample data in it
    # (with probe IDs as row names and sample ids as column names)

    # get the row names:
    firstkey = list(gse.gsms.keys())[0]
    firstsample = gse.gsms[firstkey]
    rownames = firstsample.table['ID_REF'].values.tolist()
    colnames = gse.gsms.keys()
    # create empty data frame with appropriate row & col names
    gse_df = pd.DataFrame(index=rownames, columns=colnames)
    # populate data frame one column at a time
    for samplekey, sample in gse.gsms.items():
        column = sample.table['VALUE'].tolist()
        gse_df[samplekey] = column

    # add a hugo column to gse_df
    hugo = gpl.table['Gene Symbol']

    # but first, clean up the gene symbol names by splitting
    # into segments delimited by whitespace and keeping only the first
    # segment.
    hugo = [str(x).split()[0] for x in hugo]

    # now add the column:
    gse_df['hugo']  = hugo

    # find duplicates:
    dupes_series = gse_df.duplicated('hugo') # by default gets unique ones
    which = dupes_series[dupes_series] # only those that are True
    dupes = which.index.tolist()

    # define a function to call once for each duplicate:
    def find_largest_variance(gene_symbol):
        probeset = gse_df[gse_df.hugo == gene_symbol]
        ncol = len(gse_df.columns)
        flip = probeset.iloc[:, 0:ncol-1].T
        probe_var = pd.Series(index=flip.columns)
        for idx, colname in enumerate(flip.columns):
            col = flip[colname]
            probe_var[idx] = statistics.variance(col)
        # return name of largest member of probe_var
        return probe_var.idxmax()


    keep_probe = pd.Series(index=dupes)

    import IPython;IPython.embed()


    for dupe in dupes:
        keep_probe[dupe] = find_largest_variance(dupe)




    return "something"

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("usage: python ingester.py <ACCESSION> <DISEASE>")
        sys.exit(0)
    result = ingest(sys.argv[1], sys.argv[2])
    print(result)
