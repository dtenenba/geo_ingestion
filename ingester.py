"""
The GEO datasets we are initially looking to ingest are:
GSE12102, GSE16102, GSE20196, GSE21050, GSE21122, GSE23980, GSE30929, GSE648
These all go with the 'sarc' disease (sarcoma).
"""

import os
import sys
import tempfile
import statistics
import argparse


import GEOparse
import pymongo
import pandas as pd

OK_ = all([x in os.environ for x in['MONGO_READ_URL', 'MONGO_WRITE_URL']])
if not OK_:
    print("MONGO_READ_URL and MONGO_WRITE_URL environment")
    print("variables must be set.")
    sys.exit(1)


def get_data_frame(gse):
    """
    make a new data frame that has all the sample data in it
    (with probe IDs as row names and sample ids as column names)
    """
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
        print("samplekey is {}, len is {}, columns are {}".format(samplekey, len(column), sample.table.columns))
        gse_df[samplekey] = column
    return gse_df


def delete_or_exit(db_, collection_name, force):
    """Take the appropriate action when a collection already exists."""
    if collection_name in db_.collection_names():
        if force:
            print("-f option specified, dropping '{}' collection...")
            db_.drop_collection(collection_name)
        else:
            print("The collection '{}' already exists.".format(collection_name))
            print("To overwrite, use the -f option. Exiting.")
            sys.exit(1)


def ingest(accession, disease, force): # pylint: disable=too-many-locals
    """
    Ingest data into MongoDb.
    Args:
        accession: GEO accession number
        disease: Oncoscape disease name
        force: Whether to overwrite existing collections, or exit
    """
    geodir = tempfile.TemporaryDirectory()
    print("geodir is {}".format(geodir.name))
    print("Downloading data set {} from GEO....".format(accession),
          flush=True)
    # silent=True has no effect,
    # see https://github.com/guma44/GEOparse/issues/19
    gse = GEOparse.get_GEO(geo=accession, destdir=geodir.name, silent=True)
    clinical_collection_name = "{}_geo".format(disease)
    metadata = gse.metadata
    write_client = pymongo.MongoClient(os.getenv('MONGO_WRITE_URL'))
    write_db = write_client.some_db # TODO FIXME
    # FIXME: shouldn't this have the accession # in there somewhere? :
    clinical_collection = write_db[clinical_collection_name]
    delete_or_exit(write_db, clinical_collection_name, force)
    clinical_collection.insert_one(metadata)
    # FIXME this doesn't seem right:
    print("Inserting metadata into clinical collection (1)...",
          flush=True)
    for value in gse.gsms.values():
        clinical_collection.insert_one(value.metadata)
    mol_collection_name = "{}_geo_{}".format(disease, accession)
    delete_or_exit(write_db, mol_collection_name, force)
    mol_collection = write_db[mol_collection_name]
    key0 = list(gse.gsms.keys())[0]
    numrows = len(gse.gsms[key0].table)
    docs = []
    print("Inserting into {}....".format(mol_collection_name), flush=True)
    gse_df = get_data_frame(gse)
    for idx in range(len(gse_df.index)):
    # for idx in range(numrows):
        if idx % 1000 == 0:
            print("Processed {} of {} rows...      \r".format(idx, numrows),
                  end="", flush=True)
        row = gse_df.iloc[idx, :]
        data = {}
        probe = gse_df.index[idx]
        for rowidx, cell in enumerate(row):
            data[gse_df.index[rowidx]] = float(cell)
        doc = dict(id=probe, data=data, min=float(min(row)), max=float(max(row)))
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
        print("inserting into {}....".format(gpl_collection_name), flush=True)
        numrows = len(gpl.table)
        for idx in range(numrows):
            if idx % 1000 == 0:
                print("Processed {} of {} rows...      \r".format(idx, numrows),
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

    # add a hugo column to gse_df
    hugo = gpl.table['Gene Symbol']

    # but first, clean up the gene symbol names by splitting
    # into segments delimited by whitespace and keeping only the first
    # segment.
    hugo = [str(x).split()[0] for x in hugo]

    # now add the column:
    gse_df['hugo'] = hugo

    # find duplicates:
    dupes_series = gse_df.duplicated('hugo')
    which = dupes_series[dupes_series] # only those that are True
    temp = which.index.tolist()
    # get unique values (but make it a list, not a set):
    dupes = list(set(gse_df.loc[temp, "hugo"].tolist()))

    # define a function to call once for each duplicate:
    def find_largest_variance(gene_symbol):
        """Find largest variance"""
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


    # this seems to take a while
    print("Finding largest variance for each of {} dupes...".format(len(dupes)),
          flush=True)
    for idx, dupe in enumerate(dupes):
        if idx % 100 == 0:
            print("{}           \r".format(idx), end="", flush="")
        keep_probe[dupe] = find_largest_variance(dupe)
    print("Done")

    diff = list(set(hugo) - set(dupes))
    which = gse_df.hugo.isin(diff)
    gene_gse = gse_df[which]
    # append is like rbind in R:
    gene_gse = gene_gse.append(gse_df.loc[keep_probe, :])
    gene_gse.index = gene_gse.hugo
    del gene_gse['hugo'] # remove hugo column

    mol_collection_name = "geo_{}_{}_gene".format(disease, accession)
    delete_or_exit(write_db, mol_collection_name, force)

    mol_collection = write_db[mol_collection_name]

    # TODO now insert data from gene_gse into mol_collection
    docs = []

    print("Inserting gene-oriented molecular collection {}...".format(mol_collection_name),
          flush=True)
    for idx in range(len(gene_gse.index)):
        if idx % 1000 == 0:
            print("{} of {} rows processed...\r".format(idx, len(gene_gse.index)),
                  end="", flush=True)
        gene_symbol = gene_gse.index[idx]
        row = gene_gse.iloc[idx, :]
        data = {}
        for rowidx, cell in enumerate(row):
            data[gene_gse.index[rowidx]] = float(cell)
        doc = dict(id=gene_symbol, data=data, min=float(min(row)),
                   max=float(max(row)))
        docs.append(doc)

    print("Done processing.", flush=True)
    print("Now inserting one big glob...")
    mol_collection.insert_many(docs)
    print("Done.")

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description="ingest GEO data")
    PARSER.add_argument('accession', metavar='ACCESSION', type=str,
                        help='GEO accession number')
    PARSER.add_argument('disease', metavar='DISEASE', type=str,
                        help='disease name')
    PARSER.add_argument('-f', help='overwrite existing collections',
                        action='store_true')

    ARGS = PARSER.parse_args()
    ingest(ARGS.accession, ARGS.disease, ARGS.f)
