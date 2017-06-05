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



class Wrapper(object):
    '''
    Object wrapper class.
    From https://code.activestate.com/recipes/577555-object-wrapper-class/
    This a wrapper for objects. It is initialized with the object to wrap
    and then proxies the unhandled getattribute methods to it.
    Other classes are to inherit from it.
    '''
    def __init__(self, obj, gpl_idx, disease):
        '''
        Wrapper constructor.
        @param obj: object to wrap
        '''
        # wrap the object
        self._wrapped_obj = obj
        self.gpl_idx = gpl_idx
        self.disease = disease
        gpls = list(self._wrapped_obj.gpls.keys())
        self.gpl = gpls[self.gpl_idx]

    def gsms(self):
        """Gets each sample associated with a specific platform."""
        hsh = {}
        for key, value in self._wrapped_obj.gsms.items():
            if value.metadata['platform_id'][0] == self.gpl:
                hsh[key] = value
        return hsh

    def get_molecular_collection_name(self, probe_centric=True):
        """Gets the molecular collection name"""
        accession = self._wrapped_obj.metadata['geo_accession'][0]
        if len(self._wrapped_obj.gpls) == 1:
            if probe_centric:
                return "{}_geo_{}".format(self.disease, accession)
            else:
                return "{}_geo_{}_gene".format(self.disease, accession)
        else:
            if probe_centric:
                return "{}_geo_{}-{}".format(self.disease, accession, self.gpl)
            else:
                return "{}_geo_{}-{}_gene".format(self.disease, accession, self.gpl)

    def __getattr__(self, attr):
        # see if this object has attr
        # NOTE do not use hasattr, it goes into
        # infinite recurrsion
        if attr in self.__dict__:
            # this object has it
            return getattr(self, attr)
        # proxy to the wrapped object
        return getattr(self._wrapped_obj, attr)

class GEOSeries: # pylint: disable=too-few-public-methods
    """
    Wraps a GSE object and allows you to iterate over each
    set of samples associated with a specific platform.
    """

    def __init__(self, gse, disease):
        self.current = 0
        self.gse = gse
        self.disease = disease

    def __next__(self):
        if self.current >= len(self.gse.gpls):
            raise StopIteration
        else:
            self.current += 1
            return Wrapper(self.gse, self.current - 1, self.disease)

    @property
    def series(self):
        """Return the number of series (different platforms)"""
        return len(self.gse.gpls)

    def __iter__(self):
        return self

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.gse, attr)



def get_data_frame(gse):
    """
    make a new data frame that has all the sample data in it
    (with probe IDs as row names and sample ids as column names)
    """
    # get the row names:
    firstkey = list(gse.gsms().keys())[0]
    firstsample = gse.gsms()[firstkey]
    rownames = firstsample.table['ID_REF'].values.tolist()
    colnames = gse.gsms().keys()
    # create empty data frame with appropriate row & col names
    gse_df = pd.DataFrame(index=rownames, columns=colnames)
    # populate data frame one column at a time
    for samplekey, sample in gse.gsms().items():
        column = sample.table['VALUE'].tolist()
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

def get_from_geo(accession, disease):
    """Downloads a dataset from GEO."""
    geodir = tempfile.TemporaryDirectory()
    print("geodir is {}".format(geodir.name))
    print("Downloading data set {} from GEO....".format(accession),
          flush=True)
    # silent=True has no effect,
    # see https://github.com/guma44/GEOparse/issues/19
    raw_gse = GEOparse.get_GEO(geo=accession, destdir=geodir.name, silent=True)
    return GEOSeries(raw_gse, disease)

def flatten(input_):
    output = {}
    for key, value in input_.items():
        if isinstance(value, list) and len(value) == 1:
            output[key] = value[0]
        else:
            output[key] = value
    return output


def write_metadata_to_clinical_coll(gse, disease, write_db):
    """Write metadata to clinical collection"""
    clinical_collection_name = "{}_geo_meta".format(disease)
    metadata = flatten(gse.metadata)
    clinical_collection = write_db[clinical_collection_name]
    print("Checking to see if main metadata record already exists in mongo...")
    already_exists = clinical_collection.find_one(metadata)
    if already_exists is None:
        print("Record does not exist.")
        clinical_collection.insert_one(metadata)
        print("Inserting metadata into clinical collection (1)...",
              flush=True)
    else:
        print("Record already exists, skipping...")
    for value in gse.gsms.values():
        print("Checking to see if metadata for {} already exists...".format(value.name))
        already_exists = clinical_collection.find_one(flatten(value.metadata))
        if already_exists is None:
            print("Record does not exist, inserting...")
            clinical_collection.insert_one(flatten(value.metadata))
        else:
            print("Record already exists, skipping...")

def write_molecular_collection(gse, gse_df, write_db, force, probe_centric=True):
    """Write to molecular collection"""
    mol_collection_name = gse.get_molecular_collection_name(probe_centric)
    delete_or_exit(write_db, mol_collection_name, force)
    numrows = len(gse.gsms()[list(gse.gsms().keys())[0]].table)
    docs = []
    print("Inserting into {}....".format(mol_collection_name), flush=True)
    for idx in range(len(gse_df.index)): # for each row...
        if idx % 1000 == 0:
            print("Processed {} of {} rows...      \r".format(idx, numrows),
                  end="", flush=True)
        row = gse_df.iloc[idx, :]
        data = {}
        rowname = gse_df.index[idx]
        for col_idx, cell in enumerate(row): # for each column in row...
            data[gse_df.columns[col_idx]] = float(cell)
        doc = dict(id=rowname, data=data, min=float(min(row)), max=float(max(row)))
        docs.append(doc)
    print("                              \r", flush=True)
    print("Done.", flush=True)
    print("here we go")
    write_db[mol_collection_name].insert_many(docs)
    print("really done")


def write_gpl_collection(gse, write_db):
    """Write to GPL collection"""
    gpl = gse.gpl
    gpl_collection_name = "geo_{}".format(gpl)
    if not gpl_collection_name in write_db.collection_names():
        docs = []
        gpl_collection = write_db[gpl_collection_name]
        print("inserting into {}....".format(gpl_collection_name), flush=True)
        gpl_table = get_gpl(gse).table
        numrows = len(gpl_table)
        for idx in range(numrows):
            if idx % 1000 == 0:
                print("Processed {} of {} rows...      \r".format(idx, numrows),
                      end="", flush=True)
            row = gpl_table.iloc[idx, :].to_dict()
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


def get_gpl(gse):
    """Return platform identifier"""
    return gse.gpls[gse.gpl]

def trim_gpl_table(gse, gse_df):
    """
    we are going to add the gene symbol column from the
    GPL table to our gse data frame. The gse data frame
    may have fewer rows than the GPL table, so let's
    first get rid of the corresponding (missing) rows
    in the GPL table.
    """
    gpl_table = get_gpl(gse).table
    rows_to_drop = list(set(gpl_table['ID']) - set(gse_df.index))
    # rather than calling drop, just subset the table
    # to omit those rows, and assign the result back to gpl.table:
    return gpl_table[~gpl_table.ID.isin(rows_to_drop)]


def remove_duplicates(gse_df, gpl_table):
    """
    When working with gene-centric data there may be
    multiple rows with the same gene symbol. Remove
    duplicates by only picking the row with the
    largest variance.
    """
    # add a hugo column to gse_df
    hugo = gpl_table['Gene Symbol']

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
    print("Done.                        ")

    diff = list(set(hugo) - set(dupes))
    which = gse_df.hugo.isin(diff)
    gene_gse = gse_df[which]
    # append is like rbind in R:
    gene_gse = gene_gse.append(gse_df.loc[keep_probe, :])
    gene_gse.index = gene_gse.hugo
    del gene_gse['hugo'] # remove hugo column
    return gene_gse



def ingest(accession, disease, force):
    """
    Ingest data into MongoDb.
    Args:
        accession: GEO accession number
        disease: Oncoscape disease name
        force: Whether to overwrite existing collections, or exit
    """

    geoseries = get_from_geo(accession, disease)

    write_client = pymongo.MongoClient(os.getenv('MONGO_WRITE_URL'))
    write_db = write_client.some_db # TODO FIXME database name

    write_metadata_to_clinical_coll(geoseries, disease, write_db)

    for idx, gse in enumerate(geoseries):
        gpl = gse.gpl
        print("processing series {} of {} ({})".format(idx+1, geoseries.series,
                                                       gpl))
        gse_df = get_data_frame(gse)
        write_molecular_collection(gse, gse_df, write_db, force)
        write_gpl_collection(gse, write_db)

        gpl_table = trim_gpl_table(gse, gse_df)
        gene_gse = remove_duplicates(gse_df, gpl_table)

        write_molecular_collection(gse, gene_gse, write_db, force, False)






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
