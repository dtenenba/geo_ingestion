# Oncoscape data ingestion

This is ultimately intended to be a framework for ingesting data into
Oncoscape from various sources. At present,
[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) is the only
supported data source. When other data sources are added, the code
will be generalized and refactored to accommodate them (and the repository
will probably be renamed).

## Installation

Python 3 is required.

### Clone the repo

```sh
https://github.com/dtenenba/geo_ingestion.git
cd geo_ingestion
```

### Create a virtual environment

We recommend using Python virtual environments, specifically with
[virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/).

Install and configure virtualenvwrapper, then configure an environment
for use with Python 3 as follows:

```sh
mkvirtualenv --python $(which python3) geo_ingestion
```

In future sessions working with this virtual environment, you don't
need to run that command again, but you do need to switch to the
virtual environment with:

```sh
workon geo_ingestion # tab-completion will help select the virtualenv
```

### Install dependencies

```sh
pip install -r requirements.txt
```

### Set up environment variables

```sh
cp setup_env.sh.example setup_env.sh
```

Then edit `setup_env.sh` to match your system.
At present, only `MONGO_WRITE_URL` is used,
and it defaults to `mongodb://localhost`, so if you are intending to write to the mongo DB
instance running on your local system, and you
have a `mongod` daemon running, you don't need
to edit the file further.

Then source the file with the following command:

```bash
. setup_env.sh
```

**NOTE**: The code writes to the `tcga` database.

### Run the code

Running `python ingester.py -h` gives you the following help message:

```
usage: ingester.py [-h] [-f] ACCESSION DISEASE

ingest GEO data

positional arguments:
  ACCESSION   GEO accession number
  DISEASE     disease name

optional arguments:
  -h, --help  show this help message and exit
  -f          overwrite existing collections
```

So, an example run would look like this:

```
python ingester.py GSE12102 sarc
```

Add the `-f` option if you want to overwrite
existing collections.

Given the accession number and disease from
the above example, this command will create the following collections:

* `sarc_geo_meta`  - Clinical collection (metadata)
* `sarc_geo_GSE12102_gene` - molecular collection indexed by gene symbol
* `sarc_geo_GSE12102` - molecular collection indexed by probe ID
* `geo_GPL570` the platform used by this data set

If a GEO series uses more than one platform (as does GSE16102),
then
there will be two molecular collections for each platform used:

* `sarc_geo_GSE16102-GPL96`
* `sarc_geo_GSE16102-GPL96_gene`
* `sarc_geo_GSE16102-GPL3979`
* `sarc_geo_GSE16102-GPL3979_gene`
