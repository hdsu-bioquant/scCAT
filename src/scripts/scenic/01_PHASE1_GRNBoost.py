#echo "python /icgc/dkfzlsdf/analysis/B080/crg/keding/scGRN/pyScenicGlioblastoma/pyScenicScriptDarmanis.py" |  qsub -l walltime=72:00:00,mem=250gb,nodes=1:ppn=16 -N pyScenic  -q highmem

#before starting, install pySCENIC from the treminal with:
#pip install pyscenic
#install the required data from
#http://pyscenic.readthedocs.io/en/latest/


#---------Settings------------------------------------------------------------------
print("STARTING")
runOnCluster = True
#runOnCluster = False
#RegulonsViaDask = True
RegulonsViaDask = False
#calculate regulons in Phase II with or without the intermediate dataframe of enriched features
calcRegulonsWithIntermediateDf = True

nCores = 5 

#---------import modules and declare constants--------------------------------------
import os, glob
import argparse
import pandas as pd
import numpy as np

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, save_to_yaml, load_from_yaml
from pyscenic.prune import prune, prune2df
from pyscenic.transform import df2regulons        #was missing - is this right?
from pyscenic.aucell import aucell

import seaborn as sns

from dask.distributed import Client, LocalCluster    #needed to run grnboost2
import pickle                        #for saving and loading "modules" in phase II

print("FINISHED IMPORTING MODULES")


if __name__ == '__main__':
    
    print("run SCENIC")
    
    parser = argparse.ArgumentParser(description='A python script for running SCENIC', epilog='Dependencies: pyscenic, arboreto', usage='prog <filename> [options]')
    parser.add_argument("--database_dir", default=None, help="Path to the HumanRankingDatabases and HumanMotifAnnotation Databases")
    parser.add_argument("--results_dir", default=None, help="Path to the direcotry to save SCENIC results")
    parser.add_argument("--exprs_mat", default=None, help="Path to gene expression matrix")
    parser.add_argument("-id", "-ID", dest="assayID", default=None, type=str, help="String with the assay ID")
    
    args = parser.parse_args()
    
    # Base paths
    DATABASE_FOLDER = args.database_dir
    OUT_FOLDER      = args.results_dir
    TMP_FOLDER      = os.path.join(OUT_FOLDER, "phases")
    ASSAYID         = args.assayID

    # Paths to DataBase
    FEATHER_GLOB            = os.path.join(DATABASE_FOLDER, "HumanRankingDatabases/hg19*.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(DATABASE_FOLDER, "HumanMotifAnnotation/motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
    HG_TFS_FNAME            = os.path.join(DATABASE_FOLDER, "HumanMotifAnnotation/hg_tfs.txt")
    
    # path to the expression matrix
    SC_EXP_FNAME            = args.exprs_mat
    
    # Path to results and intermediate steps
    SC_EXP_FILT_FNAME       = os.path.join(TMP_FOLDER,  "00_exMatFilt.csv")
    ADJACENCIES_FNAME       = os.path.join(TMP_FOLDER,  "01_adjacencies.csv")
    MODULES_BIN_FNAME       = os.path.join(TMP_FOLDER,  "02_modules.bin")
    MODULES_FNAME           = os.path.join(TMP_FOLDER,  "02_modules.txt")
    REGULONS_BIN_FNAME      = os.path.join(TMP_FOLDER,  "03_regulons.bin")
    REGULONS_FNAME          = os.path.join(OUT_FOLDER,  (ASSAYID + "_regulons.txt"))
    AUC_FNAME               = os.path.join(OUT_FOLDER,  (ASSAYID + "_AUC.csv"))
    CLUSTERMAP_FNAME        = os.path.join(OUT_FOLDER,  (ASSAYID + "_Clustermap.png"))
    NOMENCLATURE            = "HGNC"

    print("FINISHED DECLARING CONSTANTS")

    #----------create an expression matrix------------------------------------------------

    #Read in the expression matrix
    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep=',', header=0, index_col=0)
    
    #remove duplicate genes
    ex_matrix = ex_matrix[~ex_matrix.index.duplicated(keep='first')]
    print(ex_matrix.shape)

    #print("CREATED EXPRESSION MATRIX, head:")
    print(ex_matrix.head())
    
    #save filtered expression matrix
    ex_matrix.to_csv(SC_EXP_FILT_FNAME, sep='\t')
    
    #Load filtered expression matrix
    ex_matrix = pd.read_csv(SC_EXP_FILT_FNAME, sep='\t', header=0, index_col=0)

    ex_matrix = ex_matrix.T    
    print("LOADED ex_matrix")
    
    #load TF names
    tf_names = load_tf_names(HG_TFS_FNAME)
    print("FIRST 10 TF NAMES:")
    print(tf_names[0:10])

    #load ranking databases
    db_fnames = glob.glob(FEATHER_GLOB)
    
    #print("DATABASE FILE NAMES:")
    print(db_fnames)


    def name(fname):
        return os.path.basename(fname).split(".")[0]
    #dbs = [RankingDatabase(fname=fname, name=name(fname), nomenclature=NOMENCLATURE) for fname in db_fnames]    
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]    

    print("FEATHER RANKING DATABASE, dbs and type(dbs):")
    print(dbs)
    print(type(dbs))

    print("LOADED TF NAMES AND RANKING DATABASE")

    #------------Phase I: Inference of co-expression modules--------------------------------
    #------------GRNBoost-------------------------------------------------------------------

    print("STARTING PHASE I")
    
    # Define cluster 
    local_cluster = LocalCluster(n_workers=nCores, threads_per_worker=1)     
    client = Client(local_cluster)                
    print(client) 
    
    N_SAMPLES = ex_matrix.shape[0] # Full dataset
    print(N_SAMPLES)
    
    adjacencies = grnboost2(expression_data=ex_matrix.sample(n=N_SAMPLES, replace=False), tf_names=tf_names, seed = 123, verbose=True, client_or_address=client)        
    print("DEFINED adjacencies, type and head:")

    adjacencies.to_csv(ADJACENCIES_FNAME, sep = '\t')
    
    #load adjacencies
    adjacencies  = pd.read_csv(ADJACENCIES_FNAME, sep='\t', header=0, index_col=0)
    print("READ IN adjacencies, type and  head:")
    
    print(type(adjacencies))
    print(adjacencies.head())


    
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    print("DEFINED modules, type:")

    #write modules in a file as binary and as text
    with open(MODULES_BIN_FNAME, "wb") as f:
        pickle.dump(modules, f)

    modules_txt = open(MODULES_FNAME, "w")
    modules_txt.write(str(modules))
    modules_txt.close()

    #read in modules
    with open(MODULES_BIN_FNAME, "rb") as f:
        modules = pickle.load(f)
    print("LOADED modules, type:")

    print(type(modules))
    #print(modules)        #long list
        
    

