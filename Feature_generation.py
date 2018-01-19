from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import numpy as np
import sys

from propy.GetProteinFromUniprot import GetProteinSequence as gps

from propy.PyPro import GetProDes
import pandas as pd
from cmapPy.pandasGEXpress import setup_GCToo_logger as setup_logger


#from sklearn.svm import SVC
#from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
#from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
#from sklearn.preprocessing import StandardScaler
#from sklearn.metrics import accuracy_score, cohen_kappa_score, matthews_corrcoef
#from sklearn.externals import joblib

# read sdf file into molecule class in rdkit
def read_sdf_file(directory = "structures.sdf"):
    mols = []
    drugbank_ID = []
    INCHI_KEY_ID = []
    INCHI_KEY_ID_drugbank_ID_dict = {}

    for mol in Chem.SDMolSupplier(directory):
        if mol is not None:
            mols.append(mol)
            drugbank_ID.append(mol.GetProp("DRUGBANK_ID"))

            #a dictionary that key: drugbankID; value: InChiID ...
            INCHI_KEY_ID_drugbank_ID_dict[mol.GetProp("INCHI_KEY")] = (mol.GetProp("DRUGBANK_ID"))

    return mols, drugbank_ID, INCHI_KEY_ID, INCHI_KEY_ID_drugbank_ID_dict

# convert fingerprint BitVec to np.array
def rdkit_numpy_convert(fp):
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)

# generate fingerprint and return as a np.array
def generate_fingerprint(fingerprint_type, mols):
    # generate binary Morgan fingerprint with radius 2
    if fingerprint_type == "Morgan2":
        fp = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in mols]
        #fp = [AllChem.GetMorganFingerprintAsBitVect(mols, 2)]
    elif fingerprint_type == "Morgan3":
        fp = [AllChem.GetMorganFingerprintAsBitVect(m, 3) for m in mols]
    elif fingerprint_type == "AtomPair":
        fp = [AllChem.GetHashedAtomPairFingerprintAsBitVect(m) for m in mols]
    elif fingerprint_type == "TopologicalTorsion":
        fp = [AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(m) for m in mols]
    else:
        print "the fingerprint is not supported here!"
        
    x = rdkit_numpy_convert(fp)
    return x

#get extra chemical descriptors
def get_descriptors(mols):
    descr = []
    for m in mols:
        descr.append([Descriptors.MolLogP(m),
                    Descriptors.TPSA(m),
                    Descriptors.NHOHCount(m),
                    Descriptors.NOCount(m),
                    Descriptors.NumHAcceptors(m),
                    Descriptors.NumHDonors(m),
                    Descriptors.NumRotatableBonds(m),
                    Descriptors.NumHeteroatoms(m),
                    Descriptors.FractionCSP3(m)])
    descr = np.asarray(descr)
    return descr

#get all protein sequence descriptors by using
def get_protein_sequence_descriptor(Uniprot_ID_list):

    #write a file has non-compatable protein sequence...
    text_file = open("non_compitable_sequence.txt", "w")

    descr = []
    error_protein_list = []
    for protein in Uniprot_ID_list:

        proseq = gps(protein)
        Des = GetProDes(proseq)
        print protein

        try:
            alldes = Des.GetALL()

            PSD_list = []

            for key in alldes:
                PSD_list.append(alldes[key])
            descr.append(PSD_list)

        except (KeyError, IOError, AttributeError):
            error_protein_list.append(protein)
            text_file.write(protein)
            text_file.write("\n")
            pass

    text_file.sloe()

    descr = np.asarray(descr)
    return descr


#Read DrugBank data and get the training data
def read_DrugBank(fname = "DrugBank_drug_uniprot_links_20180117.csv"):

    protein_drug_list = []
    drug_name_list = []

    with open(fname) as fp:
        next(fp)
        for line in fp:
            line = line.rstrip()
            items = line.split(",")
            UniProt = items[3]
            type = items[2]
            drug = items[0]

            if type == "SmallMoleculeDrug": # only keep the small molecule drugs

                protein_drug = (UniProt, drug)
                protein_drug_list.append(protein_drug)
                drug_name_list.append(drug)

    fp.close()

    return protein_drug_list, drug_name_list

def main():

    # all the protein_drug from DrugBank...
    protein_drug_from_DrugBank_list, drug_name_list = read_DrugBank()

    # the path of sdf file
    fname = "/Users/lucasminghu/Desktop/Pharmacogenomics/structures.sdf"

    molecules, drugbank_ID, INCHI_KEY, INCHI_KEY_ID_drugbank_ID_dict = read_sdf_file(fname)

    # generate the lists of the legible list of drugs and corresponding proteins with legible structure can be found in DrugBank...
    drug_molecule_list = []
    corresponding_protein_list = []

    for item_pair in protein_drug_from_DrugBank_list:
        proteinID = item_pair[0]
        drugID = item_pair[1]
        if  drugID in drugbank_ID:
            index_in_molecules = drugbank_ID.index(drugID)
            drug_molecule_list.append(molecules[index_in_molecules])
            corresponding_protein_list.append(proteinID)

    # Generate protein sequence features...

    protein_sequence_features_matrix = get_protein_sequence_descriptor(corresponding_protein_list)

    print protein_sequence_features_matrix
    print type(protein_sequence_features_matrix)
    print np.shape(protein_sequence_features_matrix)

    Mogen2_matrix = generate_fingerprint("Morgan2", drug_molecule_list)
    #Mogen3_matrix = generate_fingerprint("Morgan3", drug_molecule_list)
    #AtomPair_matrix = generate_fingerprint("AtomPair", drug_molecule_list)

    # get extra chemical features...
    extra_chemical_descriptor_matrix = get_descriptors(drug_molecule_list)



    #np.savetxt('test.txt', Mogen2_matrix, delimiter='\t')

    #extra_chemical_descriptor_matrix = get_descriptors(molecules)

    #combined_matrix =np.concatenate((Mogen2_matrix, extra_chemical_descriptor_matrix), axis=1)
    #
    # print np.shape(extra_chemical_descriptor_matrix)
    # print extra_chemical_descriptor_matrix
    # print len(drugbank_ID)
    #
    # uniprotid = "P48039"
    # proseq = gps(uniprotid)
    #
    # Des = GetProDes(proseq)
    # alldes = Des.GetALL()
    #
    # print len(alldes)
    # print type(alldes)
    #
    # sys.exit()
    # return combined_matrix, drugbank_ID




if __name__ == "__main__":
    main()
