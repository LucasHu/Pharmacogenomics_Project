from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import numpy as np
import sys
import random
from sklearn.model_selection import train_test_split
import Deep_learning as DL
#import matplotlib.pyplot as plt
#import pylab as pl



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

# read drug_bank sdf file into molecule class in rdkit
def read_sdf_file_drugbank(directory = "structures.sdf"):
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
    text_file = open("non_compitable_sequence_uniprot.txt", "w")
    #text_file = open("non_compitable_sequence_drugbank.txt", "w")


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

        except (KeyError, IOError, AttributeError, ZeroDivisionError):
            error_protein_list.append(protein)
            text_file.write(protein)
            text_file.write("\n")
            pass

    text_file.close()

    descr = np.asarray(descr)
    return descr


#Read DrugBank data and get the training data
def read_DrugBank(fname = "DrugBank_drug_uniprot_links_20180117.csv"):

    protein_drug_list = []
    drug_name_list = set()

    with open(fname) as fp:
        next(fp)
        for line in fp:
            line = line.rstrip()
            items = line.split(",")
            UniProt = items[3]
            type = items[2]
            drug = items[0]

            if "SmallMoleculeDrug" in type: # only keep the small molecule drugs

                protein_drug = (UniProt, drug)
                protein_drug_list.append(protein_drug)
                drug_name_list.add(drug)

    fp.close()

    return protein_drug_list, drug_name_list

#read all uniprot human protein and return a list of the uniport ids...
def read_uniprot_database(fname = "Human_uniprot_list_reviewed_20180117.tab"):

    uniprot_id_list = []

    with open(fname) as fp:
        for line in fp:
            line = line.rstrip()
            items = line.split("\t")
            Uniprot = items[0]
            uniprot_id_list.append(Uniprot)
    fp.close()

    return uniprot_id_list

#read txt file into list...
#non_compitable_sequence_drugbank.txt
def read_txt_file(file):
    protein_list = []
    with open(file) as fp:
        for line in fp:
            line = line.rstrip()
            items = line.split()
            protein_list.append(items[0])

    return protein_list



def main():

    # all the protein_drug from DrugBank...
    protein_drug_from_DrugBank_list, drug_name_list = read_DrugBank() #14945 drug-protein pairs and 4993 drugs in those pairs

    # the path of sdf file
    #fname = "structures.sdf"

    fname = "chembl_23.sdf"

    #molecules, drugbank_ID, INCHI_KEY, INCHI_KEY_ID_drugbank_ID_dict = read_sdf_file_drugbank(fname) # in total 8738 drugs in drug bank.

    sys.exit()

    # generate the lists of the legible list of drugs and corresponding proteins with legible structure can be found in DrugBank...
    drug_molecule_list = []
    corresponding_protein_list = []
    drug_molecule_list_ID = []

    for item_pair in protein_drug_from_DrugBank_list:
        proteinID = item_pair[0]
        drugID = item_pair[1]
        if  drugID in drugbank_ID:
            index_in_molecules = drugbank_ID.index(drugID)

            if molecules[index_in_molecules] not in drug_molecule_list:
                drug_molecule_list.append(molecules[index_in_molecules])
                drug_molecule_list_ID.append(drugID) #drugbank ID used for tracking molecules list
            if proteinID not in corresponding_protein_list:
                corresponding_protein_list.append(proteinID)

    #print len(drug_molecule_list) # 4751 legible drugs
    #print len(corresponding_protein_list) # 3955 proteins here

    # Generate protein sequence features...
    #protein_sequence_features_matrix = get_protein_sequence_descriptor(corresponding_protein_list)
    #protein_sequence_features_matrix = get_protein_sequence_descriptor(read_uniprot_database())
    #np.savetxt('Protein_sequence_descriptors_matrix_for_proteins_protein_drug_interaction_from_uniprot.txt', protein_sequence_features_matrix, delimiter='\t')
    ###########the part above once done once then can be store into a txt file, no need to run it everytime.

    #read the protein list that can't be processed by propy...
    bad_proteins_drugbank = read_txt_file("non_compitable_sequence_drugbank.txt")

    #exclude the not compitable proteins names from the "corresponding_protein_list", we need the list to map the index to protein sequence features matrix...
    final_protein_list_drugbank = []

    for protein in corresponding_protein_list:
        if protein not in bad_proteins_drugbank:
            final_protein_list_drugbank.append(protein) #final_protein_list_drugbank can be used to map protein names from the final protein sequence features matrix

    Mogen2_matrix = generate_fingerprint("Morgan2", drug_molecule_list)

    protein_sequence_fingerprint = np.loadtxt("Protein_sequence_descriptors_matrix_for_proteins_protein_drug_interaction_from_drugbank.txt")

    Mogen3_matrix = generate_fingerprint("Morgan3", drug_molecule_list)
    AtomPair_matrix = generate_fingerprint("AtomPair", drug_molecule_list)

    Mogen2_matrix = np.concatenate((Mogen2_matrix, Mogen3_matrix, AtomPair_matrix), axis=1)

    print Mogen2_matrix.shape


    # get extra chemical features...
    extra_chemical_descriptor_matrix = get_descriptors(drug_molecule_list)


    # since I dont want to use inefficient "stack" function, I will build a matrix of desriable size to store the information...
    # let's count how many legible drug-protein paris first...
    final_legible_drug_protein_pair_count = 0
    for pairs in protein_drug_from_DrugBank_list:
        if pairs[0] in final_protein_list_drugbank and pairs[1] in drug_molecule_list_ID:
            final_legible_drug_protein_pair_count += 1

    column_no = Mogen2_matrix.shape[1] + protein_sequence_fingerprint.shape[1] + extra_chemical_descriptor_matrix.shape[1]

    # organize a matrix the for training_positive...
    input_matrix_positive = np.zeros((final_legible_drug_protein_pair_count, column_no))

    # read the information from Mogen2_matrix and protein_sequence_fingerprint and store in input_matrix...
    row_no = 0

    # a list to record the drug-protein pair...
    drug_protein_index_pair_exclusive_list = []

    if len(drug_molecule_list_ID) != np.shape(Mogen2_matrix)[0]:
        print "drug_molecule_list_ID and Mogen2_matrix do not match"
        sys.exit()
    elif len(drug_molecule_list_ID) != np.shape(extra_chemical_descriptor_matrix)[0]:
        print "drug_molecule_list_ID and extra_chemical_descriptor_matrix do not match"
        sys.exit()

    for pairs in protein_drug_from_DrugBank_list:
        if pairs[0] in final_protein_list_drugbank and pairs[1] in drug_molecule_list_ID:

            protein_index  = final_protein_list_drugbank.index(pairs[0])
            drug_index = drug_molecule_list_ID.index(pairs[1])
            drug_protein_pair = (drug_index, protein_index)
            drug_protein_index_pair_exclusive_list.append(drug_protein_pair)
            combined_row = np.concatenate((Mogen2_matrix[drug_index, :], protein_sequence_fingerprint[protein_index , :], extra_chemical_descriptor_matrix[drug_index, :]), axis=0)
            input_matrix_positive[row_no, :] = combined_row
            row_no += 1

    #organize a matrix the for training_negative...
    input_matrix_negative = np.zeros((final_legible_drug_protein_pair_count, column_no))

    count = 0
    while count < np.shape(input_matrix_negative)[0]:
        drug_index = random.randint(0, len(drug_molecule_list_ID) - 1)
        protein_index = random.randint(0, len(final_protein_list_drugbank) - 1)

        if (drug_index, protein_index) not in drug_protein_index_pair_exclusive_list:
            drug_protein_index_pair_exclusive_list.append((drug_index, protein_index))
            combined_row = np.concatenate((Mogen2_matrix[drug_index, :], protein_sequence_fingerprint[protein_index, :], extra_chemical_descriptor_matrix[drug_index, :]), axis=0)
            input_matrix_negative[count, :] = combined_row
            count += 1

    input_matrix = np.concatenate((input_matrix_positive, input_matrix_negative), axis = 0)
    list_positive = [1] * np.shape(input_matrix_positive)[0]
    list_negative = [0] * np.shape(input_matrix_negative)[0]
    labels = list_positive + list_negative

    x_train, x_test, y_train, y_test = train_test_split(input_matrix, labels, test_size = 0.25)

    DL.random_forest(x_train, x_test, y_train, y_test)
    #DL.deep_learning(x_train, x_test, y_train, y_test)



    sys.exit()


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
