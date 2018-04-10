from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import numpy as np
import sys
import random
from sklearn.model_selection import train_test_split
import Deep_learning as DL
import Feature_generation as Fg
import os
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

#get all protein sequence descriptors by using propy
def get_protein_sequence_descriptor(Uniprot_ID_list, no_compitable_sequence_file_name):

    #write a file has non-compatable protein sequence...
    #text_file = open("non_compitable_sequence_uniprot.txt", "w")
    #text_file = open("non_compitable_sequence_drugbank.txt", "w")
    text_file = open(no_compitable_sequence_file_name, "w")


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


#Read PubChem data and get the training data
def read_PubChem(fname = "./PubChem/BindingDB_PubChem.tsv"):

    protein_drug_list = []
    drug_name_list = set()

    with open(fname) as fp:
        next(fp)
        for line in fp:
            line = line.rstrip()
            items = line.split("\t")

            if len(items) > 41 and items[41] != "" :
                UniProt = items[41]
                drug = items[1]
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

# red fragment file...
def read_fragment_file(file = "FDB-17/FDB-17-fragmentset.smi"):
    fragment_list = []
    with open(file) as fp:
        for line in fp:
            line = line.rstrip()
            items = line.split()
            fragment_list.append(items[0])
            print items[0]
    return fragment_list

# Generate test drug input matrix
def Generate_test_drug_input_matrix(Molecules_smile_list, protein_sequence_fingerprint_file):

    drug_molecule_list = []
    for molecule_smile in Molecules_smile_list:
        drug_molecule_list.append(Chem.MolFromSmiles(molecule_smile))

    Mogen2_matrix = generate_fingerprint("Morgan2", drug_molecule_list)

    #protein_sequence_fingerprint = np.loadtxt("./PubChem/Protein_sequence_descriptors_matrix_for_proteins_protein_drug_interaction_from_PubChem.txt")
    protein_sequence_fingerprint = (np.loadtxt(protein_sequence_fingerprint_file))[0:10,:] #test the first ten proteins

    #print "debug here"
    #print np.shape(protein_sequence_fingerprint)
    #sys.exit()

    Mogen3_matrix = generate_fingerprint("Morgan3", drug_molecule_list)
    AtomPair_matrix = generate_fingerprint("AtomPair", drug_molecule_list)

    combined_matrix = np.concatenate((Mogen2_matrix, Mogen3_matrix, AtomPair_matrix), axis=1)

    exp_matrix = np.zeros((np.shape(combined_matrix)[0] * np.shape(protein_sequence_fingerprint)[0], np.shape(combined_matrix)[1] + np.shape(protein_sequence_fingerprint)[1]))

    ########the following has bugs
    i = np.shape(protein_sequence_fingerprint)[0] # i is the number of items in protein_sequence_fingerprint matrix (the no. of proteins)...
    for j in range(np.shape(combined_matrix)[0]):
        for k in range(np.shape(protein_sequence_fingerprint)[0]):
            #exp_matrix[(j + 1) * (k + 1) - 1,:] = np.concatenate((combined_matrix[j, :], protein_sequence_fingerprint[k, :]), axis=0)
            exp_matrix[j * i + k, :] = np.concatenate((combined_matrix[j, :], protein_sequence_fingerprint[k, :]), axis=0)

    return exp_matrix


# Get Frament library
def get_fragment(file = "./Fragment_Library/Enamine_Ro3_Fragments.sdf"):
    file_items = file.split("/")
    actual_file = file_items.pop()
    path = "/".join(file_items)
    cwd = os.getcwd() # get the current directory
    os.chdir(path)
    fragment_library = Chem.SDMolSupplier(actual_file)
    smile_list = []
    for fragment in fragment_library:
        smile = Chem.MolToSmiles(fragment)
        smile_list.append(smile)

    print "the number of fragment molecules is:" + str(len(smile_list))

    os.chdir(cwd) # need to change back

    return smile_list

# Exclude non_compatible list proteins from the whole list...
# the order shouldnt be changed...
def get_final_protein_list(whole_list, non_compatible_list):
    new_list = []
    for item in whole_list:
        if item not in non_compatible_list:
            new_list.append(item)
    return new_list

#read non_compatible protein files into a list...
def get_non_compatible_list(file):
    new_list = []
    with open(file) as fp:
        for line in fp:
            line = line.rstrip()
            items = line.split()
            new_list.append(items[0])
    return new_list

def main():

    # read Enamine_Ro3_Fragments.sdf database...
    fragment_list = get_fragment()

    # read DRUGBANK targets information and exluce the non-reading targets and generate the final targets list
    DRUG_bank_targets = Fg.read_targets_DrugBank()
    non_compatible_DRUG_bank_targets = get_non_compatible_list("non_compitable_sequence_DRUGBANK_targets.txt")
    final_DRUG_bank_targets = get_final_protein_list(DRUG_bank_targets, non_compatible_DRUG_bank_targets)

    new_fragments_matrix = Generate_test_drug_input_matrix(fragment_list, "Protein_sequence_descriptors_matrix_for_proteins_from_DRUGBANK_targets.txt")
    #print new_fragments_matrix
    print "the size of matrix is:"
    print np.shape(new_fragments_matrix)

    # all the protein_drug from DrugBank...
    protein_drug_from_PubChem_list, drug_name_list = read_PubChem() #14945 drug-protein pairs and 4993 drugs in those pairs

    # the path of sdf file
    #fname = "structures.sdf"
    #molecules, drugbank_ID, INCHI_KEY, INCHI_KEY_ID_drugbank_ID_dict = read_sdf_file_drugbank(fname) # in total 8738 drugs in drug bank.
    #molecules_set = set()
    protein_set = set()
    molecules_smiles = set()

    for pairs in protein_drug_from_PubChem_list:
        drug = pairs[1]
        protein = pairs[0]
        #drug_mol = Chem.MolFromSmiles(drug)
        #molecules_set.add(drug_mol)
        protein_set.add(protein)
        molecules_smiles.add(drug)
    print "molecules_set and protein_set are finished!"

    # generate the lists of the legible list of drugs and corresponding proteins with legible structure can be found in PubChem...
    corresponding_protein_list = list(protein_set)
    molecules_smiles_list = list(molecules_smiles)

    print len(corresponding_protein_list) # 533 proteins here
    print len(molecules_smiles_list) # 52802

    # Generate protein sequence features...
    #protein_sequence_features_matrix = get_protein_sequence_descriptor(corresponding_protein_list, "non_compitable_sequence_drugbank.txt")
    #protein_sequence_features_matrix = get_protein_sequence_descriptor(read_uniprot_database(), "non_compitable_sequence_uniprot.txt")
    #protein_sequence_features_matrix = get_protein_sequence_descriptor(corresponding_protein_list, "non_compitable_sequence_PubChem.txt")
    #np.savetxt('Protein_sequence_descriptors_matrix_for_proteins_protein_drug_interaction_from_uniprot.txt', protein_sequence_features_matrix, delimiter='\t')
    #np.savetxt('Protein_sequence_descriptors_matrix_for_proteins_protein_drug_interaction_from_PubChem.txt', protein_sequence_features_matrix, delimiter='\t')
    ###########the part above once done once then can be store into a txt file, no need to run it everytime.


    #read the protein list that can't be processed by propy...
    bad_proteins_drugbank = read_txt_file("./PubChem/non_compitable_sequence_PubChem.txt")

    #exclude the not compitable proteins names from the "corresponding_protein_list", we need the list to map the index to protein sequence features matrix...
    final_protein_list_PubChem = []

    for protein in corresponding_protein_list:
        if protein not in bad_proteins_drugbank:
            final_protein_list_PubChem.append(protein) #final_protein_list_drugbank can be used to map protein names from the final protein sequence features matrix

    # generate drug_molecule_list...
    drug_molecule_list = []
    for molecule_smile in molecules_smiles_list:
        drug_molecule_list.append(Chem.MolFromSmiles(molecule_smile))

    Mogen2_matrix = generate_fingerprint("Morgan2", drug_molecule_list)

    protein_sequence_fingerprint = np.loadtxt("./PubChem/Protein_sequence_descriptors_matrix_for_proteins_protein_drug_interaction_from_PubChem.txt")

    Mogen3_matrix = generate_fingerprint("Morgan3", drug_molecule_list)
    AtomPair_matrix = generate_fingerprint("AtomPair", drug_molecule_list)

    Mogen2_matrix = np.concatenate((Mogen2_matrix, Mogen3_matrix, AtomPair_matrix), axis=1)

    print Mogen2_matrix.shape
    print "molecules fingerprint calculation is done!"


    # get extra chemical features...
    #extra_chemical_descriptor_matrix = get_descriptors(drug_molecule_list)
    #print "chemical descriptors is done."


    # since I dont want to use inefficient "stack" function, I will build a matrix of desriable size to store the information...
    # let's count how many legible drug-protein paris first...
    final_legible_drug_protein_pair_count = 0
    for pairs in protein_drug_from_PubChem_list:
        if pairs[0] in final_protein_list_PubChem:
            final_legible_drug_protein_pair_count += 1

    #column_no = Mogen2_matrix.shape[1] + protein_sequence_fingerprint.shape[1] + extra_chemical_descriptor_matrix.shape[1]
    column_no = Mogen2_matrix.shape[1] + protein_sequence_fingerprint.shape[1]

    # organize a matrix the for training_positive...
    input_matrix_positive = np.zeros((final_legible_drug_protein_pair_count, column_no))

    print "the size of positive training matrix is: "
    print np.shape(input_matrix_positive)

    # read the information from Mogen2_matrix and protein_sequence_fingerprint and store in input_matrix...
    row_no = 0

    # a list to record the drug-protein pair...
    drug_protein_index_pair_exclusive_list = []

    for pairs in protein_drug_from_PubChem_list:

        if pairs[0] in final_protein_list_PubChem:

            protein_index  = final_protein_list_PubChem.index(pairs[0])
            drug_index = molecules_smiles_list.index(pairs[1])
            drug_protein_pair = (drug_index, protein_index)
            drug_protein_index_pair_exclusive_list.append(drug_protein_pair)
            #combined_row = np.concatenate((Mogen2_matrix[drug_index, :], protein_sequence_fingerprint[protein_index , :], extra_chemical_descriptor_matrix[drug_index, :]), axis=0)
            combined_row = np.concatenate((Mogen2_matrix[drug_index, :], protein_sequence_fingerprint[protein_index, :]), axis=0)
            input_matrix_positive[row_no, :] = combined_row
            row_no += 1

    #organize a matrix the for training_negative...
    input_matrix_negative = np.zeros((final_legible_drug_protein_pair_count, column_no))

    count = 0
    while count < np.shape(input_matrix_negative)[0]:
        drug_index = random.randint(0, len(molecules_smiles_list) - 1)
        protein_index = random.randint(0, len(final_protein_list_PubChem) - 1)

        if (drug_index, protein_index) not in drug_protein_index_pair_exclusive_list:
            drug_protein_index_pair_exclusive_list.append((drug_index, protein_index))
            #combined_row = np.concatenate((Mogen2_matrix[drug_index, :], protein_sequence_fingerprint[protein_index, :], extra_chemical_descriptor_matrix[drug_index, :]), axis=0)
            combined_row = np.concatenate((Mogen2_matrix[drug_index, :], protein_sequence_fingerprint[protein_index, :]), axis=0)
            input_matrix_negative[count, :] = combined_row
            count += 1

    input_matrix = np.concatenate((input_matrix_positive, input_matrix_negative), axis = 0)

    list_positive = [1] * np.shape(input_matrix_positive)[0]
    list_negative = [0] * np.shape(input_matrix_negative)[0]

    labels = list_positive + list_negative

    x_train, x_test, y_train, y_test = train_test_split(input_matrix, labels, test_size = 0.25)

    print "we are doing machine leanring now."
    #NADH_list = read_uniprot_database('NADH_dehydrogenase.txt')

    #DL.random_forest(x_train, x_test, y_train, y_test)



    #final_scores = DL.random_forest_prediction(input_matrix, labels, new_fragments_matrix, fragment_list, final_DRUG_bank_targets[0:10], "fragment_prediction_result.txt")
    #print final_scores
    #print type(final_scores)
    #print len(final_protein_list_PubChem)
    #print len(final_scores)

    y_train_DL = np.asarray(y_train).reshape((len(y_train), 1))
    y_test_DL = np.asarray(y_test).reshape((len(y_test), 1))

    DL.deep_learning(x_train, x_test, y_train_DL, y_test_DL)





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
