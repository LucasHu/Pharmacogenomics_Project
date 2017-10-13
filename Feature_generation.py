from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import numpy as np
#from sklearn.svm import SVC
#from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
#from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
#from sklearn.preprocessing import StandardScaler
#from sklearn.metrics import accuracy_score, cohen_kappa_score, matthews_corrcoef
#from sklearn.externals import joblib

# read sdf file into molecule class in rdkit
def read_sdf_file(directory):
    mols = []
    drugbank_ID = []
    for mol in Chem.SDMolSupplier(directory):
        if mol is not None:
            mols.append(mol)
            drugbank_ID.append(mol.GetProp("DRUGBANK_ID"))
    return mols, drugbank_ID

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


def main():

    # the path of sdf file
    fname = "/Users/lucasminghu/Desktop/Pharmacogenomics/structures.sdf"

    molecules, drugbank_ID = read_sdf_file(fname)

    Mogen2_matrix = generate_fingerprint("Morgan2", molecules)
    
    extra_chemical_descriptor_matrix = get_descriptors(molecules)

    combined_matrix =np.concatenate((Mogen2_matrix, extra_chemical_descriptor_matrix), axis=1)

    return combined_matrix, drugbank_ID





if __name__ == "__main__":
    main()
