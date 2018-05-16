import pandas as pd, numpy as np, scipy, random, statistics, sys, os
from statistics import mean
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.tree import export_graphviz

# Last edit: May 3 2018

# This script takes in a data frame where 1st column is gene names, last 7 columns are the asymp samples
# Training set -> columns between the first and last 7 columns
# See process_df.py -> to see what the df will look like

# You can choose to take top genes or random genes, shuffle labels, match group sizes
# You can input # genes, # of trees generated, training group size, # iterations
# Script will return accuracy and std from # iterations, # trees

def readfile(dataframe):
    df = pd.read_csv(dataframe, sep=",", comment='#')
    return df

def create_forest(data, labels, size, n):
    train_features, test_features, train_labels, test_labels = train_test_split(data, labels, train_size = size, stratify=labels)
    train_con = len([ _ for _ in train_labels if _.startswith('C')])
    train_HD = len([ _ for _ in train_labels if _.startswith('H')])
    test_con = len([ _ for _ in test_labels if _.startswith('C')])
    test_HD = len([ _ for _ in test_labels if _.startswith('H')])
    clf = RandomForestClassifier(n_estimators = n)
    clf.fit(train_features, train_labels)
    pred = clf.predict(test_features)
    train_acc = accuracy_score(train_labels, clf.predict(train_features))
    test_acc = accuracy_score(test_labels, pred)
    conf_matrix = confusion_matrix(test_labels, clf.predict(test_features))
    spec = conf_matrix[0][0] / (conf_matrix[0][0] + conf_matrix[0][1])
    sens = conf_matrix[1][1] / (conf_matrix[1][0] + conf_matrix[1][1])
    false_pos = 1 - sens
    false_neg = 1 - spec
    return train_acc, test_acc, sens, spec, false_pos, false_neg, clf, train_con, train_HD, test_con, test_HD

def test_asymp(data, labels, clf):
    labels, clf_pred = list(labels), clf.predict(data)
    hdba9, ccap, hdcap = labels.count('HDpos_BA9'), labels.count('C_CAP'), labels.count('HDpos_CAP')
    aBA9_CBA9, aBA9_HDBA9, aCAP_CBA9, aCAP_HDBA9, CCAP_CBA9, CCAP_HDBA9 = 0, 0, 0, 0, 0, 0
    for i in range(len(data)):
        if labels[i] == 'HDpos_BA9' and  clf_pred[i] == 'C_BA9':
            aBA9_CBA9 += 1
        elif labels[i] == 'HDpos_BA9' and  clf_pred[i] == 'HD_BA9':
            aBA9_HDBA9 += 1
        elif labels[i] == 'HDpos_CAP' and  clf_pred[i] == 'C_BA9':
            aCAP_CBA9 += 1
        elif labels[i] == 'HDpos_CAP' and  clf_pred[i] == 'HD_BA9':
            aCAP_HDBA9 += 1
        elif labels[i] == 'C_CAP' and  clf_pred[i] == 'C_BA9':
            CCAP_CBA9 += 1
        elif labels[i] == 'C_CAP' and  clf_pred[i] == 'HD_BA9':
            CCAP_HDBA9 += 1
    aBA9_CBA9, aCAP_HDBA9, CCAP_CBA9 = aBA9_CBA9/hdba9, aCAP_HDBA9/ccap, CCAP_CBA9/hdcap
    aBA9_HDBA9, aCAP_CBA9, CCAP_HDBA9 = aBA9_HDBA9/hdba9, aCAP_CBA9/ccap, CCAP_HDBA9/hdcap
    return aBA9_CBA9, aCAP_HDBA9, CCAP_CBA9, aBA9_HDBA9, aCAP_CBA9, CCAP_HDBA9

def create_set(df, n_genes, e):
    df = df.drop_duplicates(keep='last')
    # Saves genes to list
    top_n = df['Gene name'].head(n=n_genes).tolist()
    all_genes = df['Gene name'].tolist()
    rand_genes = random.sample(all_genes, n_genes)
    assert(len(rand_genes) == len(top_n) == n_genes)
    
    # Set gene name as index
    df = df.set_index(df['Gene name'])
    df.index.name = None
    df = df.drop(['Gene name'], 1)
    # Save sample IDs
    hdpos_cap, train_set = list(df)[-7:], list(df)[:-7]
    # Match sample group sizes
    if e == 'equal':
        control = [ _ for _ in train_set if _.startswith('C')]
        HD = [ _ for _ in train_set if _.startswith('H')]
        if len(HD) < len(control):
            rand_control = random.sample(control, len(HD))
            train_set = rand_control + HD
        if len(control) < len(HD):
            rand_HD = random.sample(HD, len(control))
            train_set = control + rand_HD
        if len(control) == len(HD):
            train_set = control + HD

    # Make labels
    rf_train, rf_predict = df[train_set].T.copy(), df[hdpos_cap].T.copy()
    rf_train.loc[:,'Label'] = rf_train.index.map(lambda x:"HD_BA9" if x.startswith('H') else 'C_BA9')
    rf_predict.loc[:,'Label'] = rf_predict.index.map(lambda x:"HDpos_BA9" if 'BA9' in x and x.startswith('H') else ('HDpos_CAP') if x.startswith('H') else 'C_CAP')
    rf_train_labels, rf_predict_labels = rf_train.pop('Label'), rf_predict.pop('Label')
    # Saves df with the n genes
    rf_train_n, rf_predict_n = rf_train[top_n], rf_predict[top_n]
    rf_train_rand, rf_predict_rand = rf_train[rand_genes], rf_predict[rand_genes]
    return rf_train_n, rf_predict_n, rf_train_rand, rf_predict_rand, rf_train_labels, rf_predict_labels
        
def repeat(n, all_df, p, est, n_genes, top_rand, shuff, eq):
    train_acc, test_acc, sens, spec, false_pos, false_neg = [], [], [], [], [], []
    aBA9_CBA9, aCAP_HDBA9, CCAP_CBA9, aBA9_HDBA9, aCAP_CBA9, CCAP_HDBA9 = [], [], [], [], [], []
    for i in range(n):
        cs = create_set(all_df, n_genes, eq)
        if shuff == 'shuffle':
            random.shuffle(cs[4])
        if top_rand == 'top':
            cf = create_forest(cs[0], cs[4], p, est)
            ta = test_asymp(cs[1], cs[5], cf[6])
        if top_rand == 'random':
            cf = create_forest(cs[2], cs[4], p, est)
            ta = test_asymp(cs[3], cs[5], cf[6])
        train_acc += [cf[0]]
        test_acc += [cf[1]]
        sens += [cf[2]]
        spec += [cf[3]]
        false_pos += [cf[4]]
        false_neg += [cf[5]]
        aBA9_CBA9 += [ta[0]]
        aCAP_HDBA9 += [ta[1]]
        CCAP_CBA9 += [ta[2]]
        aBA9_HDBA9 += [ta[3]]
        aCAP_CBA9 += [ta[4]]
        CCAP_HDBA9 += [ta[5]]
#    return(mean(train_acc), mean(test_acc), mean(sens), mean(spec), mean(false_pos), \
#            mean(false_neg), mean(aBA9_CBA9), mean(aCAP_HDBA9), mean(CCAP_CBA9), \
#            mean(aBA9_HDBA9), mean(aCAP_CBA9), mean(CCAP_HDBA9), statistics.stdev(train_acc), \
#            statistics.stdev(test_acc), statistics.stdev(sens), statistics.stdev(spec), \
#            statistics.stdev(false_pos), statistics.stdev(false_neg), statistics.stdev(aBA9_CBA9), \
#            statistics.stdev(aCAP_HDBA9), statistics.stdev(CCAP_CBA9), statistics.stdev(aBA9_HDBA9), \
#            statistics.stdev(aCAP_CBA9), statistics.stdev(CCAP_HDBA9), cf[7], cf[8], cf[9], cf[10])
    return round(mean(train_acc),3), round(mean(test_acc),3), round(mean(sens),3), round(mean(spec),3), \
            round(mean(false_pos),3), round(mean(false_neg),3), round(mean(aBA9_CBA9),3), \
            round(mean(aCAP_HDBA9),3), round(mean(CCAP_CBA9),3), round(mean(aBA9_HDBA9),3), \
            round(mean(aCAP_CBA9),3), round(mean(CCAP_HDBA9),3), round(statistics.stdev(train_acc),3), \
            round(statistics.stdev(test_acc),3), round(statistics.stdev(sens),3), \
            round(statistics.stdev(spec),3), round(statistics.stdev(false_pos),3), \
            round(statistics.stdev(false_neg),3), round(statistics.stdev(aBA9_CBA9),3), \
            round(statistics.stdev(aCAP_HDBA9),3), round(statistics.stdev(CCAP_CBA9),3), \
            round(statistics.stdev(aBA9_HDBA9),3), round(statistics.stdev(aCAP_CBA9),3), \
            round(statistics.stdev(CCAP_HDBA9),3), cf[7], cf[8], cf[9], cf[10]

def write_f(res, n_genes, set_size, n_trees, n_repeats, top_rand, shuff, eq):
    add_name = '_shuffle'
    add_eq = '_equal'
    if shuff == 'ns':
        shuff = 'original'
        add_name = ''
    if eq != 'equal':
        eq = 'unequal'
        add_eq = ''
    
    name = '../../samples/Analysis_Results/Random_forest_results/rf_' + str(top_rand) + '_' + str(n_genes) + '_genes_' + str(n_trees)+ '_trees_' + \
            str(n_repeats) + '_repeats' + add_name + add_eq + '.txt'
    file = open(os.path.abspath(name), 'w')
    file.write('Random forest: ' + eq + ' group sizes, '+ str(set_size*100) + '% training set, ' + \
            top_rand + ' ' + str(n_genes) + ' genes, ' + str(n_trees) + ' trees, ' + \
            str(n_repeats) + ' iterations \n')
    file.write('Train control: ' + str(res[24]) + ', Train HD: ' + str(res[25]) + \
            ', Test control: ' + str(res[26]) + ', Test HD: ' + str(res[27]) + '\n')
    file.write('\n')
    file.write('Results for ' + top_rand + ' '  + str(n_genes) + ' genes, ' + shuff + ' labels: \n')
    file.write('Mean train accuracy: ' +  str(res[0]) + ' ,std: '+ str(res[12]) + '\n')
    file.write('Mean test accuracy: ' +  str(res[1]) + ' ,std: '+ str(res[13]) + '\n')
    file.write('Mean sensitivity: ' +  str(res[2]) + ' ,std: '+ str(res[14]) + '\n')
    file.write('Mean specificity: ' +  str(res[3]) + ' ,std: '+ str(res[15]) + '\n')
    file.write('Mean false positive: ' +  str(res[4]) + ' ,std: '+ str(res[16]) + '\n')
    file.write('Mean false negative: ' +  str(res[5]) + ' ,std: '+ str(res[17]) + '\n')
    file.write('\n')
    file.write('Mean hdpos BA9 -> control BA9: ' +  str(res[6]) + ' ,std: '+ str(res[18]) + '\n')
    file.write('Mean hdpos CAP -> HD BA9: ' + str(res[7]) + ' ,std: '+ str(res[19]) + '\n')
    file.write('Mean control CAP -> control BA9: ' + str(res[8]) + ' ,std: '+ str(res[20]) + '\n')
    file.write('Mean hdpos BA9 -> HD BA9: ' + str(res[9]) + ' ,std: '+ str(res[21]) + '\n')
    file.write('Mean hdpos CAP -> control BA9: ' + str(res[10]) + ' ,std: '+ str(res[22]) + '\n')
    file.write('Mean control CAP -> HD BA9: ' + str(res[11]) + ' ,std: '+ str(res[23]) + '\n')
    file.close()

def main():
    if len(sys.argv) < 9:
        print("you must call program as:  ")
        print("rf_sympBA9.py <CSV file> <n genes> <train size> <n trees> <n iterations> <top or random> <shuffle or ns> <equal or ne>")
        print("argument 6: enter top or random [top n genes or random n genes]")
        print("argument 7: enter shuffle or ns [shuffle labels or ns: original labels]")
        print("argument 8: enter equal or ne [match sample sizes between groups or ne: unequal size between groups]")
        print("e.g. python rf_sympBA9.py counts.csv 250 0.75 20000 1000 random shuffle ne")
        print("")
        sys.exit(1)

    fl = sys.argv[1]
    n_genes = int(sys.argv[2])
    set_size = float(sys.argv[3])
    n_trees = int(sys.argv[4])
    n_repeats = int(sys.argv[5])
    top_rand = str(sys.argv[6]).lower()
    shuff = str(sys.argv[7]).lower()
    eq = str(sys.argv[8]).lower()
    assert(n_genes != 0 or n_trees != 0 or n_repeats != 0)
    
    df = readfile(fl)
    res = repeat(n_repeats, df, set_size, n_trees, n_genes, top_rand, shuff, eq)
    wf = write_f(res, n_genes, set_size, n_trees, n_repeats, top_rand, shuff, eq)

if __name__ == "__main__":
    main()
