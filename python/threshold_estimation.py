import ROOT as r
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
from scipy import stats
import pickle
from scipy.stats import iqr
from sklearn.ensemble import IsolationForest
import gc

def main():
    # Reading run parameters
    input_files = args.INFILES
    output_dir = args.OUTDIR
    output_file = args.OUTFILE
    maxnfile = args.MAX_NFILES
    ploton = bool(int(args.PLOTON))

    # list of dict in which cut values will be stored
    z_score_file_cont = []

    for i_file,filename in enumerate(input_files):

        root_file = r.TFile(filename,"READ")
        tree_sig = root_file.Get("events")
        if i_file > maxnfile:
            break

        if ".root" not in filename:
            continue

        print('Processing file '+filename)
        # TChain structure in place since in princle we can estimate the cut 2,3,4 days (now 1 day)
        nentries = tree_sig.GetEntries()

        # This step is necessary when we do not now the exact value of each L,alpha,energy bin.
        # It results in a extra loop over the entries, if we already know bin center values skip it.
        L, alpha, energy, day, tot_count = [],[],[],[],[]

        for i,__ in enumerate(np.arange(0,nentries)):

            tree_sig.GetEntry(i)
            L.append(tree_sig.L)
            alpha.append(tree_sig.alpha)
            energy.append(tree_sig.energy)
            day.append(float(str(tree_sig.day)[-2:]))
            tot_count.append(tree_sig.counts)

        L, alpha, energy, tot_count = np.array(L), np.array(alpha), np.array(energy), np.array(tot_count)
        
        # name of the plot collection, if PLOTON == 1
        name = filename.replace('./','').replace('/','_').replace('.root','')

        if ploton: 
            print('Wrting threshold plots here '+output_dir+'/'+name+'.pdf')
            pdf = PdfPages(output_dir+'/'+name+'.pdf')

        for L_i in np.unique(L):
            
            for i_ene,ene_i in enumerate(np.unique(energy)):
                print('L-shell = '+str(L_i)+', Energy = '+str(ene_i)+' [MeV] ... calculating threshold')
                if ploton: 
                    fig = plt.figure(figsize=(12,7))  
                    plt.suptitle('L-shell = '+str(L_i)+', Energy = '+str(ene_i)+' [MeV]',fontsize = 15)
                
                for i_alpha,alpha_i in enumerate(np.unique(alpha)):

                    count = np.array(tot_count[(L == L_i) & (alpha == alpha_i) & (energy == ene_i)])
                    max_count=1
                    if len(count)>0:
                        max_count = max(count)

                    # skim events with z_score<2 or z_score<3
                    z_score_more2 = count[stats.zscore(count)>2.]
                    z_score_more3 = count[stats.zscore(count)>3.]
                    # create a IsoForest classifier
                    # clf=IsolationForest(n_estimators=200, max_samples='auto', contamination='auto', max_features=1, n_jobs=-1, random_state=42, verbose=0)
                    
                    # if less then 10 events cuts are not reliable, set cuts = 0 (uneffective)
                    if len(count) <= 10:
                        z_score_more2_val = 0
                        z_score_more3_val = 0
                        iso_score = 0
                        perc99 = 0
                    else:
                        # set z_score_more2 cut at the lowest count number with z_score > 2
                        if len(z_score_more2)>0:
                            z_score_more2_val = min(z_score_more2)
                        # discard all bin events = z_score_more2_val = maximum counts value
                        else:
                            z_score_more2_val = max(count)

                        # set z_score_more2 cut at the lowest count number with z_score > 3
                        if len(z_score_more3)>0:
                            z_score_more3_val = min(z_score_more3)
                        # discard all bin events = z_score_more3_val = maximum counts value
                        else:
                            z_score_more3_val = max(count)                        

                        # IsoForest training
                        #clf.fit(count.reshape(len(count),1))
                        #y_pred_train = clf.score_samples(count.reshape(len(count),1))
                        # cut value to select the 1% most anomalous events
                        #if len(count[y_pred_train<=np.quantile(y_pred_train,0.01)])>0:
                        #    iso_score = 0 #min(count[y_pred_train<=np.quantile(y_pred_train,0.01)])
                        #else:
                        #    iso_score = max(count)

                        # 99% value calculation
                        perc99 = np.quantile(count,0.99)

                    if ploton: 
                        plt.subplot(3,3,i_alpha+1)
                        plt.title('Alpha = '+str(alpha_i),fontsize = 15)
                        plt.hist(count,bins = 100, range = [0,max_count])
                        plt.xlabel('counts',fontsize = 15)
                        plt.ylabel('entries',fontsize = 15)
                        plt.plot([z_score_more2_val,z_score_more2_val],[0,plt.gca().get_ylim()[1]],'r--',label='zscore=2')
                        plt.plot([z_score_more3_val,z_score_more3_val],[0,plt.gca().get_ylim()[1]],'b--',label='zscore=3')
                        plt.plot([perc99,perc99],[0,plt.gca().get_ylim()[1]],'g--',label='perc99')
                        #plt.plot([iso_score,iso_score],[0,plt.gca().get_ylim()[1]],'y--',label='iso_score')
                        plt.yscale('log')
                        plt.tight_layout()
                        plt.legend()

                    # dict line with day/bin information and cut values
                    line = {'day':day[-1],'L':L_i,'alpha':alpha_i,'energy':ene_i,'99perc':perc99,'z_score_more2':z_score_more2_val,'z_score_more3':z_score_more3_val,'dummy_cut':0.}
                    z_score_file_cont.append(line)

                if ploton:
                    pdf.savefig()
                    fig.clear()
                    plt.close(fig)
                    # save RAM: important for long jobs
                    del fig

                # save RAM: important for long jobs
                del count

        # save RAM: important for long jobs
        if ploton:
            pdf.close()

    # save a pickle file with the output cuts binned in day/L/alpha/energy
    with open(output_dir+"/"+output_file,"wb") as z_score_file:
        pickle.dump(z_score_file_cont,z_score_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--INFILES', help='input files', default=[], nargs='+')
    parser.add_argument('--OUTDIR', help='output dir', default = "./")
    parser.add_argument('--OUTFILE', help='output file', default = "test.pickle")
    parser.add_argument('--MAX_NFILES', help='maximum number of file', default = 100, type = int)
    # slow the processing a lot (x10), use it only in debug mode for 1/2 days
    parser.add_argument('--PLOTON', help='plotting flag: if 1 save counts histogram + cuts', default = 0, type = int)
    args = parser.parse_args()
    main()
