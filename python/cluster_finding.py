import ROOT as r
import argparse
import numpy as np
from utils import * 
from datetime import datetime
import seaborn as sns
import os
from sklearn.cluster import DBSCAN
import pandas as pd
from scipy import stats
import pickle as pkl
import csv
from array import array

def main():
    input_file = args.IN
    output = args.OUT
    window_size = args.WINDOW
    cut_name = args.CUT
    min_seed_number = args.MINNSEED
    parameters = args.ONLYIN
    twoseed = args.DOUBLESEED
    det = args.DET
    integral = args.INTEGRAL

    if cut_name!='99perc' and cut_name!='weights' and cut_name!='dummy_cut':
        cut_file = pkl.load(open(args.CUTfile,'rb'))
        df_cut = pd.DataFrame(cut_file)

    file_root = r.TFile(args.IN,"READ")
    tree_sig = file_root.Get("events")

    nentries = tree_sig.GetEntries()

    print('Entries before '+str(cut_name), nentries)

    n_energy_bins, energy_bins, en_max = getEnergyBins(det, True)
    time_bin = getTimeBins(det)
    if integral:
        time_bin = integral
    _,alpha_bins = getPitchBins()
    _,L_bins = getLbins()
    time, energy = [], []

    test_energies = set()
    for ev in np.arange(0,nentries):
         tree_sig.GetEntry(ev)
         if len(test_energies) < n_energy_bins:
             test_energies.add(tree_sig.energy)
         else:
             break
    energy_bins = test_energies

    time, alpha, pitch, L, energy, event_idx, cut = [], [], [], [], [], [], []

    for i in np.arange(0,nentries):
        tree_sig.GetEntry(i)

        day = float(str(tree_sig.day)[-2:])
        
        if cut_name=='99perc':
            cut_val = tree_sig.rms99_of_99
            cut.append(cut_val)
            if tree_sig.counts < cut_val:
                continue

        elif cut_name=='weights':
            cut_val = tree_sig.weight
            cut.append(cut_val)
            if cut_val < 2:
                continue

        elif cut_name=='dummy_cut':
            cut_val = 0
            cut.append(cut_val)
            
        else:
            cut_val = df_cut[(df_cut.day.values==day) & (df_cut.L.values== tree_sig.L) & (df_cut.alpha.values == tree_sig.alpha) & (df_cut.energy.values==tree_sig.energy)][cut_name].values[0]
            cut.append(cut_val)
            if tree_sig.counts < cut_val:
                continue

        time_hour = float(tree_sig.time)+24.*(day-1)

        time_sec = time_hour*3600.

        # variales for clustering
        time.append(time_sec)
        alpha.append(tree_sig.alpha)
        L.append(tree_sig.L)
        energy.append(tree_sig.energy)
        event_idx.append(i)

    print('Entries after '+str(cut_name), len(time))

    if twoseed is True:

        time = np.array(time)
        alpha = np.array(alpha)
        L = np.array(L)
        energy = np.array(energy)
        event_idx = np.array(event_idx)

        time_seed, alpha_seed, pitch_seed, L_seed, energy_seed, event_idx_seed, cut_seed = [], [], [], [], [], [], []

        if parameters==['L','Pitch','Energy']:
            for i_ene in energy_bins:
                for i_L in L_bins:
                    for i_alpha in alpha_bins:
                        time_temp = time[(energy == i_ene) & (L == i_L) & (alpha == i_alpha)]
                        alpha_temp = alpha[(energy == i_ene) & (L == i_L) & (alpha == i_alpha)]
                        L_temp = L[(energy == i_ene) & (L == i_L) & (alpha == i_alpha)]
                        energy_temp = energy[(energy == i_ene) & (L == i_L) & (alpha == i_alpha)]
                        event_idx_temp = event_idx[(energy == i_ene) & (L == i_L) & (alpha == i_alpha)]
                        i_seed = 0

                        while i_seed < len(time_temp)-1:
                            if (time_temp[i_seed+1] - time_temp[i_seed]) < 2*time_bin:
                                time_seed.append(time_temp[i_seed])
                                alpha_seed.append(alpha_temp[i_seed])
                                L_seed.append(L_temp[i_seed])
                                energy_seed.append(energy_temp[i_seed])
                                event_idx_seed.append(event_idx_temp[i_seed])

                                time_seed.append(time_temp[i_seed+1])
                                alpha_seed.append(alpha_temp[i_seed+1])
                                L_seed.append(L_temp[i_seed+1])
                                energy_seed.append(energy_temp[i_seed+1])
                                event_idx_seed.append(event_idx_temp[i_seed+1])
                                i_seed += 2
                            else:
                                i_seed += 1

        elif parameters==['L','Energy']:
            for i_ene in energy_bins:
                for i_L in L_bins:
                    time_temp = time[(energy == i_ene) & (L == i_L)]
                    time_temp = time[(energy == i_ene) & (L == i_L)]
                    alpha_temp = alpha[(energy == i_ene) & (L == i_L)]
                    L_temp = L[(energy == i_ene) & (L == i_L)]
                    energy_temp = energy[(energy == i_ene) & (L == i_L)]
                    event_idx_temp = event_idx[(energy == i_ene) & (L == i_L)]
                    i_seed = 0

                    while i_seed < len(time_temp)-1:
                        if (time_temp[i_seed+1] - time_temp[i_seed]) < 2*time_bin:
                            time_seed.append(time_temp[i_seed])
                            alpha_seed.append(alpha_temp[i_seed])
                            L_seed.append(L_temp[i_seed])
                            energy_seed.append(energy_temp[i_seed])
                            event_idx_seed.append(event_idx_temp[i_seed])

                            time_seed.append(time_temp[i_seed+1])
                            alpha_seed.append(alpha_temp[i_seed+1])
                            L_seed.append(L_temp[i_seed+1])
                            energy_seed.append(energy_temp[i_seed+1])
                            event_idx_seed.append(event_idx_temp[i_seed+1])
                            i_seed += 2
                        else:
                            i_seed += 1

        elif parameters==['Energy']:
            for i_ene in energy_bins:
                time_temp = time[(energy == i_ene)]
                alpha_temp = alpha[(energy == i_ene)]
                L_temp = L[(energy == i_ene)]
                energy_temp = energy[(energy == i_ene)]
                event_idx_temp = event_idx[(energy == i_ene)]
                i_seed = 0

                while i_seed < len(time_temp)-1:
                    if (time_temp[i_seed+1] - time_temp[i_seed]) < 2*time_bin:
                        time_seed.append(time_temp[i_seed])
                        alpha_seed.append(alpha_temp[i_seed])
                        L_seed.append(L_temp[i_seed])
                        energy_seed.append(energy_temp[i_seed])
                        event_idx_seed.append(event_idx_temp[i_seed])

                        time_seed.append(time_temp[i_seed+1])
                        alpha_seed.append(alpha_temp[i_seed+1])
                        L_seed.append(L_temp[i_seed+1])
                        energy_seed.append(energy_temp[i_seed+1])
                        event_idx_seed.append(event_idx_temp[i_seed+1])
                        i_seed += 2
                    else:
                        i_seed += 1
        else:
            time_temp = time
            i_seed = 0

            while i_seed < len(time_temp)-1:
                if (time_temp[i_seed+1] - time_temp[i_seed]) < 2*time_bin:
                    time_seed.append(time_temp[i_seed])
                    alpha_seed.append(alpha[i_seed])
                    L_seed.append(L[i_seed])
                    energy_seed.append(energy[i_seed])
                    event_idx_seed.append(event_idx[i_seed])

                    time_seed.append(time_temp[i_seed+1])
                    alpha_seed.append(alpha[i_seed+1])
                    L_seed.append(L[i_seed+1])
                    energy_seed.append(energy[i_seed+1])
                    event_idx_seed.append(event_idx[i_seed+1])
                    i_seed += 2
                else:
                    i_seed += 1

        print('Entries after seed finding ', len(time_seed))

        time = np.array(time_seed)
        alpha = np.array(alpha_seed)
        L = np.array(L_seed)
        energy = np.array(energy_seed)
        event_idx = np.array(event_idx_seed)

        time_idx = time.argsort()
        time = time[time_idx]
        alpha = alpha[time_idx]
        L = L[time_idx]
        energy = energy[time_idx]
        event_idx = event_idx[time_idx]
    

    # clustering algorithm lines
    X = np.stack([np.array(time),np.array(alpha)*10000.,np.array(L)*10000.,np.array(energy)*10000],axis=1)
    Xfit = np.stack([np.array(time)],axis=1)
    if parameters==['L','Pitch','Energy']:
        print("Clustering is run in time/L/alpha/energy only.")
        Xfit = np.stack([np.array(time),np.array(alpha)*10000.,np.array(L)*10000.,np.array(energy)*10000],axis=1)
    elif parameters==['L','Energy']:
        print("Clustering is run in time/L/energy only.")
        Xfit = np.stack([np.array(time),np.array(L)*10000.,np.array(energy)*10000],axis=1)
    elif parameters==['Energy']:
        print("Clustering is run in time/energy only.")
        Xfit = np.stack([np.array(time),np.array(energy)*10000],axis=1)
    else:
        print("Clustering is run in time only.")

    clustering = DBSCAN(eps=window_size, metric='euclidean', min_samples=1, n_jobs=-1).fit(Xfit)
    y_temp = clustering.labels_

    y = []
    for i_y in y_temp:
        if i_y != -1:
            if len(y_temp[y_temp==i_y]) < min_seed_number:
                y.append(-1)
            else:
                y.append(i_y)
        else:
            y.append(i_y)

    y = np.array(y)
    n_clusters = len(set(y)) - (1 if -1 in y else 0)
    print('Nclusters ',n_clusters, np.unique(y))
    y = y.reshape([len(X),1])


    Xy = np.concatenate((X,np.array(event_idx).reshape([len(X),1])),axis=1)
    Xy = np.concatenate((Xy,y),axis=1)

    # nnumber of good clusters
    n_clusters = len(set(clustering.labels_)) - (1 if -1 in clustering.labels_ else 0)
    n_noise = list(clustering.labels_).count(-1)

    good_cluster_list = np.unique(Xy[Xy[:,-1]!=-1][:,-1])
    good_cluster_index = np.arange(0,len(good_cluster_list))
    cluster_dict = dict(zip(good_cluster_list, good_cluster_index))

    start_cluster = []
    end_cluster = []
    L_cluster = []
    alpha_cluster = []

    cluster_index = -1*np.ones(nentries,dtype=int)

    for cls_i in good_cluster_list:
        cluster_entries = Xy[Xy[:,-1] == cls_i]
        start_cluster = cluster_entries[0,4]
        end_cluster = cluster_entries[-1,4]
        alpha_cluster = cluster_entries[0,1]/10000
        L_cluster = cluster_entries[0,2]/10000
        energy_cluster = cluster_entries[0,3]/10000
        for cls_ev in np.arange(start_cluster,end_cluster+1):
            tree_sig.GetEntry(int(cls_ev))
            if parameters==['L','Pitch','Energy']:
                if (tree_sig.L == L_cluster) and (tree_sig.alpha == alpha_cluster) and (tree_sig.energy == energy_cluster):
                    cluster_index[int(cls_ev)] = int(cluster_dict[int(cls_i)])
            elif parameters==['L','Energy']:
                if (tree_sig.L == L_cluster) and (tree_sig.energy == energy_cluster):
                    cluster_index[int(cls_ev)] = int(cluster_dict[int(cls_i)])
            elif parameters==['Energy']:
                if (tree_sig.energy == energy_cluster):
                    cluster_index[int(cls_ev)] = int(cluster_dict[int(cls_i)])
            else:
                cluster_index[int(cls_ev)] = int(cluster_dict[int(cls_i)])

    file_root.Close()

    clsnr_b = array( 'i', [ -1 ] )
    thr_b = array( 'd', [ 0. ] )
    newroot = r.TFile(input_file,"update")
    t = newroot.Get("events")
    clsnr_new = t.Branch('cls_idx', clsnr_b, 'cls_idx/I' )
    thr_new = t.Branch('thr_cut', thr_b, 'thr_cut/D' )

    for i in np.arange(0,nentries):
        t.GetEntry(i)

        clsnr_b[0] = cluster_index[i]
        thr_b[0] = cut[i]

        clsnr_new.Fill()
        thr_new.Fill()

    newroot.Write("", r.TObject.kOverwrite)
    newroot.Close()  

    '''
    for i in np.arange(0,nentries):
        t.GetEntry(i)
        if i in Xy[Xy[:,-1]!=-1][:,3]:
            clsnr_b[0] = cluster_dict[y[Xy[:,3]==i][0][0]]
        else:
            clsnr_b[0] = -1
        thr_b[0] = cut[i]
        clsnr_new.Fill()
        thr_new.Fill()
    newroot.Write("", r.TObject.kOverwrite)
    newroot.Close()
    '''

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--IN', type=str, help='input file')
    parser.add_argument('--OUT', type=str, default='./', help='output file')
    parser.add_argument('--CUT', type=str, help='type of cut 99perc/z_score_more2/z_score_more3/iso_forest/dummy')
    parser.add_argument('--CUTfile', type=str, default='99perc', help='file with thresholds per L/L/alpha/energy')
    parser.add_argument('--WINDOW', type=float, default=10, help='window size [sec]')
    parser.add_argument('--MINNSEED', type=int, default=2, help='minimum number of seeds to form a cluster [#]')
    parser.add_argument('--ONLYIN', nargs='+', choices=['L','Pitch','Energy'], help='cluster in specific parameters.')
    parser.add_argument('--DOUBLESEED', action='store_true', help='Require two subsequent seeds.')
    parser.add_argument('--DET', type=str, choices=['hepp_l','noaa'], help='Require two subsequent seeds.')
    parser.add_argument('--INTEGRAL', type=int, help='Fluxes integrated in seconds.')


    args = parser.parse_args()
    main()
