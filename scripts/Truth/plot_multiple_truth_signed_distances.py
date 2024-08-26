#  plot_multiple_truth_signed_distances.py
# this script is used to plot the sign distance of
# multiple different muon KE ranges.
# Each distribution is from a ROOT file.

# Note: user may need to install some python packages:
# pip install --user matplotlib uproot
# (--force-reinstall if you have a conflict)

# Aug. 2024.
#  M. Dolce, mdolce@fnal.gov.

# For reference, the ROOT file used as input here is created from:
#   dune-tms/scripts/Truth/truth_signed_distance_mu_momentum_slices.py
# That output ROOT file is used here.

# TODO: scale the muon DOWN to match the a-muon, so they normalized to same area.

import argparse
import matplotlib.pyplot as plt
import os
import uproot as up

# parse the arguments
parser = argparse.ArgumentParser(description='Plot the sign distance of muons and anti-muons together.')
parser.add_argument('--infile', type=str, help='The input ROOT file.')
parser.add_argument('--outdir', type=str, help='The output dir name.')
args = parser.parse_args()

def process_args(args_from_argparse):
    if not os.path.exists(args_from_argparse.infile):
        raise FileNotFoundError(f'File {args_from_argparse.infile} not found.')
    print('Infile.......', args_from_argparse.infile)
    out_dir = args_from_argparse.outdir or f"/exp/dune/data/users/{os.environ['USER']}/dune-tms_hists/plot_multiple_truth_signed_distances"
    os.makedirs(out_dir, exist_ok=True)
    print('Outdir.......', out_dir)

    return args_from_argparse.infile, out_dir

infile, outdir = process_args(args)

hists = {}
with up.open(infile) as file:
    print('Ignore anything other than TH1s in the ROOT file...')
    for key in file.keys():

        # clean up the key names
        key = key.split(';')[0]
        hists[key] = file[key]

        # I saved the canvases to the ROOT file, but I don't want them here
        if not hists[key].classname.startswith("TH1"):
            del hists[key]

print('Here are the keys in the file:')
print(hists.keys())

hists_tms_muon, hists_tms_amuon = {}, {}
hists_lar_muon, hists_lar_amuon = {}, {}
print('now go through this dictionary and divide them into groups...')
for key in hists.keys():
    if key.startswith("amuon") and "_tms_ke" in key:
        hists_tms_amuon[key] = hists[key]
    elif key.startswith("muon") and "_tms_ke" in key:
        hists_tms_muon[key] = hists[key]
    elif key.startswith("amuon") and "lar" in key:
        hists_lar_amuon[key] = hists[key]
    elif key.startswith("muon") and "_lar_ke" in key:
        hists_lar_muon[key] = hists[key]
    else:
        print(f'I have a key I don\'t know what to do with: {key}')


# convert to lists of tuples, it is easier to work with
list_hists_tms_muon = [(key, hists_tms_muon[key]) for key in hists_tms_muon.keys()]
list_hists_tms_amuon = [(key, hists_tms_amuon[key]) for key in hists_tms_amuon.keys()]
list_hists_lar_muon = [(key, hists_lar_muon[key]) for key in hists_lar_muon.keys()]
list_hists_lar_amuon = [(key, hists_lar_amuon[key]) for key in hists_lar_amuon.keys()]

# these are the KE ranges of the Muons from truth_signed_distance_mu_momentum.py
# And TH1 names have an index, like: muon_tms_ke_3. The index is the KE range.
ke_ranges = [(0, 250), (250, 500), (500, 750), (750, 1000),
                       (1000, 2000), (2000, 3000), (3000, 4000),
                       (4000, 4250), (4250, 4500), (4500, 4750), (4750, 5000)]

# We have 11 KE ranges, so we need 11 colors.
list_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red']

for hists_muon, hists_amuon in [(list_hists_tms_muon, list_hists_tms_amuon), (list_hists_lar_muon, list_hists_lar_amuon)]:
    assert len(hists_muon) == len(hists_amuon)  # should be the same length
    title_name = 'KE Entering TMS' if hists_muon == list_hists_tms_muon else 'KE at birth in LAr'
    print('plotting', title_name)
    fig = plt.figure(figsize=(10, 6))
    idx = 0
    # loop through muons, but also needs to access the a-muons via their own list
    for (key, value) in hists_muon:
        assert (key == hists_amuon[idx][0].lstrip("a"))  # needs the "a" in front.

        ke = ke_ranges[idx]
        plt.step(value.axis().edges()[:-1], value.values(), where='post', label=rf'{ke[0]} < KE$_\mu^-$ < {ke[1]} MeV', color=list_colors[idx])
        # hists_muon[key].plot(title=key, xtitle='sign distance [cm]', ytitle='counts', color=idx)
        idx += 1

    # now do the a-muons
    idx = 0
    for (key, value) in hists_amuon:
        assert (key.startswith("a") and key == "a" + hists_muon[idx][0])  # needs the "a" in front.
        ke = ke_ranges[idx]
        plt.step(value.axis().edges()[:-1], value.values(), where='post', label=rf'{ke[0]} < KE$_\mu^+$ < {ke[1]} MeV', linestyle='dotted', color=list_colors[idx])
        # hists_amuon[key].plot(title=key, xtitle='sign distance [cm]', ytitle='counts', color=idx)
        idx += 1

    plt.title('TMS Sign Distance')
    plt.xlabel('Sign Distance (mm)')
    plt.ylabel('Events')
    plt.legend(title=title_name, loc='upper left', ncol=2, fontsize=8, bbox_to_anchor=(0.0, 1.0), columnspacing=0.5)

    plt.grid()
    for ext in ['png', 'pdf']:
        plt.savefig(f'{outdir}/plot_multiple_truth_signed_distances_{title_name.strip()[-3:].lower()}_truth_muon_vs_amuon.{ext}', dpi=300, bbox_inches='tight')
    # plt.show()

print('Done.')
