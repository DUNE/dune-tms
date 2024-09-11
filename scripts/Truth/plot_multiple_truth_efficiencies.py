#  plot_multiple_truth_signed_distances.py
# this script is used to plot the efficiencies from
# multiple different B-fields.
# Reads in the output ROOT file from
# dune-tms/scripts/Truth/sign_distance_efficiency_vs_muon_true_ke.py.
# So, that script has to be run first for any B field efficiency you want to plot here.

# Note: user may need to install some python packages:
# pip install --user matplotlib uproot
# (--force-reinstall if you have a conflict)

# Sept. 2024.
# M. Dolce, mdolce@fnal.gov.


import argparse
import matplotlib.pyplot as plt
import os
import uproot as up

# parse the arguments
parser = argparse.ArgumentParser(description='Plot the sign distance of muons and anti-muons together.')
parser.add_argument('--indir', type=str, help='The input dir of the ROOT files.')
parser.add_argument('--outdir', type=str, help='The output dir name.')
args = parser.parse_args()

def process_args(args_from_argparse):
    """
    Process the arguments from the command line.
    :param args_from_argparse: the arguments from argparse
    :rtype: dictionary of files, string of outdir
    :returns: infiles, out_dir
    """
    assert os.path.exists(args_from_argparse.indir), f'Directory {args_from_argparse.indir} not found.'
    d_infiles = {}
    for f in os.listdir(args_from_argparse.indir):
        print('file.......', args_from_argparse.indir + '/' + f)
        if not f.endswith('.root'):
            continue
        if not os.path.exists(args_from_argparse.indir + '/' + f):
            raise FileNotFoundError(f'File {args_from_argparse.indir}/{f} not found.')
        # find the B Field from the filename.
        s = f.split('_')
        if 'bfield' in s:
            bfield = s[s.index('bfield') + 1].strip('.root')
            d_infiles[bfield] = args_from_argparse.indir + '/' + f
    out_dir = args_from_argparse.outdir or f"/exp/dune/data/users/{os.environ['USER']}/dune-tms_hists/plot_multiple_truth_efficiencies"
    os.makedirs(out_dir, exist_ok=True)
    print('Outdir.......', out_dir)

    return d_infiles, out_dir

infiles, outdir = process_args(args)


# Four separate dictionaries for {LAr,TMS} x {muon, a-muon}.
d_bfield_th1s_muon_lar, b_field_th1s_muon_tms, b_field_th1s_amuon_lar, b_field_th1s_amuon_tms = {}, {}, {}, {}
for bfield, infile in infiles.items():
    print('Opening file...', infile)
    with up.open(infile) as file:
        for key in file.keys():

            # clean up the key names
            key = key.split(';')[0]

            # I saved the canvases to the ROOT file, but I don't want them here (prob. a canvas)
            if not file[key].classname.startswith("TH1"):
                continue

            if 'lar' in key and key.startswith('muon'):
                d_bfield_th1s_muon_lar[bfield] = file[key]
            elif 'tms' in key and key.startswith('muon'):
                b_field_th1s_muon_tms[bfield] = file[key]
            elif 'lar' in key and key.startswith('amuon'):
                b_field_th1s_amuon_lar[bfield] = file[key]
            elif 'tms' in key and key.startswith('amuon'):
                b_field_th1s_amuon_tms[bfield] = file[key]
            else:
                print('I don\'t knw what to do with this key:', key)


# we have 4 histograms for each B field: muon and amuon in TMS and LAr.
list_colors = ['blue', 'red', 'green', 'orange']

# first two are muons, second two are anti-muons.
for i, d in enumerate([d_bfield_th1s_muon_lar, b_field_th1s_muon_tms, b_field_th1s_amuon_lar, b_field_th1s_amuon_tms]):
    print(i, ' -- ', d)
    fig = plt.figure(figsize=(10, 6))
    color_idx = 0
    for key_bfield, th1 in d.items():
        det = 'tms' if (i % 2 == 0) else 'lar'
        particle = 'muon' if (i <= 1) else 'amuon'
        pdg_string = fr'$\mu^-$' if particle == "muon" else fr'$\mu^+$'
        plt.step(th1.axis().edges()[:-1], th1.values(), where='post', label=f'{pdg_string} @ {key_bfield.replace("p", ".")}', color=list_colors[color_idx])

        plt.title('Signed Distance Efficiency')
        plt.xlabel(f'Muon Kinetic Energy (MeV), in {det.upper().replace("LAR", "LAr")}')
        plt.ylabel('Efficiency')
        plt.legend(title='B Fields', loc='lower left', ncol=3, fontsize=8, columnspacing=1.0)
        plt.xlim(0, 5000)
        plt.grid()
        for ext in ['png', 'pdf']:
            plt.savefig(f'{outdir}/plot_multiple_truth_efficiencies_{particle}_{det}.{ext}', dpi=300, bbox_inches='tight')
        # plt.show()
        color_idx += 1

print('Done.')
