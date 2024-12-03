echo "Setting up dune_plot_style"
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dune_plot_style
if [ $? -eq 0 ]; then echo "Setup complete!"
else echo "Setup failed!"
fi
