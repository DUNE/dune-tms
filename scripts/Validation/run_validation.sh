infile="$1"
shift
base_filename=$(basename "$infile")
outfile="/exp/dune/data/users/${USER}/dune-tms/Validation/Tracking_Validation/${base_filename%.root}.root"
make && ./Tracking_Validation "$infile" $@ && python simply_draw_everything.py "${outfile}"
echo "Output should be in"
echo "${outfile/.root/_images}"
