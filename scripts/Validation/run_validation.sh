infile="$1"
n="$2"
base_filename=$(basename "$infile")
outfile="/exp/dune/data/users/${USER}/dune-tms/Validation/Tracking_Validation/${base_filename}"
make && ./Tracking_Validation "$infile" $n && python simply_draw_everything.py "${outfile}"
echo "Output should be in"
echo "${outfile/.root/_images}"
