set sel_id [[atomselect top "protein and alpha"] get resid]
set rmsd "residue_rmsd.dat"
set fid [open $rmsd r]
set length [exec wc -l $rmsd]
set length [lindex $length 0]
puts "number of residues: $length"
for {set i 1} {$i<=$length} {incr i} {
gets $fid line
set sel [atomselect top "protein and resid $i"]
$sel set beta [lindex $line 1]
}
close $fid


