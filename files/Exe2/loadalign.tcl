mol delete all
foreach species {whale-1MBC aplysia-1MBA horse-1WLA seal-1MBS tuna-1MYT turtle-1LHS} {
  mol new $species.psf waitfor all
  mol addfile $species.pdb waitfor all
  mol modstyle 0 top Tube 0.300000 6.000000
  mol modcolor 0 top Molecule
  mol addrep top
  mol modselect 1 top resname HEME
  mol modstyle 1 top Licorice 0.300000 8.000000 6.000000
  mol modcolor 1 top Molecule
}
