load prot.pdb
load watmd_out.pdb
show lines
hide everything, resn CVX
hide everything, resn SUR
show spheres, (resn CVX and b>2.0)
alter (resn CVX), vdw=(b/10)
spectrum q, red_white_blue, (resn CVX)
show lines, (resn SUR)
spectrum b, blue_white_green, (resn SUR)
show cartoon
