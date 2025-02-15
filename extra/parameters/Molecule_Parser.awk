function reset_all_arrays(dummy){
    split("",header_set,":")         # Set of lines from the header of a mol2 file
    split("",atom_label,":")         # Array containing the atom labels, index by atom number in the molecule
    split("",coordinates,":")        # Array containing the atomic coordinates, index by atom number then x=1, y=2 and z=3
    split("",atom_type,":")          # Atom type, initial as given in the source file - can be modified, index by atom number in the molecule
    split("",res_num,":")            # Residue number as given in the source file, index by atom number in the molecule
    split("",res_name,":")           # Residue name as given in the source file, index by atom number in the molecule
    split("",atomic_charge,":")      # Atomic charge, initial as given in the source file - can be modified, index by atom number in the molecule
    split("",atomic_tag,":")         # Atom tag, initial as given in the source file - can be modified, index by atom number in the molecule
    split("",connected,":")          # Connection array, indexed byatom number, feild 0 is the number of connection, feilds 1-n are the atom numbers connected
    split("",bondorder,":")          # Bond order array, indexed byatom number, feild 0 is the number of connection, feilds 1-n are the bond orders to the atoms in the connected array -2=un -1=du 0=nc 1=single 1.5=ar 2=double 3=triple
    split("",substructure,":")       # Set of lines from the substructure of a mol2 file
    split("",a1,":")                 # temp array used for splitting : seperated strings into individual values
    split("",stripped_connected,":") # Connection array but with all chain atoms removed, only rings and linker between rings remain
    split("",linker_connected,":")   # Connection array but with all chain and ring atoms removed, only linkers between rings remain
    split("",chain_atom,":")         # Array of chain atom tags, index by atom number in the molecule (C.0 terminal atom, C.1 next to terminal...)
    split("",random_path_array,":")  # Array containing atom numbers defining a random walk through the stripped molecule (no chain atoms)
    split("",Found_Rings,":")        # Array containing a list of all rings found n1:n2:n3:n4... where n1 is an atom number in the ring - contains duplicates
    split("",Clean_Rings,":")        # Array containing a list of rings found n1:n2:n3:n4... where n1 is an atom number in the ring - contains fewwer duplicates
    split("",Unique_Rings,":")       # Array containing a list of unique rings found n1:n2:n3:n4... where n1 is an atom number in the ring - contains NO duplicates
    split("",Clean_Ring_Systems,":") # Array containing list of unique rings that make up fused ring systems
    split("",PiString,":")           # Array containing the final string of Pi electrons for each Unique ring
    split("",RingOverlaps,":")       # Array containg the overlaps (0 or 1) between unique rings r1,r2 (a fused ring system is made of at least three unique rings)
    split("",Linker_atoms,":")       # Array containing a list of linker atoms, 1-m where m is the number of linker atoms, each element is an atom number
    split("",Ring_atoms,":")         # Array containing a list of ring atoms, 1-m where m is the number of ring atoms, each element is an atom number
    split("",Ring_to_Ring_link,":")  # Array containing a list of ring atoms connected to a DIFFERENT ring, 1-m where m is the number of ring-to-ring atoms, each element is an atom number
    split("",Internal_atoms,":")     # Array containing atoms within a ring system which have bonds internal to the ring (ring_system_number, internal_atom_number) element 0 is the number of internal atoms. 
    split("",Ring_System_Properties,":") # num_core_rings, ARcount, number_of_internal_atoms, Total_AR_atom_list, Internal_atoms_list" :  "Largest ring atom list,Bond_Order_String (for largest ring)
    split("",Total_Valance,":")      # Array containing the total number of connected atoms to a given atom, index by atom number in the molecule
    split("",Atomic_Symbol,":")      # Array containing the Atomic Symbol of the atoms, index by atom number in the molecule
    split("",Ring_Type,":")          # Array containing the Ring Type for the Unique_Rings indexed by 1-n where n is the number of unique rings
    split("",Ring_Type_Atom,":")     # Array containing the Ring Type String for each atom in rings, index is the atom number
    split("",desc,":")               # Array containing the Final atom descriptors to be used in patterns
    split("",main,":")               # Array containing main[i,j]:i is atom number j=0 number of elements(4) j=1: atom number (i) j=2: descriptors for atom i j=3 descriptors for connected atoms (separated by /) j=4 connection string (atom numbers separated by :)
}
#
function read_single_mol2(file_name){
    flag=0
    while ((getline < file_name) > 0){
	if ((substr($1,1,1)=="@")&&(substr($1,3,6)=="TRIPOS")&&(substr($1,10,8)=="MOLECULE")){
	    flag=1
	    pos1=1
	    getline < file_name}
	if ((substr($1,1,1)=="@")&&(substr($1,3,6)=="TRIPOS")&&(substr($1,10,4)=="ATOM")){
	    flag=2
	    pos2=1
	    getline < file_name}
	if ((substr($1,1,1)=="@")&&(substr($1,3,6)=="TRIPOS")&&(substr($1,10,4)=="BOND")){
	    flag=3
	    minatom=0
	    maxatom=0
	    getline < file_name}
	if ((substr($1,1,1)=="@")&&(substr($1,3,6)=="TRIPOS")&&(substr($1,10,12)=="SUBSTRUCTURE")){
	    flag=4
	    pos4=1
	    getline < file_name}
	if (flag==1){header_set[pos1]=$0;pos1++}
	if (flag==2){
	    atom_label[$1]=$2
	    coordinates[$1,1]=$3
	    coordinates[$1,2]=$4
	    coordinates[$1,3]=$5
	    atom_type[$1]=$6
	    res_num[$1]=$7
	    res_name[$1]=$8
	    atomic_charge[$1]=$9
	    atomic_tag[$1]=$10
	    pos2++}
	if (flag==3){
	    bond_order=$4
	    if (bond_order==1){bond_order=1}
	    else if (bond_order==2){bond_order=2}
	    else if (bond_order==3){bond_order=3}
	    else if (bond_order=="am"){bond_order=1}
	    else if (bond_order=="ar"){bond_order=1.5}
	    else if (bond_order=="nc"){bond_order=0}
	    else if (bond_order=="du"){bond_order=-1}
	    else if (bond_order=="un"){bond_order=-2}
	    if (connected[$2,0]==""){connected[$2,0]=1;connected[$2,1]=$3;bondorder[$2,1]=bond_order}
	    else {point1=(connected[$2,0])+1;connected[$2,point1]=$3;connected[$2,0]=point1;bondorder[$2,point1]=bond_order}
	    if (connected[$3,0]==""){connected[$3,0]=1;connected[$3,1]=$2;bondorder[$3,1]=bond_order}
	    else {point2=(connected[$3,0])+1;connected[$3,point2]=$2;connected[$3,0]=point2;bondorder[$3,point2]=bond_order}
	    if (minatom==0){minatom=$2}
	    if (maxatom==0){maxatom=$2}
	    if ($2 < minatom){minatom=$2}
	    if ($3 < minatom){minatom=$3}
	    if ($2 > maxatom){maxatom=$2}
	    if ($3 > maxatom){maxatom=$3}}
	if (flag==4){
	    substructure[pos4,1]=$1
	    substructure[pos4,2]=$2
	    substructure[pos4,3]=$3
	    substructure[pos4,4]=$4
	    substructure[pos4,5]=$5
	    substructure[pos4,6]=$6
	    substructure[pos4,7]=$7
	    substructure[pos4,8]=$8
	    pos4++}}
    minmax_return=""minatom":"maxatom":"pos1":"pos4
    return minmax_return}
#
function terminal_atoms(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    i=minatom
    while (i<=maxatom){
	j=1
	stripped_connected[i,0]=connected[i,0]
	while (j<=connected[i,0]){
	    stripped_connected[i,j]=connected[i,j]
	    j++}
	i++}
    loop=1
    terminal_atom_count=0
    still_terminal=1
    while (still_terminal==1){
	i=minatom
	still_terminal=0
	while (i<=maxatom){
	    if (stripped_connected[i,0]==1){
		chain_atom[i]="C."loop
		connected_to=stripped_connected[i,1]
		chain_atom[connected_to]="CC."loop
		stripped_connected[i,0]=0
		terminal_atom_count++}
	    i++}
	loop++
	i=minatom
	while (i<=maxatom){
	    j=1
	    sum=0
	    list_atoms=""
	    while (j<=stripped_connected[i,0]){
		check_atom=stripped_connected[i,j]
		if (substr(chain_atom[check_atom],1,2)!="C."){
		    sum++
		    list_atoms=list_atoms""check_atom":"}
		j++}
	    stripped_connected[i,0]=sum
	    if (sum==1){still_terminal=1}
	    split(list_atoms,a3,":")
	    j=1
	    while (j<=sum){
		stripped_connected[i,j]=a3[j]
		j++}
	    i++}}
    return terminal_atom_count}
#
#
function make_random_path(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    i=minatom
    num_atoms_not_chain=0
    while (i<=maxatom){
	if (stripped_connected[i,0]>0){
	    Start_atom=i
	    num_atoms_not_chain++}
	i++}
    maxloop=500*num_atoms_not_chain
    i=1
    old_step=""
    random_path_list=""
    while (i<=maxloop){
	accepted=0
	while (accepted==0){
	    conn=int((stripped_connected[Start_atom,0]*rand())+1)
	    next_step=stripped_connected[Start_atom,conn]
	    if (next_step != old_step){
		random_path_array[i]=Start_atom
		old_step=Start_atom
		Start_atom=next_step
		accepted=1}}
	i++}
    return maxloop}
#
function find_all_rings_in_one_pass(minmax, maxloop){
    i=1
    ring_count=0
    split("",found_atoms,":")  # initialize found_atoms array
    ring_string=""
    while (i<=maxloop){
	current_atom=random_path_array[i]
	if ((current_atom in found_atoms)==1){
	    loop_end=i
	    loop_start=found_atoms[current_atom]
	    if ((loop_end-loop_start)>2){
	       j=loop_start
	       while (j<loop_end){
		   ring_string=ring_string""random_path_array[j]":"
		   j++}
	       ring_count++
	       Found_Rings[ring_count]=ring_string
	       ring_string=""
	       i=loop_end
	       split("",found_atoms,":")}}  # initialize found_atoms array
	else if ((current_atom in found_atoms)==0){
	    found_atoms[current_atom]=i}
	i++}
    return ring_count}
# 
#
function remove_straight_dupe_rings(dummy){
    n=asort(Found_Rings)
    i=1
    old=""
    clean_ring_count=0
    while (i<=n){
	if (Found_Rings[i]!=old){
	    clean_ring_count++
	    Clean_Rings[clean_ring_count]=Found_Rings[i]
	    old=Found_Rings[i]}
	i++}
    return clean_ring_count}
#
function compare_rings(ring1, ring2){
    size_ring_1=split(ring1,RG1,":")
    size_ring_2=split(ring2,RG2,":")
    if (size_ring_1 != size_ring_2){
	matched_rings=0}
    else {
	ii=1
	match_count=0
	while (ii<=size_ring_1){
	    check_atom=RG1[i]
	    jj=1
	    while (jj<=size_ring_2){
		if (RG1[ii]==RG2[jj]){
		    match_count++}
		jj++}
	    ii++}
	if (match_count==size_ring_1){
	    matched_rings=1}
	else {matched_rings=0}}
    return matched_rings}
#
function remove_complex_dupe_rings(clean_ring_count){
    i5=1
    unique_rings=0
    while (i5<=clean_ring_count){
	if (Clean_Rings[i5]!=""){
	    unique_rings++
	    Unique_Rings[unique_rings]=Clean_Rings[i5]}
	j5=i5+1
	while (j5<=clean_ring_count){
	    check_rings=compare_rings(Clean_Rings[i5], Clean_Rings[j5])
	    if (check_rings==1){
		Clean_Rings[j5]=""}
	    j5++}
	i5++}
    return unique_rings}
#
function is_in_ring(atom_to_check, unique_rings){
    i6=1
    is_in_a_ring=0
    while (i6<=unique_rings){
	ringsize=split(Unique_Rings[i6],testring,":") - 1
	j6=1
	while (j6<=ringsize){
	    if (atom_to_check==testring[j6]){
		is_in_a_ring++}
	    j6++}
	i6++}
    return is_in_a_ring}
#
function find_linker_atoms(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    i=1
    split("",Linker_atoms,":")
    split("",Ring_atoms,":")
    number_of_linker_atoms=0
    number_of_ring_atoms=0
    while (i<=maxatom){
	number_of_connections=stripped_connected[i,0]
	if (stripped_connected[i,0]>0){
	    ring_check=is_in_ring(i, unique_rings)
	    if (ring_check==0){
		number_of_linker_atoms++
		Linker_atoms[number_of_linker_atoms]=i}
	    if (ring_check!=0){
		number_of_ring_atoms++
		Ring_atoms[number_of_ring_atoms]=i}}
	i++}
    number_of_ring_and_linker_atoms=number_of_ring_atoms":"number_of_linker_atoms
    return number_of_ring_and_linker_atoms}
#
function which_rings_atom_in(atom_to_check, unique_rings){
    i7=1
    set_of_rings=""
    while (i7<=unique_rings){
	ringsize=split(Unique_Rings[i7],testring,":") - 1
	j7=1
	while (j7<=ringsize){
	    if (atom_to_check==testring[j7]){
		set_of_rings=set_of_rings""i7":"}
	    j7++}
	i7++}
    return set_of_rings}
#
function are_atoms_in_same_ring_system(atom_1, atom_2, unique_rings){
    ring_set1=which_rings_atom_in(atom_1, unique_rings)
    ring_set2=which_rings_atom_in(atom_2, unique_rings)
    ring_set_size_1=split(ring_set1,RS1,":") - 1
    ring_set_size_2=split(ring_set2,RS2,":") - 1
    i8=1
    in_same_system=0
    while (i8<=ring_set_size_1){
	j8=1
	while (j8<=ring_set_size_2){
	    if (RS1[i8]==RS2[j8]){
		in_same_system=1}
	    j8++}
	i8++}
    return in_same_system}
#
function are_atoms_connected(atom_1, atom_2){
    number_connected_to_1=connected[atom_1,0]
    i9=1
    atoms_connected=0
    while (i9<=number_connected_to_1){
	if (atom_2==connected[atom_1,i9]){
	    atoms_connected=1
	    pointofconnection=i9}
	i9++}
    returnstring=atoms_connected":"pointofconnection
    return returnstring}
#
function find_ring_to_ring_linkers(number_of_ring_and_linker_atoms, unique_rings){
    split(number_of_ring_and_linker_atoms,MM,":")
    num_ring_atoms=MM[1]
    i10=1
    while (i10<=num_ring_atoms){
	j10=i10
	while (j10<=num_ring_atoms){
	    InSameSystem=are_atoms_in_same_ring_system(Ring_atoms[i10], Ring_atoms[j10], unique_rings)
	    AreConnectedAtomsString=are_atoms_connected(Ring_atoms[i10], Ring_atoms[j10])
	    split(AreConnectedAtomsString,AreConnectedAtoms,":")
	    if ((InSameSystem==0)&&(AreConnectedAtoms[1]==1)){
		Ring_to_Ring_link[Ring_atoms[i10]]=1
		Ring_to_Ring_link[Ring_atoms[j10]]=1}
	    j10++}
	i10++}}
#
function define_total_valance(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    i11=minatom
    while (i11<=maxatom){
	Total_Valance[i11]=connected[i11,0]
	i11++}}
#
function define_atomic_symbol(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    i12=minatom
    while (i12<=maxatom){
	current_type=atom_type[i12]
	split(current_type,ATS,".")
	atomic_rep=ATS[1]
	ll=length(atomic_rep)
	if (ll==1){Atomic_Symbol[i12]=atomic_rep}
	if (ll>1){
	    if ((substr(atomic_rep,1,1)=="C")){
		if (substr(atomic_rep,2,1)=="T"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="A"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="W"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="*"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="R"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="C"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="J"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="B"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="M"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="P"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="V"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="1"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="2"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="3"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)==""){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="L"){Atomic_Symbol[i12]="CL"}
		if (substr(atomic_rep,2,1)=="l"){Atomic_Symbol[i12]="CL"}
		if (substr(atomic_rep,2,1)=="a"){Atomic_Symbol[i12]="Ca"}
		if (substr(atomic_rep,2,1)=="u"){Atomic_Symbol[i12]="C"}
		if (substr(atomic_rep,2,1)=="s"){Atomic_Symbol[i12]="C"}}
	    if ((substr(atomic_rep,1,1)=="N")){
		if (substr(atomic_rep,2,1)=="*"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="2"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="3"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="A"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="B"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="C"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="D"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="J"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="L"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="s"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="u"){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)==""){Atomic_Symbol[i12]="N"}
		if (substr(atomic_rep,2,1)=="a"){Atomic_Symbol[i12]="IP"}}
	    if ((substr(atomic_rep,1,1)=="H")){
		if (substr(atomic_rep,2,1)=="1"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="2"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="3"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="4"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="5"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="A"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="C"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="P"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="S"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="W"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="O"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)==""){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="X"){Atomic_Symbol[i12]="H"}
		if (substr(atomic_rep,2,1)=="e"){Atomic_Symbol[i12]="He"}}
	    if ((substr(atomic_rep,1,1)=="O")){
		if (substr(atomic_rep,2,1)=="2"){Atomic_Symbol[i12]="O"}
		if (substr(atomic_rep,2,1)=="H"){Atomic_Symbol[i12]="O"}
		if (substr(atomic_rep,2,1)=="W"){Atomic_Symbol[i12]="O"}
		if (substr(atomic_rep,2,1)=="u"){Atomic_Symbol[i12]="O"}
		if (substr(atomic_rep,2,1)==""){Atomic_Symbol[i12]="O"}
		if (substr(atomic_rep,2,1)=="S"){Atomic_Symbol[i12]="O"}}
	    if ((substr(atomic_rep,1,1)=="S")){
		if (substr(atomic_rep,2,1)=="D"){Atomic_Symbol[i12]="S"}
		if (substr(atomic_rep,2,1)=="O"){Atomic_Symbol[i12]="S"}
		if (substr(atomic_rep,2,1)=="R"){Atomic_Symbol[i12]="S"}
		if (substr(atomic_rep,2,1)==""){Atomic_Symbol[i12]="S"}
		if (substr(atomic_rep,2,1)=="H"){Atomic_Symbol[i12]="S"}
		if (substr(atomic_rep,2,1)=="u"){Atomic_Symbol[i12]="S"}
		if (substr(atomic_rep,2,1)=="c"){Atomic_Symbol[i12]="Sc"}}
	    if (atomic_rep=="P"){Atomic_Symbol[i12]="P"}
	    if (atomic_rep=="Br"){Atomic_Symbol[i12]="Br"}
	    if (atomic_rep=="BR"){Atomic_Symbol[i12]="Br"}
	    if (atomic_rep=="Li"){Atomic_Symbol[i12]="Li"}
	    if (atomic_rep=="Rb"){Atomic_Symbol[i12]="Rb"}
	    if (atomic_rep=="Mg"){Atomic_Symbol[i12]="MG"}
	    if (atomic_rep=="MG"){Atomic_Symbol[i12]="MG"}
	    if (atomic_rep=="Fe"){Atomic_Symbol[i12]="FE"}
	    if (atomic_rep=="FE"){Atomic_Symbol[i12]="FE"}}
	i12++}}
#
function define_NonH_valance(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    i13=minatom
    while (i13<=maxatom){
	NonH_val=0
	number_connected_to_this_atom=connected[i13,0]
	j13=1
	while (j13<=number_connected_to_this_atom){
	    if ((Atomic_Symbol[connected[i13,j13]])!="H"){
		NonH_val++}
	    j13++}
	NonH_Valance[i13]=NonH_val
	i13++}}
#
function define_Pi_electrons(number_of_ring_and_linker_atoms){
    split(number_of_ring_and_linker_atoms,a4,":")
    number_of_ring_atoms=a4[1]
    i14=1
    while (i14<=number_of_ring_atoms){
	Piele=""
	atom_to_determine_pi=Ring_atoms[i14]
	AS=Atomic_Symbol[atom_to_determine_pi]
	TV=Total_Valance[atom_to_determine_pi]
	if ((AS=="C")&&(TV==2)){Piele=1}
	if ((AS=="C")&&(TV==4)){Piele=0}
	if ((AS=="N")&&(TV==2)){Piele=1}
	if ((AS=="N")&&(TV==3)){Piele=2}
	if ((AS=="N")&&(TV==4)){Piele=0}
	if ((AS=="O")&&(TV==2)){Piele=2}
	if ((AS=="S")&&(TV==2)){Piele=2}
	if ((AS=="S")&&(TV==4)){Piele=0}
	if ((AS=="B")&&(TV==3)){Piele=0}
	if ((AS=="P")){Piele=0}
	if (((AS=="C")||(AS=="S"))&&(TV==3)){
	    carbonyl_count=0
	    j14=1
	    while (j14<=TV){
		check_connected_atom=connected[atom_to_determine_pi,j14]
		if (((Atomic_Symbol[check_connected_atom]=="O")||(Atomic_Symbol[check_connected_atom]=="O"))&&(Total_Valance[check_connected_atom]==1)){
		    carbonyl_count++}
		j14++}
	    if (carbonyl_count>0){Piele=0}
	    else {Piele=1}}
	if (Piele!=""){
	    Pielectrons[atom_to_determine_pi]=Piele}
	else {print "FAILED TO ASSIGN PI electron count to: "AS" "atom_to_determine_pi" "TV}
	i14++}}
#
function get_pi_electrons_for_ring(specific_ring){
	workingring=Unique_Rings[specific_ring]
	ring_size=split(workingring,WR,":") - 1
	j15=1
	pi_elec_string=""
	while (j15<=ring_size){
	     pi_elec_string=pi_elec_string""Pielectrons[WR[j15]]":"
	    j15++}
	return pi_elec_string}
#
function does_ring_contain_SP3_atom(specific_ring){
	workingring=Unique_Rings[specific_ring]
	ring_size=split(workingring,WWR,":") - 1
	jj15=1
	SP3atom=0
	while (jj15<=ring_size){
             AS=Atomic_Symbol[WWR[jj15]]
	     TV=Total_Valance[WWR[jj15]]
	     if ((AS="C")&&(TV==4)){SP3atom=1}
	     if ((AS="S")&&(TV==4)){SP3atom=1}
	    jj15++}
	return SP3atom}
#
function get_pi_electron_pairs_from_pi_elec_string(pi_elec_string){
    ring_size=split(pi_elec_string,PS,":") - 1
    PS[(ring_size+1)]=PS[1]
    II17=1
    while (II17 <=2){
      upper_limit=(ring_size - (2-II17))
      i17=II17
      piset[II17]=""
      while (i17<=upper_limit){
	if ((PS[i17]+PS[(i17+1)])==0){
	    piset[II17]=piset[II17]"0:0:"
	    i17=i17+2}
	else if ((PS[i17]+PS[(i17+1)])==1){
	    if (PS[i17]==1){
		piset[II17]=piset[II17]"1:"
		i17=i17+1}
	    if (PS[i17+1]==1){
		piset[II17]=piset[II17]"0:"
		i17=i17+1}}
	else if ((PS[i17]+PS[(i17+1)])==2){
	    piset[II17]=piset[II17]"2:0:"
	    i17=i17+2}
	else if ((PS[i17]+PS[(i17+1)])==4){
	    piset[II17]=piset[II17]"2:2:"
	    i17=i17+2}
	else if ((PS[i17]+PS[(i17+1)])==3){
	    if ((PS[i17]==1)&&(PS[(i17+1)]==2)){
		piset[II17]=piset[II17]"1:2:"
		i17=i17+2}
	    if ((PS[i17]==2)&&(PS[(i17+1)]==1)){
		piset[II17]=piset[II17]"2:"
		i17=i17+1}}}
      if (i17==(ring_size+(II17-1))){
	piset[II17]=piset[II17]""PS[ring_size]":"}
      II17++}
    elec_pairs=piset[1]"/"piset[2]
    return elec_pairs}

#
function select_best_pi_elec_pairs(elec_pairs){
    split(elec_pairs,PI_SETS,"/")
    ring_size=split(PI_SETS[1],PS_1,":")
    split(PI_SETS[2],PS_2,":")
    k15=1
    number_of_double_zeros_1=0
    number_of_pairs_1=0
    number_of_singles_1=0
    total_number_of_electrons_1=0
    number_of_double_zeros_2=0
    number_of_pairs_2=0
    number_of_singles_2=0
    total_number_of_electrons_2=0
    Ring_type="DU"
    while (k15<=ring_size){
	total_number_of_electrons_1=total_number_of_electrons_1+PS_1[k15]
	total_number_of_electrons_2=total_number_of_electrons_2+PS_2[k15]
	if ((PS_1[k15]==0)&&(PS_1[(k15+1)])==0){number_of_double_zeros_1++}
	if ((PS_2[k15]==0)&&(PS_2[(k15+1)])==0){number_of_double_zeros_2++}
	if (PS_1[k15]==1){number_of_singles_1++}
	if (PS_2[k15]==1){number_of_singles_2++}
	if (PS_1[k15]==2){number_of_pairs_1++}
	if (PS_2[k15]==2){number_of_pairs_2++}
	k15++}
    if (number_of_singles_1 < number_of_singles_2){
        number_of_double_zeros=number_of_double_zeros_1
	number_of_singles=number_of_singles_1
	number_of_pairs=number_of_pairs_1
        BPS=PI_SETS[1]}
    else if (number_of_singles_1 > number_of_singles_2){
        number_of_double_zeros=number_of_double_zeros_2
	number_of_singles=number_of_singles_2
	number_of_pairs=number_of_pairs_2
        BPS=PI_SETS[2]}
    else if (number_of_singles_1 == number_of_singles_2){
        if (number_of_double_zeros_1 < number_of_double_zeros_2){
	  number_of_double_zeros=number_of_double_zeros_1
	  number_of_singles=number_of_singles_1
	  number_of_pairs=number_of_pairs_1
          BPS=PI_SETS[1]}
	else if (number_of_double_zeros_1 > number_of_double_zeros_2){
	  number_of_double_zeros=number_of_double_zeros_2
	  number_of_singles=number_of_singles_2
	  number_of_pairs=number_of_pairs_2
          BPS=PI_SETS[2]}
	else if (number_of_double_zeros_1 == number_of_double_zeros_2){
	  if (number_of_pairs_1 > number_of_pairs_2){
	    number_of_double_zeros=number_of_double_zeros_1
	    number_of_singles=number_of_singles_1
	    number_of_pairs=number_of_pairs_1
            BPS=PI_SETS[1]}
	  else if (number_of_pairs_1 < number_of_pairs_2){
	    number_of_double_zeros=number_of_double_zeros_2
	    number_of_singles=number_of_singles_2
	    number_of_pairs=number_of_pairs_2
            BPS=PI_SETS[2]}
	  else {
	    number_of_double_zeros=number_of_double_zeros_1
	    number_of_singles=number_of_singles_1
	    number_of_pairs=number_of_pairs_1
            BPS=PI_SETS[1]}}}
    best_pi_pairs=BPS"/"number_of_double_zeros"/"number_of_singles"/"number_of_pairs"/"total_number_of_electrons_1
    return best_pi_pairs}
#
function determine_ring_type(pi_elec_string){
    elec_pairs=get_pi_electron_pairs_from_pi_elec_string(pi_elec_string)
    best_pairs_data=select_best_pi_elec_pairs(elec_pairs)
    split(best_pairs_data,BPD,"/")
    ring_size=(split(pi_elec_string,JUNK,":")-1)
    Bpiset=BPD[1]
    number_of_double_zeros=BPD[2]
    number_of_singles=BPD[3]
    number_of_pairs=BPD[4]
    total_number_of_electrons=BPD[5]
    if (total_number_of_electrons==0){
	Ring_type="SR"ring_size}
#    else if ((number_of_singles==0)&&(total_number_of_electrons>0)&&(number_of_double_zeros==0)){
    else if ((number_of_singles==0)&&(total_number_of_electrons>0)){
	n_rawAR=((total_number_of_electrons-2)/4)     # Check for 4n+2 electrons, aromatic ring AR
	n_intAR=int((total_number_of_electrons-2)/4)  # Check for 4n+2 electrons, aromatic ring AR
        n_rawAAR=((total_number_of_electrons)/4)      # Check for 4n electrons, anti-aromatic ring AAR
        n_intAAR=int((total_number_of_electrons)/4)   # Check for 4n electrons, anti-aromatic ring AAR
	if (n_rawAR==n_intAR){
	    Ring_type="AR"ring_size}
	else if ((n_rawAAR==n_intAAR)&&(total_number_of_electrons>=4)){
	    Ring_type="AAR"ring_size}
	else {Ring_type="CR"ring_size}}
    else if ((number_of_singles>0)&&(number_of_pairs>0)){
	Ring_type="MR"ring_size}
    else {Ring_type="RR"ring_size}
    return Ring_type}
        
#
function fix_ring_charge(specific_ring, pi_elec_string){
	workingring=Unique_Rings[specific_ring]
	ring_size=split(workingring,WR,":") - 1
	ring_electrons=split(pi_elec_string,PIES,":") - 1
	total_elec=0
	ring_charge=0
	i16=1
	while (i16<=ring_size){
	    total_elec=total_elec+PIES[i16]
	    i16++}
	n_rawAR=((total_elec - 2)/4)     # Check for 4n+2 electrons, aromatic ring AR
	n_intAR=int((total_elec - 2)/4)  # Check for 4n+2 electrons, aromatic ring AR
	n_rawAAR=((total_elec)/4)      # Check for 4n electrons, anti-aromatic ring AAR
	n_intAAR=int((total_elec)/4)   # Check for 4n electrons, anti-aromatic ring AAR
	n_rawAR_an=((total_elec - 1)/4)     # Check for 4n+2 electrons, anionic aromatic ring AR
	n_intAR_an=int((total_elec - 1)/4)  # Check for 4n+2 electrons, anionic aromatic ring AR
	n_rawAR_cat=((total_elec - 3)/4)     # Check for 4n+2 electrons, cationic aromatic ring AR
	n_intAR_cat=int((total_elec - 3)/4)  # Check for 4n+2 electrons, cationic aromatic ring AR
	if (n_rawAR==n_intAR){
	    ring_charge=0}
	else if (n_rawAAR==n_intAAR){
	    ring_charge=0}
	else if (n_rawAR_an==n_intAR_an){
	    ring_charge=-1}
	else if (n_rawAR_cat==n_intAR_cat){
	    ring_charge=1}
	if (ring_charge==-1){
	    flag=0
	    i16=1
	    while (i16<=ring_size){
		if ((Atomic_Symbol[WR[i16]]=="N")&&(PIES[i16]==1)&&(flag==0)){
		    flag=1
		    PIES[i16]=2
		    FCharge[WR[i16]]=-1}
		i16++}}
	if (ring_charge==1){
	    flag=0
	    i16=1
	    while (i16<=ring_size){
		if ((Atomic_Symbol[WR[i16]]=="N")&&(PIES[i16]==2)&&(flag==0)){
		    flag=1
		    PIES[i16]=1
		    FCharge[WR[i16]]=1}
		i16++}}
	new_pi_elec_string=""
	i16=1
	while (i16<=ring_size){
	    new_pi_elec_string=new_pi_elec_string""PIES[i16]":"
	    i16++}
	return new_pi_elec_string}
#
#
function get_ring_types(unique_rings){
    i18=1
    while (i18<=unique_rings){
        SP3atom=0
	pi_elec_string=get_pi_electrons_for_ring(i18)
	pi_elec_string=fix_ring_charge(i18, pi_elec_string)
	PiString[i18]=pi_elec_string
	Ring_type=determine_ring_type(pi_elec_string)
	SP3atom=does_ring_contain_SP3_atom(i18) 
	if (SP3atom==0){Ring_Type[i18]=Ring_type}
	else {
	  if (substr(Ring_type,1,2)=="AR"){gsub("AR","MR",Ring_type)}
	  if (substr(Ring_type,1,3)=="AAR"){gsub("AAR","MR",Ring_type)}
	  if (substr(Ring_type,1,2)=="CR"){gsub("CR","MR",Ring_type)}
	  Ring_Type[i18]=Ring_type}
	i18++}}
#
function define_ring_types_for_atoms(unique_rings){
    i19=1
    while (i19<=unique_rings){
	temp_ring=Unique_Rings[i19]
	ring_type=Ring_Type[i19]
	ring_size=split(temp_ring,TR,":") - 1
	j19=1
	while (j19<=ring_size){
	    specific_atom=TR[j19]
	    Ring_Type_Atom[specific_atom]=Ring_Type_Atom[specific_atom]""ring_type":"
	    j19++}
	i19++}}
#
function remove_ring_atoms(unique_rings,minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    i20=minatom
    while (i20<=maxatom){
	j20=1
	check_atom=is_in_ring(i20, unique_rings)
	if ((check_atom==0)&&(stripped_connected[i20,0]!=0)){
	    linker_connected[i20,0]=stripped_connected[i20,0]
	    while (j20<=stripped_connected[i20,0]){
		conn_atom=stripped_connected[i20,j20]
		check_atom1=is_in_ring(conn_atom, unique_rings)
		if (check_atom1!=0){
		    linker_connected[i20,0]=linker_connected[i20,0] - 1}
		else {linker_connected[i20,j20]=stripped_connected[i20,j20]}
		j20++}}
	i20++}}
#
function define_linker_atoms(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    flag=0
    position_counter=1
    while (flag==0){
	flag=1
	i21=minatom
	terminal_set=""
	while (i21<=maxatom){
	    num_connected=linker_connected[i21,0]
	    if (num_connected==1){
		Linker_Set[i21]="L."position_counter
		terminal_set=terminal_set""i21":"
	        linker_connected[i21,0]=0}
	    if (num_connected>1){
		flag=0}
	    i21++}
	i21=minatom
	while (i21<=maxatom){
	    found_term=split(terminal_set,TS,":")
	    num_connected=linker_connected[i21,0]
	    j21=1
	    while (j21<=found_term){
		k21=1
		while (k21<=num_connected){
		    if (linker_connected[i21,k21]==TS[j21]){
			linker_connected[i21,0]=linker_connected[i21,0]-1
			temp_pos=k21+1
			while (temp_pos<=num_connected){
			    linker_connected[i21,(temp_pos-1)]=linker_connected[i21,(temp_pos)]
			    temp_pos++}}
		    k21++}
		j21++}
	    i21++}
	position_counter++}}
#
function Bond_String_from_Pi_string(pistring){
    Npi=split(pistring, PiStrg, ":") - 1
    i22=1
    flag1=0
    BondOrderString=""
    while (i22 <= Npi){
	if (PiStrg[i22]!=1){
	    flag1=1}
	i22++}
    if (flag1==0){
	i22=1
	flip_flag=0
	while (i22 <= Npi){
	    if (flip_flag == 0){
		BondOrderString=BondOrderString"2:"
		flip_flag=1}
	    else {
		BondOrderString=BondOrderString"1:"
		flip_flag=0}
	    i22++}}
    else {
	flag2=0
	i22=1
	while (i22 < Npi){
	    if ((PiStrg[i22]==1)&&(PiStrg[(i22+1)]==1)&&(flag2==0)){
		BondOrderString=BondOrderString"2:"
		flag2=1}
	    else if ((PiStrg[i22]==1)&&(PiStrg[(i22+1)]==1)&&(flag2==1)){
		BondOrderString=BondOrderString"1:"
		flag2=0}
	    else {
		BondOrderString=BondOrderString"1:"
		flag2=0}
	    i22++}
	if ((PiStrg[1]==1)&&(PiStrg[(Npi)]==1)&&(flag2==0)){
	    BondOrderString=BondOrderString"2:"
	    flag2=1}
	else if ((PiStrg[1]==1)&&(PiStrg[(Npi)]==1)&&(flag2==1)){
	    BondOrderString=BondOrderString"1:"
	    flag2=0}
	else {
	    BondOrderString=BondOrderString"1:"
	    flag2=0}}
    return BondOrderString}
#
function Get_all_desc(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    i13=minatom
    while (i13<=maxatom){
	if (Ring_to_Ring_link[i13]==1){Linker_Set[i13]="L.0"}
	desc[i13]=""
	desc[i13]=desc[i13]""Atomic_Symbol[i13]":V."Total_Valance[i13]":v."NonH_Valance[i13]":"Ring_Type_Atom[i13]":"chain_atom[i13]":"Linker_Set[i13]
	gsub("::",":",desc[i13])
	gsub("::",":",desc[i13])
	i13++}}
#
function Rebuild_desc(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    L13=minatom
    while (L13<=maxatom){
      desc_string=main[L13,2]
      desc_length=split(desc_string,des1,":")
      LL13=1
      new_desc=""
      while (LL13 <=desc_length){
	LLL13=LL13+1
	while (LLL13 <= desc_length){
	  if (des1[LL13]==des1[LLL13]){
	    if ((match(des1[LL13],"RR"))||(match(des1[LL13],"AR"))||(match(des1[LL13],"CR"))||(match(des1[LL13],"MR"))||(match(des1[LL13],"SR"))){
	      des1[LL13]=des1[LL13]}
	    else {
	      des1[LLL13]=""}}
	  LLL13++}
	if (des1[LL13]!=""){
	  new_desc=new_desc":"des1[LL13]}
	LL13++}
      desc[L13]=new_desc
      L13++}}
#
function Prep_System(Mol2_File){
    reset_all_arrays(0)
    minmax=read_single_mol2(Mol2_File)
    split(minmax,a1,":")
    minmax=a1[1]":"a1[2]
    pos1=a1[3]
    pos4=a1[4]
    var=terminal_atoms(minmax)
    maxloop=make_random_path(minmax)
    num_rings=find_all_rings_in_one_pass(minmax, maxloop) 
    clean_ring_count=remove_straight_dupe_rings(0)
    unique_rings=remove_complex_dupe_rings(clean_ring_count)
    number_of_ring_and_linker_atoms=find_linker_atoms(minmax)
    find_ring_to_ring_linkers(number_of_ring_and_linker_atoms, unique_rings)
    define_total_valance(minmax)
    define_atomic_symbol(minmax)
    define_NonH_valance(minmax)
    define_Pi_electrons(number_of_ring_and_linker_atoms)
    get_ring_types(unique_rings)
    define_ring_types_for_atoms(unique_rings)
    remove_ring_atoms(unique_rings,minmax)
    define_linker_atoms(minmax)
    Get_all_desc(minmax)
    Return_string=""
    Return_string=Return_string""minmax":"unique_rings":"pos1":"pos4
    Build_main(minmax)
    return Return_string}
#
function Build_main(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    i31=minatom
    while (i31<=maxatom){
	main[i31,0]=4
	main[i31,1]=i31
	main[i31,2]=desc[i31]
	j31=1
	num_connected=connected[i31,0]
	adjacent_atoms=""
	connected_list=""
	while (j31<=num_connected){
	    adjacent_atoms=adjacent_atoms""desc[connected[i31,j31]]"/"
	    connected_list=connected_list""connected[i31,j31]":"
	    j31++}
	main[i31,3]=adjacent_atoms
	main[i31,4]=connected_list
	atom_types[i31]="Du"
	formal_charges[i31]=0.0
	i31++}}
#
function Rebuild_main(minmax){
    split(minmax,a1,":")
    minatom=a1[1]
    maxatom=a1[2]
    L31=minatom
    while (L31<=maxatom){
	main[L31,0]=4
	main[L31,1]=L31
	main[L31,2]=desc[L31]
	K31=1
	num_connected=connected[L31,0]
	adjacent_atoms=""
	connected_list=""
	while (K31<=num_connected){
	    adjacent_atoms=adjacent_atoms""desc[connected[L31,K31]]"/"
	    K31++}
	main[L31,3]=adjacent_atoms
	L31++}}
#
# Now the Pattern Matching Functions:
#
function match_wild(desc_with_wild, primary_desc){
    num_characters_wild=(split(desc_with_wild,DWW,""))
    num_characters_prim=(split(primary_desc,PW,""))
    num_matching_characters=0
    they_match=0
    i0=1
    while (i0<=num_characters_wild){
	if (DWW[i0]==PW[i0]){
	    num_matching_characters++}
	i0++}
    if (DWW[num_characters_wild]=="?"){
	if (num_matching_characters==(num_characters_wild-1)){
	    they_match=1}}
    else {if ((num_matching_characters==num_characters_wild)&&(num_matching_characters==num_characters_prim)){
	    they_match=1}}
    return they_match}
#
function match_wild_withrings(desc_with_wild, primary_desc, current_atom, pattern_list){
    num_characters_wild=(split(desc_with_wild,DWW,""))
    num_characters_prim=(split(primary_desc,PW,""))
    num_matching_characters=0
    they_match=0
    i0=1
    if (substr(desc_with_wild,1,1)=="@"){ # Ring closuer token}
	atom_number_to_link=desc_with_wild
	sub("@","",atom_number_to_link)
	k11=split(pattern_list,PL,":")
	matching_atom=PL[atom_number_to_link]
	kk11=match_connect(current_atom, matching_atom)
	if (kk11=1){they_match=1}}
    else{
    while (i0<=num_characters_wild){
	if (DWW[i0]==PW[i0]){
	    num_matching_characters++}
	i0++}
    if (DWW[num_characters_wild]=="?"){
	if (num_matching_characters==(num_characters_wild-1)){
	    they_match=1}}
    else {if ((num_matching_characters==num_characters_wild)&&(num_matching_characters==num_characters_prim)){
	    they_match=1}}}
    return they_match}
#
function Match_Count(desc, target_string, current_atom, pattern_list){  # Returns the number of times desc appears in target string
    num_desc=split(target_string,a1,":")
    i1=1
    count=0
    while (i1<=num_desc){
	k1=match_wild_withrings(desc, a1[i1], current_atom, pattern_list)
	if (k1==1){count++}
	i1++}
    return count}
#
function match_desc(desc, target_string, current_atom, pattern_list){ # Returns 1 if the descriptor (desc) is found in the target string, else returns 0
    num_of_ors=split(desc,a2,"|")
    i2=1
    num_required1=1
    match_made1=0
    while (i2<=num_of_ors){
	repeats=substr(a2[i2],1,1)
	if (repeats ~ /[0-9]/){
	    num_required1=repeats
	    actual_match=a2[i2]
	    sub(repeats,"",actual_match)}
	else {actual_match=a2[i2]
	    num_required1=-1}
	k1=Match_Count(actual_match, target_string, current_atom, pattern_list)
	if ((num_required1==0)&&(k1==0)){match_made1=1}
	else if ((num_required1<0)&&(k1>=1)){match_made1=1}
	else if ((num_required1>=1)&&(k1==num_required1)){match_made1=1}
	i2++}
    return match_made1}
#
function match_strings(pattern, target_string, current_atom, pattern_list){ # Returns 1 if the pattern is found in the target string, else returns 0
    num_of_desc3=split(pattern,a3,":")
    i3=1
    total_match3=0
    match_made3=0
    while (i3<=num_of_desc3){
	k3=match_desc(a3[i3], target_string, current_atom, pattern_list)
	if (k3==1){total_match3++}
	i3++}
    if (total_match3==num_of_desc3){
	match_made3=1}
    else if(total_match3!=num_of_desc3){
	match_made3=0}
    return match_made3}
#
function match_nestedpattern(pattern, target_atom, pattern_list){
    matchnested=1
    sub_match=0
    num_sub_patterns=(split(pattern,NSP,"/") - 1)
    num_of_connected_atoms=(split(main[target_atom,4],CA,":") - 1)
    if (num_sub_patterns>(num_of_connected_atoms)){matchnested=0}
    if (num_sub_patterns==0){
	k9=match_strings(NSP[1], main[target_atom,2], target_atom, pattern_list)
	if (k9==0){
	    matchnested=0}
	else if (k9==1){
	    matchnested=1}}
    if ((num_sub_patterns>0)&&(num_sub_patterns<=num_of_connected_atoms)){
	k9=match_strings(NSP[1], main[target_atom,2], target_atom, pattern_list)
	if (k9==0){
	    matchnested=0}
	else if (k9==1){
	    matchnested=1
	    i9=1
	    while (i9<=num_of_connected_atoms){
		j9=1
		while (j9<=(num_sub_patterns)){
		    kk9=match_strings(NSP[j9+1], main[CA[i9],2])
		    if (kk9==1){
			CA[i9]="YY"
			NSP[j9+1]="XX"
			sub_match++}
		    j9++}
		i9++}
	if (sub_match>=(num_sub_patterns)){
	    matchnested=1}
	else{
	    matchnested=0}}}
    if (matchnested==1){
	match_made_9=1}
    else{
	match_made_9=0}
    return match_made_9}
#
function match_primary(pattern, number_of_atoms, pattern_list){
    i4=1
    match_set4=""
    while (i4<=number_of_atoms){
	test_string=main[i4,2]
	k4=match_nestedpattern(pattern, i4, pattern_list)
	if (k4==1){match_set4=match_set4""i4":"}
	i4++}
    return match_set4}
#
function match_full_pattern(test_desc, num_atoms){
    num_atoms_in_pattern=split(test_desc,desc_set,"-")
    i6=1
    list_3=""
    while (i6<=num_atoms_in_pattern){
	if (i6==1){
	    list_1=match_primary(desc_set[i6], num_atoms, list_3)
	    gsub(":",":/",list_1)}
	if (i6>1){
	    num_matches=(split(list_1,list_2,"/") - 1)
	    j6=1
	    new_Match_Count=0
	    while (j6<=num_matches){
		ee=(split(list_2[j6],temp,":") - 1)
		atom_number=temp[ee]
		gsub(":","",atom_number)
		next_set=main[atom_number,4]
		num_in_next_step=(split(next_set,zzz,":") -1)
		zz=1
		while (zz<=num_in_next_step){
		    k=match_nestedpattern(desc_set[i6], zzz[zz], list_2[j6])
		if (k==1){
		    new_Match_Count++
		    num_matches_so_far=split(list_2[j6],ZZZZ,":")
		    if (num_matches_so_far==i6){
			list_3=list_3""list_2[j6]""zzz[zz]":/"}}
		zz++}
		j6++}
	    list_1=list_3}
	i6++}
    num_total_and_partial_matches=(split(list_1,list_2,"/") - 1)
    i6=1
    list_3=""
    while (i6<=num_total_and_partial_matches){
	k=(split(list_2[i6],HH,":") - 1)
	if (k==num_atoms_in_pattern){
	    repeated_atom=0
	    p6=1
	    while (p6<num_atoms_in_pattern){
		jj6=p6+1
		while (jj6<=num_atoms_in_pattern){
		    if (HH[p6]==HH[jj6]){repeated_atom=1}
		    jj6++}
		p6++}
	    if (repeated_atom==0){
		list_3=list_3""list_2[i6]"/"}}
	i6++}
    return list_3}
#
function modify_desc(match_list, new_desc_list, minmax){
    number_of_match_sets=(split(match_list,NMS,"/") - 1)
    number_of_new_desc=(split(new_desc_list,NDL,":") - 1)
    i7=1
    while (i7<=number_of_match_sets){
	speific_match_set=(split(NMS[i7],SMS,":") - 1)
	ii7=1
	while (ii7<=speific_match_set){
	    k7=match_wild(NDL[ii7], SMS[ii7])
	    if (k7==0){
		if ((NDL[ii7]!=".")&&(NDL[ii7]!="")){
		    current_desc=main[SMS[ii7],2]
		    if (substr(current_desc,length(current_desc),1)==":"){
		      main[SMS[ii7],2]=main[SMS[ii7],2]""NDL[ii7]}
		    else {
		      main[SMS[ii7],2]=main[SMS[ii7],2]":"NDL[ii7]}}}
	    ii7++}
	i7++}}
#
function modify_atom_types(match_list, new_AT_list){
    number_of_match_sets=(split(match_list,NMS,"/") - 1)
    number_of_new_AT=(split(new_AT_list,NATL,":") - 1)
    i8=1
    while (i8<=number_of_match_sets){
	speific_match_set=(split(NMS[i8],SMS,":") - 1)
	ii8=1
	while (ii8<=speific_match_set){
	    if ((NATL[ii8]!=".")&&(NATL[ii8]!="")){
		atom_types[SMS[ii8]]=NATL[ii8]}
	    ii8++}
	i8++}}
#
function match_connect(primary_atom, secondary_atom){
    num_connected=(split(main[primary_atom,4],CL,":")-1)
    i12=1
    conn_match_made=0
    while (i12<=num_connected){
	if (secondary_atom==CL[i12]){
	    conn_match_made=1}
	i12++}
    return conn_match_made}
#
function modify_formal_charges(match_list, new_Fcharge_list){
    number_of_match_sets=(split(match_list,NMS,"/") - 1)
    number_of_new_FC=(split(new_Fcharge_list,NFCL,":") - 1)
    i13=1
    while (i13<=number_of_match_sets){
	speific_match_set=(split(NMS[i13],SMS,":") - 1)
	ii13=1
	while (ii13<=speific_match_set){
	    if ((NFCL[ii13]!=".")&&(NFCL[ii13]!="")){
		formal_charges[SMS[ii13]]=NFCL[ii13]}
	    ii13++}
	i13++}}
#
#
# Generate output
#
#
function print_atom_numbers(match_list){
    number_of_match_sets=(split(match_list,NMS,"/") - 1)
    i14=1
    while (i14<=number_of_match_sets){
	speific_match_set=(split(NMS[i14],SMS,":") - 1)
	ii14=1
	while (ii14<=speific_match_set){
	    print SMS[ii14] 
	    ii14++}
	print " "
	i14++}}
#
function print_atom_centered_charges(test_desc, match_list){
    number_of_match_sets=(split(match_list,NMS,"/") - 1)
    i14=1
    while (i14<=number_of_match_sets){
	speific_match_set=(split(NMS[i14],SMS,":") - 1)
	ii14=1
	printf test_desc" "
	while (ii14<=speific_match_set){
	    printf atomic_charge[SMS[ii14]]" " 
	    ii14++}
	print " "
	i14++}}
#
function out_mol2(minatom,maxatom,pos1,pos4,outtype){
    i35=minatom
    j35=1
    print "@<TRIPOS>MOLECULE"
    while (j35<pos1){
	print header_set[j35]
	j35++}
    print "@<TRIPOS>ATOM"
    while (i35<=maxatom){
      if (outtype=="atom_type"){
	print i35, atom_label[i35], coordinates[i35,1], coordinates[i35,2], coordinates[i35,3], atom_types[i35], res_num[i35], res_name[i35], atomic_charge[i35], atomic_tag[i35]}
      else if  (outtype=="atom_desc"){
	print i35, atom_label[i35], coordinates[i35,1], coordinates[i35,2], coordinates[i35,3], main[i35,2], res_num[i35], res_name[i35], atomic_charge[i35], atomic_tag[i35]}
      else if  (outtype=="atom_chrg"){
	print i35, atom_label[i35], coordinates[i35,1], coordinates[i35,2], coordinates[i35,3], atomic_charge[i35,2], res_num[i35], res_name[i35], atomic_charge[i35], atomic_tag[i35]}
      else {
	print i35, atom_label[i35], coordinates[i35,1], coordinates[i35,2], coordinates[i35,3], atom_types[i35], res_num[i35], res_name[i35], atomic_charge[i35], atomic_tag[i35]}
      i35++}
    print "@<TRIPOS>BOND"
    k35=1
    i35=minatom
    while (i35<=maxatom){
	j35=1
	while (j35<=connected[i35,0]){
	    atom1=i35
	    atom2=connected[i35,j35]
	    bond_order=bondorder[i35,j35]
	    if (bond_order==1){bond_order=1}
	    else if  (bond_order==2){bond_order=2}
	    else if  (bond_order==3){bond_order=3}
	    else if  (bond_order==1.5){bond_order="ar"}
	    else if  (bond_order==-1){bond_order="du"}
	    else if  (bond_order==-2){bond_order="uc"}
	    else {bond_order="nc"}
	    if (atom2!=0){
		print k35, atom1, atom2, bond_order
	    k35++
	    l35=1
	    while (l35<=connected[atom2,0]){
		if (connected[atom2,l35]==atom1){
		    connected[atom2,l35]=0}
		l35++}}
	    j35++}
	i35++}
    print "@<TRIPOS>SUBSTRUCTURE"
#    print " "
    i35=minatom
    while (i35<=pos4){
	print substructure[i35,1], substructure[i35,2], substructure[i35,3], substructure[i35,4], substructure[i35,5], substructure[i35,6], substructure[i35,7], substructure[i35,8]
	i35++}}
#
function find_core_only_atoms(main_ring_list,all_atom_list){
    number_atoms_in_main=split(main_ring_list,main_ring_list_array,":") - 1
    number_atoms_in_system=split(all_atom_list,all_atom_list_array,":") - 1
    core_only_atoms=""
    core_only_exist=0
    i37=1
    while (i37 <= number_atoms_in_main){
	main_ring_atoms[main_ring_list_array[i37]]=1
	i37++}
    i37=1
    while (i37 <= number_atoms_in_system){
	if (all_atom_list_array[i37] in main_ring_atoms){
	    donothing=1}
	else{
	    core_only_exist=1
	    core_only_atoms=core_only_atoms""all_atom_list_array[i37]":"}
	i37++}
    if (core_only_exist == 1){
	return core_only_atoms}
    else {
	return "0"}}
#
function find_highest_occurances_in_list(atom_list,list_of_occurances){
    number_of_atoms=split(list_of_occurances,list_of_atom_count_pairs,"/") - 1
    number_of_atoms_in_specific_list=split(atom_list,speicific_atoms,":") - 1
    max=0
    return_string=""
    i36=1
    while (i36 <= number_of_atoms){
	split(list_of_atom_count_pairs[i36],pair,":")
	occurance_count[pair[1]]=pair[2]
	i36++}
    i36=1
    while (i36 <= number_of_atoms_in_specific_list){
	if (occurance_count[speicific_atoms[i36]] > max){
	    count=1
	    max=occurance_count[speicific_atoms[i36]]
	    return_string=speicific_atoms[i36]":"}
	else if (occurance_count[speicific_atoms[i36]] == max){
	    count++
	    return_string=return_string""speicific_atoms[i36]":"}
	i36++}
    if (count < number_of_atoms_in_specific_list){
	return return_string}
    else{
	return "0"}}
#
function find_number_of_times_entries_occurs_in_list(full_atom_list_repeats){
    atom_list=unique_only(full_atom_list_repeats)
    number_unique_entries=split(atom_list,list_of_atoms,":") -1
    number_nonunique_entries=split(full_atom_list_repeats,list_of_atoms_with_repeats,":") -1
    i35=1
    split("",count_list,":")
    while (i35 <= number_nonunique_entries){
	if (list_of_atoms_with_repeats[i35] in count_list){
	    count_list[list_of_atoms_with_repeats[i35]] = count_list[list_of_atoms_with_repeats[i35]] +1}
	else {
	    count_list[list_of_atoms_with_repeats[i35]] = 1}
	i35++}
    return_string=""
    i35=1
    while (i35 <= number_unique_entries){
	return_string=return_string""list_of_atoms[i35]":"count_list[list_of_atoms[i35]]"/"
	i35++}
    return return_string}
#
function find_largest_ring(unique_AR_ring_string,list_of_occurances){
  i33=1
    number_of_rings_in_string=split(unique_AR_ring_string,rings_in_string,":") - 1
    number_of_atoms_in_system=split(list_of_occurances,Apairs,"/") - 1
    while (i33 <= number_of_atoms_in_system){
	split(Apairs[i33],pair,":")
	atom_score[pair[1]]=pair[2]
	i33++}
    maxsize=0
    minscore=9999999999999
    i33=1
    while (i33 <= number_of_rings_in_string){
      ring_length=split(Unique_Rings[rings_in_string[i33]],dummy_list,":") - 1 
      j33=1
      rings_in_string_score=0
      while (j33 <= ring_length){
	  rings_in_string_score=rings_in_string_score + atom_score[dummy_list[j33]]
	  j33++}
      if ((ring_length > maxsize)){
	  minscore=9999999999999
	  if (rings_in_string_score < minscore){
	      maxsize=ring_length
	      minscore=rings_in_string_score
	      largest_ring=rings_in_string[i33]}}
      if ((ring_length == maxsize)){
	  if (rings_in_string_score < minscore){
	      maxsize=ring_length
	      minscore=rings_in_string_score
	      largest_ring=rings_in_string[i33]}}
      i33++}
    return largest_ring}
#
function common_entries(ring1, ring2){
    Natoms1=split(ring1,Ring1Set,":") - 1
    Natoms2=split(ring2,Ring2Set,":") - 1
    num_common_atoms=0
    i23=1
    while (i23 <= Natoms1){
	j23=1
	while (j23 <= Natoms2){
	    if (Ring1Set[i23]==Ring2Set[j23]){
		num_common_atoms++}
	    j23++}
	i23++}
    return num_common_atoms}
#
function find_connected_rings(test_ring){
    i25=1
    listofrings=test_ring":"
    while (i25 <= unique_rings){
	if (RingOverlaps[test_ring,i25]>=2){
	    listofrings=listofrings i25":"}
	i25++}
    return listofrings}
#
function find_connected_ARrings(test_ring){
    ii25=1
    listofrings=test_ring":"
    while (ii25 <= unique_rings){
	if (ARRingOverlaps[test_ring,ii25]>=2){
	    listofrings=listofrings ii25":"}
	ii25++}
    return listofrings}
#
function unique_only(string){
    n26=split(string,string_list,":") - 1
    split("",counterlist,":")
    new_string=""
    i26=1
    while (i26 <= n26){
	counterlist[(string_list[i26])]=0
	i26++}
    for (elementofstring in counterlist){
	if (elementofstring!=""){
	    new_string=new_string elementofstring ":"}}
    return new_string}
#
function modify_atomic_charges(match_list, chrg_mods){
    number_of_match_sets=(split(match_list,NMS,"/") - 1)
    number_of_chrg_mods=(split(chrg_mods,NATL,":") - 1)
    i27=1
    while (i27<=number_of_match_sets){
	speific_match_set=(split(NMS[i27],SMS,":") - 1)
	ii27=1
	while (ii27<=speific_match_set){
	    if ((NATL[ii27]!=".")&&(NATL[ii27]!="")){
		atomic_charge[SMS[ii27]]=atomic_charge[SMS[ii27]]+NATL[ii27]}
	    ii27++}
	i27++}}
#
function modify_atomic_label(match_list, label_mods){
    number_of_match_sets=(split(match_list,NMS,"/") - 1)
    number_of_label_mods=(split(label_mods,NATL,":") - 1)
    i28=1
    while (i28<=number_of_match_sets){
	speific_match_set=(split(NMS[i28],SMS,":") - 1)
	ii28=1
	while (ii28<=speific_match_set){
	    if ((NATL[ii28]!=".")&&(NATL[ii28]!="")){
		atom_label[SMS[ii28]]=NATL[ii28]}
	    ii28++}
	i28++}}
#
#
BEGIN{setout=Prep_System(mol_file)
    split(setout,a1,":")
    minmax=a1[1]":"a1[2]
    minatom=a1[1]
    maxatom=a1[2]
    unique_rings=a1[3]
    pos1=a1[4]
    pos4=a1[5]
    MOL2_OUT=1}
#
{test_desc=$1
 work=$2
 mod_list=$3
 if (substr($1,1,1)!="#"){
   if (work=="desc"){
     match_list=match_full_pattern(test_desc, maxatom)
     modify_desc(match_list, mod_list, minmax)
     Rebuild_desc(minmax)
     Rebuild_main(minmax)}
   else if (work=="AT"){
     match_list=match_full_pattern(test_desc, maxatom)
     modify_atom_types(match_list, mod_list)}
   else if (work=="FC"){
     match_list=match_full_pattern(test_desc, maxatom)
     modify_formal_charges(match_list, mod_list)}
   else if (work=="PR"){
     match_list=match_full_pattern(test_desc, maxatom)
     print_atom_numbers(match_list)
     MOL2_OUT=1}
   else if (work=="PRCHRGS"){
     match_list=match_full_pattern(test_desc, maxatom)
     print_atom_centered_charges(test_desc, match_list)
     MOL2_OUT=0}
   else if (work=="MC"){
     match_list=match_full_pattern(test_desc, maxatom)
     modify_atomic_charges(match_list, mod_list)}
   else if (work=="label"){
     match_list=match_full_pattern(test_desc, maxatom)
     modify_atomic_label(match_list, mod_list)}
   else if (work=="outmol2"){
       out_mol2(minatom,maxatom,pos1,pos4,mod_list)}}}
#END{if (MOL2_OUT==1){out_mol2(minatom,maxatom,pos1,pos4)}}
#  LLL=1;while (LLL <=maxatom){print "ONE ",main[LLL,1], "TWO ",main[LLL,2], "THREE ", main[LLL,3], "FOUR ", main[LLL,4];LLL++}}