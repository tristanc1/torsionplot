#
# Ramachandran plot generator
#
# $Id: torsionplot.tcl,v 1.9 2006/08/16 04:05:45 johns Exp $
# 
# Version history:
#

package provide torsionplot 1.0

namespace eval ::TorPlot:: {
  #namespace export torsionplot

  variable torsel ""	;# atom selections for current molecule 
  variable molid ""	;# molid of current molecule
  variable seltext all  ;# selection text in entry box
#   variable hlresid -1   ;# highlight resid
  variable w 		;# handle to window
#   variable data		;# you know, data
#   variable box		;# all those little boxes I draw
  variable lastsavedfile torsionplot.tga
#  variable moldispstat "";# store displayed status before hiding for mk3drama
  variable tempdir "$env(TMPDIR)"   ;# directory to which all temporary files will be written
  variable datadir "$env(TORPLOTDIR)"         ;# directory containing all the probability maps
  variable torfilename [file join $tempdir torplot_temp.pdb]
  variable torplotmol "-1"                    ;# molid handle for the molecule that will contain the data points
  variable torplotall ""
  
  
  # variable connecting each case to a filename containing the probability map for that case
  # format is:
  # 0. Case Name
  # 1. atomselect text
  # 2. map file
  # 3. is3D (boolean: 0 for most cases, 1 for glycan 1-6 links)
  # 4. list of cutoffs to use for isosurface representations
  # 5. Boolean: 1 for protein torsions, 0 for everything else
  # 6. (only needed for non-protein torsions)  list of atoms from second residue contributing to phi/psi/omega
  #       in order from linking O or N
  # Glycan selections should uniquely identify the residue containing the C1 atom in the glycosidic bond.
  # The rest will be worked out below.
  
  variable tortypes {
      {{Protein General Case} {protein and not resname GLY PRO and not (same residue as (name C and within 1.5 of (resname PRO and name N)))} rama-general.dx 0 {0.0005 0.02} 1} 
      {Glycine {resname GLY} rama-gly.dx 0 {0.002 0.02} 1}
      {{Preceding Proline} {protein and not resname GLY PRO and (same residue as (name C and within 1.5 of (resname PRO and name N)))} rama-prepro.dx 0 {0.002 0.02} 1 } 
      {Proline {resname PRO} rama-pro.dx 0 {0.002 0.02} 1}
      {{Glycan N-link (fucosylated)} {resname BGLN and same residue as (name C1 and within 2 of name ND2) and same residue as (name O6 and within 2 of (resname AFUC and name C1))} NGLB_FUC.dx 0 {0.01 0.05} 0 {ND2 CG CB}}
      {{Glycan N-link (no core fucose)} {resname BGLN and same residue as (name C1 and within 2 of name ND2) and same residue as name HO6} NGLB_no_FUC.dx 0 {0.01 0.05} 0 {ND2 CG CB}}
      {{Man (alpha-1,3) Man} {resname AMAN and same residue as (name C1 and within 2 of (resname AMAN BMAN and name O3))} MAN_13_AMAN.dx 0 {0.01 0.05} 0 {O3 C3 C2}}
      {{alpha-Man (alpha-1,6) Man} {resname AMAN and same residue as (name C1 and within 2 of (resname AMAN and name O6))} AMAN_16_AMAN.dx 1 {0.01 0.05} 0 {O6 C6 C5 O5}}
      {{beta-Man (alpha-1,6) Man} {resname AMAN and same residue as (name C1 and within 2 of (resname BMAN and name O6))} BMAN_16_AMAN.dx 1 {0.01 0.05} 0 {O6 C6 C5 O5}}
      {{(GlcNAc or alpha-Man) (beta-1,4) GlcNAc} {resname BGLN and same residue as (name C1 and within 2 of (resname BGLN AMAN and name O4))} {BGLN_or_AMAN_14_BGLN.dx} 0 {0.01 0.05} 0 {O4 C4 C3}}
      {{GlcNAc (alpha-1,6) Fuc} {resname AFUC and same residue as (name C1 and within 2 of (resname BGLN and name O6))} {BGLN_16_AFUC.dx} 1 {0.01 0.05} 0 {O6 C6 C5 O5}}
      {{GlcNAc (beta-1,4) Man} {resname BMAN and same residue as (name C1 and within 2 of (resname BGLN and name O4))} {BGLN_14_BMAN.dx} 0 {0.01 0.05} 0 {O4 C4 C3}}
      {{Man (alpha-1,2) Man} {resname AMAN and same residue as (name C1 and within 2 of (resname AMAN and name O2))} {AMAN_12_AMAN.dx} 0 {0.01 0.05} 0 {O2 C2 C1}}
      {{Man (alpha-1,2) GlcNAc} {resname BGLN and same residue as (name C1 and within 2 of (resname AMAN and name O2))} {AMAN_12_BGLN.dx} 0 {0.01 0.05} 0 {O2 C2 C1}}
      {{Man (alpha-1,4) GlcNAc} {resname BGLN and same residue as (name C1 and within 2 of (resname AMAN and name O4))} {AMAN_14_BGLN.dx} 0 {0.01 0.05} 0 {O4 C4 C3}}
      
  }
  
  variable torrepids ""                         ;# list of rep indices for the different torsion cases
  variable molsellist ""                    ;# list of atomselections in the molecule you're analysing
  variable torsellist ""                    ;# list of atomselections in torplot molecule corresponding to $molsellist
  variable torselnamelist ""                ;# list of names associated with each atomselection in $torsellist
  variable torplotdir $env(TORPLOTDIR)
  variable visualisation ""
  variable 2Daxis "-1"
  variable 3Daxis "-1"
  variable 2D3D ""
  
}

proc torsionplot {} {
  ::TorPlot::torsionplot
}

proc ::TorPlot::torsionplot {} {
  # Just create the window and initialize data structures
  # No molecule has been selected yet
  # Also set up traces on VMD variables
  
  variable torsel
  variable w
    
  global vmd_frame
  global vmd_initialize_structure

  # If already initialized, just turn on 
  if [winfo exists .rama] {
    wm deiconify $w
    return
  }

  set w [toplevel ".rama"]
  wm title $w "TorsionPlot - Ramachandran plots and similar metrics for biomolecules"
  wm resizable $w 0 0
  bind $w <Destroy> ::TorPlot::destroy
   
  frame $w.top

  # Create menubar
  frame $w.top.menubar -relief raised -bd 2
  pack $w.top.menubar -padx 1 -fill x -side top
  menubutton $w.top.menubar.file -text "File   " -underline 0 -menu $w.top.menubar.file.menu
  $w.top.menubar.file config -width 5
  pack $w.top.menubar.file -side left

  # File menu
  menu $w.top.menubar.file.menu -tearoff no
  $w.top.menubar.file.menu add command -label "Save snapshot" \
        -command [namespace code takeSnapshot]


  # Help menu
  menubutton $w.top.menubar.help -text "Help   " -menu $w.top.menubar.help.menu
  $w.top.menubar.help config -width 5
  pack $w.top.menubar.help -side right 
  menu $w.top.menubar.help.menu -tearoff no
  $w.top.menubar.help.menu add command -label "Torsionplot Help..." -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/ramaplot"


## 
 #   ramaCanvas $w.fr -width 370 -height 370
 #   pack $w.fr -in $w.top -fill both -expand true -side left
 # 
 #   ramaZones $w.fr.canvas
 ##

  frame $w.data

  # Create molecule selection menu
  frame $w.data.mol
  label $w.data.mol.l -text Molecule -anchor w
  menubutton $w.data.mol.m -relief raised -bd 2 -direction flush \
	-textvariable ::TorPlot::molid \
	-menu $w.data.mol.m.menu
  menu $w.data.mol.m.menu 
   # trace variable ::TorPlot::molid w [namespace code torplotChangeMolecule]

  pack $w.data.mol.l -side left
  pack $w.data.mol.m -side right
  pack $w.data.mol -side top

  # Create atom selection entry
  frame $w.data.sel
  label $w.data.sel.l -text Selection -anchor w
  entry $w.data.sel.e -relief sunken -width 18 -bg White \
	-textvariable [namespace current]::seltext
  # bind $w.data.sel.e <Return> [namespace code torplotChangeMolecule]
  pack $w.data.sel.l -side left
  pack $w.data.sel.e -side right
  pack $w.data.sel -side top
  
  # Create menu for selection of torsions to visualise
  frame $w.data.torvis
  label $w.data.torvis.l -text "Torsion to inspect" -anchor w
  menubutton $w.data.torvis.m -relief raised -bd 2 -direction flush \
    -textvariable ::TorPlot::visualisation \
    -menu $w.data.torvis.m.menu
  menu $w.data.torvis.m.menu
  
  pack $w.data.torvis.l -side left
  pack $w.data.torvis.m -side right
  pack $w.data.torvis -side top
  trace variable ::TorPlot::visualisation w [namespace code torplotChangeVis]
  
  # Buttons to initialise and update plots
  frame $w.data.plot
  button $w.data.plot.initialize -text "Initialise plot for current selection" -command [namespace current]::torplotInitialize
  button $w.data.plot.update -text "plot torsions for current frame" -command [namespace current]::torplotUpdate
  button $w.data.plot.reset -text "Reset TorsionPlot and remove all visualisations" -command [namespace current]::torplotReset
  pack $w.data.plot.initialize -side left
  pack $w.data.plot.update -side top
  pack $w.data.plot.reset -side top
  pack $w.data.plot -side top
  
  
  pack $w.data -in $w.top -side right
  pack $w.top

  ## 
   # # Create fields for displaying info about last clicked residue
   # frame $w.data.info
   # foreach field {Segid Resname Resid Phi Psi} {
   #   label $w.data.info.l$field -text $field -anchor w
   #   entry $w.data.info.e$field -relief sunken -width 10
   #   grid $w.data.info.l$field $w.data.info.e$field -sticky news
   # }
   # pack $w.data.info -side top
   ##

  # 3d rama plot buttons.
## 
 #   frame  $w.data.space -height 50
 #   button $w.data.mk3d  -text {Create 3-d Histogram} -command [namespace code mk3drama]
 #   button $w.data.del3d -text {Delete 3-d Histogram} -command [namespace code del3drama]
 #   pack $w.data.space $w.data.mk3d $w.data.del3d -side top
 # 
 #   pack $w.data -in $w.top -side right
 # 
 #   pack $w.top
 #   # Draw grid lines on the x and y axes
 #   $w.fr.canvas create line 5 185 365 185 -tag grid
 #   $w.fr.canvas create line 185 5 185 365 -tag grid
 #   $w.fr.canvas bind residue <Button-1> [namespace code { ramaHighlight %x %y}]
 #   $w.fr.canvas bind line <Button-1> [namespace code { ramaGoto %W}]
 # 
 #   # Update the marks every time there's a new frame
 #   trace variable vmd_frame w [namespace code ramaUpdate]
 ##

  # Update the molecules when molecules are deleted or added
  trace variable vmd_initialize_structure w [namespace code torplotUpdateMolecules]
  trace variable ::TorPlot::torselnamelist w [namespace code torplotUpdateVisList]

  # Set up the molecule list
  ::TorPlot::torplotUpdateMolecules
}

## 
 # # Finds the 
 # proc ::TorPlot::ramaGoto { w } {
 #   set id [ $w find withtag current]
 #   set taglist [$w gettags $id]
 #   set listindex [lsearch -glob $taglist frame*]
 #   if { $listindex < 0 } {
 #     return
 #   }
 #   set frametag [lindex $taglist $listindex]
 #   lassign [split $frametag :] foo frame
 #   animate goto $frame
 # }
 ##


proc ::TorPlot::takeSnapshot { args } {
  variable w
  variable lastsavedfile

  set filename [tk_getSaveFile \
        -initialfile $lastsavedfile \
        -title "take snapshot of current display" \
        -parent $w \
        -filetypes [list {{TARGA files} {.tga}} {{All files} {*}}]]
  if {$filename != ""} {
    render snapshot $filename
    set lastsavedfile $filename
  }
  return
}

proc ::TorPlot::destroy { args } {
  # Delete traces
  # Delete remaining selections

  variable torsel
  variable molid
  global vmd_frame
  global vmd_initialize_structure
  
  trace vdelete ::TorPlot::torselnamelist w [namespace code torplotUpdateVisList]
  trace vdelete vmd_initialize_structure w [namespace code torplotUpdateMolecules]
  trace vdelete ::TorPlot::visualisation w [namespace code torplotChangeVis]

  ## 
   # trace vdelete molid w [namespace code torplotChangeMolecule]
   # trace vdelete vmd_frame w [namespace code torplotUpdate]
   # trace vdelete vmd_initialize_structure w [namespace code torplotUpdateMolecules]
   ##

  catch {$torsel delete}
}

proc ::TorPlot::torplotUpdateMolecules { args } {
  variable torsel
  variable w
  variable molid

  set mollist [molinfo list]
  
  # Invalidate the selection if necessary
  if { [lsearch $mollist $molid] < 0 } {
    catch {$torsel delete}
    set torsel ""
  }    

  # Update the molecule browser
  $w.data.mol.m.menu delete 0 end
  $w.data.mol.m configure -state disabled
  if { [llength $mollist] != 0 } {
    foreach id $mollist {
      if {[molinfo $id get filetype] != "graphics"} {
        $w.data.mol.m configure -state normal 
        $w.data.mol.m.menu add radiobutton -value $id \
	  -label "$id [molinfo $id get name]" \
	  -variable ::TorPlot::molid 
      }
    }
  }
}


proc ::TorPlot::torplotUpdateVisList { args } {
    variable w
    variable torselnamelist
    variable visualisation
    
    $w.data.torvis.m.menu delete 0 end
    $w.data.torvis.m configure -state disabled
    if { [llength $torselnamelist] !=0 } {
	for {set i 0} {$i < [llength $torselnamelist]} {incr i} {
	    $w.data.torvis.m.menu add radiobutton -value $i \
	      -label "[lindex $torselnamelist $i]" \
	      -variable ::TorPlot::visualisation
	}
	$w.data.torvis.m configure -state normal
    }
}



proc ::TorPlot::torplotChangeVis { args } {
    variable w
    variable torrepids
    variable visualisation
    variable torplotmol
    variable 2D3D
    variable 2Daxis
    variable 3Daxis
	variable tortypes
	variable torsellist
    
    # hide all visualisations
    
    if {[llength $torrepids] > 0} {
	    foreach repset $torrepids {
	        foreach rep $repset {
		        mol showrep $torplotmol $rep off
	        }
	    }
	    mol off $2Daxis
	    mol off $3Daxis
	
	    # turn on visualisations for selected rep
	
	    foreach rep [lindex $torrepids $visualisation] {
	        mol showrep $torplotmol $rep on
	    }
	    set is3D [lindex $2D3D $visualisation]
	    if {$is3D} {
	        mol on $3Daxis
	    } else {
			mol on $2Daxis
		}
		
	    #list outliers (best put in its own proc later)
		set cutoff [lindex [lindex [lindex $tortypes $visualisation] 4] 0]
        catch {set outliers [atomselect $torplotmol "index [[lindex $torsellist $visualisation] get index] and interpvol$visualisation < $cutoff"]
		  puts "Outliers at $cutoff cutoff: [$outliers get "segname resid interpvol$visualisation"]"
		  $outliers delete
        }
    }
}

proc ::TorPlot::torplotReset { args } {
    variable torplotmol
    variable molsellist
    variable torsellist
    variable torrepids
    variable torselnamelist
    variable phiaxis
    variable psiaxis
    variable omegaaxis
    variable 2Daxis
    variable 3Daxis
    variable 2D3D
  

    mol delete $torplotmol
    mol delete $2Daxis
    mol delete $3Daxis
    set torplotmol ""
    foreach sel $molsellist {
       catch {$sel delete}
    }
    set molsellist ""
    foreach sel $torsellist {
       catch {$sel delete}
    }
    set torsellist ""
    set torrepids ""
    set torselnamelist ""
    set 2D3D ""

}


	    

proc ::TorPlot::torplotInitialize { args } {
  variable w
  variable torsel
  variable seltext
  variable molid
  variable tempdir
  variable torfilename
  variable torplotmol
  variable torplotall
  variable tortypes
  variable torrepids
  variable molsellist
  variable torsellist
  variable torselnamelist
  variable torplotdir
  variable 2Daxis
  variable 3Daxis
  variable 2D3D


  if { $molid == "" || [lsearch [molinfo list] $molid] < 0} {
    return
  }
   
#  wm title $w "Ramachandran plot for molecule $molid [molinfo $molid get name]"
  if { $seltext == "" } {
    set seltext all
  }
  if {![catch {set sel [atomselect $molid "((protein and name CA and not (same residue as name HT1 HN1)) or (glycan and name C1)) and $seltext"]}]} {
    catch {$torsel delete}
    set torsel $sel
    $torsel global
    $sel delete
  } else {
    puts "Unable to create new selection!"
    return
  }
  
  # Write C-alpha atoms of protein and C1 atoms of glycans to a new PDB file. These will become the points in the Phi-Psi plot
  $torsel writepdb $torfilename


  # set up new molecules for visualisations. The probability distribution for each case is saved in a 3D or pseudo-2D .dx volumetric map
  display update off
  
  # load the PDB file we just saved. This will be used to write the Phi/Psi coordinates 
  mol new $torfilename
  set torplotmol [molinfo top]
  mol delrep 0 $torplotmol
  set allreps ""
  set repcount 0
  # Set up representations for each case, and atomselections for easy handling later. Repids are stored in $torrepids in blocks of 4
  
  for {set i 0} {$i < [llength $tortypes]} {incr i} {
      # set up atomic representations
      set torseltext "name CA C1 and (none"
      set thisseltext [lindex [lindex $tortypes $i] 1]
      set sel [atomselect $molid "name CA C1 and not (protein and same residue as name HT1 HN1) and ($thisseltext) and ($seltext)"]
      foreach seg [lsort -unique [$sel get segname]] {
		set sel2 [atomselect $molid "segname $seg and index [$sel get index]"]
	  	if {[$sel2 num] > 0} {
	      append torseltext " or (segname $seg and resid [$sel2 get resid])"
	  	}
	  	$sel2 delete
      }
      append torseltext ")"

      # create atomselections for this case in both the molecule being analysed and the torsion plot atoms
      lappend molsellist $sel
      $sel global
      
      set thistorsel [atomselect $torplotmol $torseltext]
      lappend torsellist $thistorsel
      $thistorsel global

      set thisreplist ""
      set is3D [lindex [lindex $tortypes $i] 3]

      if {$is3D} {
	  mol representation VDW 4.0 8.0
      } else {
	  mol representation VDW 2.0 8.0
      }
      
  
      mol color ColorID 4
      mol selection $torseltext
      mol material Opaque
      mol addrep $torplotmol
      lappend thisreplist $repcount
      mol showrep $torplotmol $repcount off
      incr repcount
      
      # load associated map and set up representations in 2 or 3 dimensions
      mol addfile [file join $torplotdir [lindex [lindex $tortypes $i] 2]]
      
      
      if {!$is3D} {
	  mol representation VolumeSlice 0 $i 2 2
	  mol color Name
	  mol selection all
	  mol material Opaque
	  mol addrep $torplotmol
	  lappend thisreplist $repcount
	  mol showrep $torplotmol $repcount off
	  incr repcount
      }
      
      set cutoffs [lindex [lindex $tortypes $i] 4]
      set colorlist {1 4 7 0}
      set c 0
      foreach cutoff $cutoffs {
	  if {$is3D} {
	      	  mol representation Isosurface $cutoff $i 2 1 2 1
	      } else {
	  mol representation Isosurface $cutoff $i 2 1 1 4
      }
      
	      
      mol color colorID [lindex $colorlist $c]
      mol selection all
      mol material Opaque
      mol addrep $torplotmol
      lappend thisreplist $repcount
      mol showrep $torplotmol $repcount off
      incr repcount
      incr c
  }
      	  
      lappend torrepids $thisreplist
      lappend torselnamelist [lindex [lindex $tortypes $i] 0]
      lappend 2D3D $is3D
      
     
  }
  
      # Process data and boxes
  # data has unique residues as keys and {{segid resid resname} id phi psi}
  # as values.
  # Populate data array with names.  Since the names presumably won't change
  # we can do this now and save some time in ramaHighlight.  We store the id
  # in the data array as well so that when we update, we can easily move the
  # correct box.  We use the box array only to look up selections.


  ## 
   # foreach residue [$selection get residue] \
   #         namedata [$selection get {segid resid resname}] {
   #   set x1 [expr -$xmin - $recsize]
   #   set y1 [expr $ymax + $recsize]
   #   set x2 [expr $x1 + 2 * $recsize]
   #   set y2 [expr $y1 - 2 * $recsize]
   #   set id [$w.fr.canvas create rectangle $x1 $y1 $x2 $y2 \
   #               -fill yellow -tags residue]
   #   set box($id) $residue
   #   set data($residue) [list $namedata $id 0 0]
   # }
   ##
    # fill the display with the plot
  molinfo $torplotmol set center_matrix {{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}
  molinfo $torplotmol set scale_matrix {{{0.008 0 0 0} {0 0.008 0 0} {0 0 0.008 0} {0 0 0 1}}}
  molinfo $torplotmol set rotate_matrix {{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}
  molinfo $torplotmol set global_matrix {{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}
  
  set 2Daxis [mol new]
  mol rename $2Daxis "2D axes"
  
  set 3Daxis [mol new]
  mol rename $3Daxis "3D axes"
  
  ::TorPlot::draw_axes
  
 
  TorPlot::torplotUpdate
  
}

proc ::TorPlot::draw_axes { args } {
    variable 2Daxis
    variable 3Daxis
    
	graphics $2Daxis color red
	graphics $2Daxis cylinder {-185 -185 0} {185 -185 0} radius 2 resolution 8
	for {set i -180} {$i <= 180} {incr i 30} {
	    graphics $2Daxis cylinder "$i -185 0" "$i -180 0" radius 0.5 resolution 8
	    graphics $2Daxis text "$i -200 0" "$i" size 0.7 thickness 1
	}
	graphics $2Daxis text {-10 -220 0} "Phi" size 2 thickness 2

	graphics $2Daxis cylinder {-185 -185 0} {-185 185 0} radius 2 resolution 8
	for {set i -180} {$i <= 180} {incr i 30} {
	    graphics $2Daxis cylinder "-185 $i 0" "-180 $i 0" radius 0.5 resolution 8
	    graphics $2Daxis text "-220 $i 0" "$i" size 0.7 thickness 1
	}
	graphics $2Daxis text {-260 0 0} "Psi" size 2 thickness 2



	graphics $3Daxis color red
	graphics $3Daxis cylinder {-185 -185 -185} {185 -185 -185} radius 2 resolution 8
	for {set i -180} {$i <= 180} {incr i 30} {
	    graphics $3Daxis cylinder "$i -185 -185" "$i -180 -180" radius 0.5 resolution 8
	    graphics $3Daxis text "$i -200 -200" "$i" size 0.7 thickness 1
	}
	graphics $3Daxis text {-10 -220 -210} "Phi" size 2 thickness 2

	graphics $3Daxis cylinder {-185 -185 -185} {-185 185 -185} radius 2 resolution 8
	for {set i -180} {$i <= 180} {incr i 30} {
	    graphics $3Daxis cylinder "-185 $i -185" "-180 $i -180" radius 0.5 resolution 8
	    graphics $3Daxis text "-220 $i -200" "$i" size 0.7 thickness 1
	}
	graphics $3Daxis text {-260 0 -210} "Psi" size 2 thickness 2
	
    graphics $3Daxis cylinder {-185 -185 -185} {-185 -185 185} radius 2 resolution 8
    for {set i -180} {$i <= 180} {incr i 30} {
	graphics $3Daxis cylinder "-185 -185 $i" "-180 -180 $i" radius 0.5 resolution 8
	graphics $3Daxis text "-200 -200 $i" "$i" size 0.7 thickness 1
    }
    graphics $3Daxis text {-250 -220 0} "Omega" size 2 thickness 2
    
    mol off $3Daxis
    mol off $2Daxis
    
}

   
   
proc ::TorPlot::torplotUpdate { args } {
  variable w
  variable torsel
  variable seltext
  variable molid
  variable tempdir
  variable torfilename
  variable torplotmol
  variable torplotall
  variable tortypes
  variable torrepids
  variable molsellist
  variable torsellist
  
  display update off
  
  if { ![string compare $torsel ""] } {
    return
  }

  for {set i 0} {$i < [llength $tortypes]} {incr i} {
      set isProtein [lindex [lindex $tortypes $i] 5]
      set is3D [lindex [lindex $tortypes $i] 3]
      set thismolsellist [lindex $molsellist $i]
      if {$thismolsellist == ""} {
         continue
      }
      set thismolsel [lindex $molsellist $i]
      set thistorsel [lindex $torsellist $i]

      if {$isProtein} {	  
	  
	  # then everything is easy since phi and psi are predefined
	  $thistorsel set {x y} [$thismolsel get {phi psi}]
	  $thistorsel set z 0
      } else {     
	  
	  # things get a little more complex
	  
	  set coorlist {}
      
	  foreach res [$thismolsel get residue] {
	      set bondedatomlist [lindex [lindex $tortypes $i] 6]
	      set bondedres [atomselect $molid "not (residue $res) and same residue as within 2 of (residue $res and name C1)"]
	      set bondedresidue [$bondedres get residue]
	      $bondedres delete
	      set indexlist ""
	      set thisO5 [atomselect $molid "residue $res and name O5"]
	      lappend indexlist [$thisO5 get index]
	      $thisO5 delete
	      set thisC1 [atomselect $molid "residue $res and name C1"]
	      lappend indexlist [$thisC1 get index]
	      $thisC1 delete
	      foreach atom $bondedatomlist {
		  set thisatom [atomselect $molid "residue $bondedresidue and name $atom"]
		  lappend indexlist [$thisatom get index]
		  $thisatom delete
	      }
	      set phi [measure dihed [lrange $indexlist 0 3] molid $molid]
	      set psi [measure dihed [lrange $indexlist 1 4] molid $molid]
	      if {$is3D} {
		  set omega [measure dihed [lrange $indexlist 2 5] molid $molid]
		  lappend coorlist "$phi $psi $omega"
	      } else {
		  lappend coorlist "$phi $psi 0"
	      }
	  }
	  
	  $thistorsel set {x y z} $coorlist
      }
  }
  
	      
      

  display update on
}


proc torsionplot_tk {} {
  ::TorPlot::torsionplot
  return $::TorPlot::w
}

