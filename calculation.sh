#!/bin/bash
# This script was developed to calculate:
# 1. Geometry
# 2. Band structure
# 3. Density of states

# It's a basic but useful script. It was developed based on the input and output of the Quantum Espresso 6.8V.
# Developed by Igor Peixoto. Email: igorpeixoto_@outlook.com.br
#____________________________________________________________________________________________________
# 						GLOBAL VARIABLES
#____________________________________________________________________________________________________
# File name
input=si_n_sitio_a_

# Pseudopotential used, if there is more than one pseudopotential, simply create the corresponding variable and add it to the other parts of the script
pseudopotencial1="Si   28.08500  Si.pbe-nl-rrkjus_psl.1.0.0.UPF" 
pseudopotencial2="N    14.00674  N.pbe-n-rrkjus_psl.1.0.0.UPF"

#Types of atoms
atomo1=Si
atomo2=N
atomo3=Li

# Number of different atoms
tipo_atomos=3

# Wave function energy value
ewfc=80
# Energy value of energy density. For pseudopotentials of the type: #PAW= 4*ewfc and USPP= 8*ewfc
erho=640  

# Number of points ks
k_points="K_POINTS {automatic}
 12 12 12   0 0 0
"
# Simple cubic network 
# Points ks in the Brillouin zone, modify for each type of bravais cell"
k_points_bands="K_POINTS {crystal_b}
6
0.000000000   0.000000000   0.000000000 30 !GAMMA
0.500000000   0.000000000   0.000000000 30 !X
0.500000000   0.500000000   0.000000000 30 !M
0.000000000   0.000000000   0.000000000 30 !GAMMA
0.500000000   0.500000000   0.500000000 30 !R
0.500000000   0.000000000   0.000000000 30 !X
"
# Number of bands to be calculated
num_bands=32

# type of electron occupation, for DoS calculating 
occupations_dos="'tetrahedra'"

#Mixing electronic densities
mixbeta=0.7d0

#_______________________________________________________________________________________________________
# 						GEOMETRY CALCULATION
#_______________________________________________________________________________________________________

# Starting the test to check if the geometry calculation has already been performed
# The code is searching for the contents of the file with the specified name ($input$tipo.out), looking for the word "JOB" in that content and, if found, cutting a specific part of that line. The final result is stored in the test variable.

tipo="vc-relax"
teste=`grep JOB $input$tipo.out | cut -b4-6`
     check="JOB"
        if [[ $teste == $check ]];
           then
               echo "Geometry calculation already carried out"
	       else	
               tipo="vc-relax"

# This part od the code is creating or overwriting a file with the specified name ($input$tipo.in) and allowing you to input the contents of the file.
# Remember to replace all the information about your system

cat > $input$tipo.in <<EOF
&CONTROL
        title          = '$input$tipo' ,
        calculation    = '$tipo' ,
        outdir         ='temp_file',
        prefix         ='$input',
        pseudo_dir     ='/home/iprodrigues/pseudopotenciais',
        verbosity      ='low',
        restart_mode   ='from_scratch',
        tprnfor        =.true.,
        tstress        =.true.,
        forc_conv_thr  =1.0d-6
        etot_conv_thr  =1.0d-8
        nstep          =1000
 /

 &SYSTEM
        ibrav       = 0,
	A           = 5.43094100,
        nat         = 8,
        ntyp        = $tipo_atomos,
        ecutwfc     = $ewfc ,
        ecutrho     = $erho ,
        occupations = 'smearing',
        smearing    = 'mv',
        degauss     = 0.005d0,
        lspinorb    = .false.,
        noncolin    = .false.,
        nbnd        = $num_bands,
 /

&ELECTRONS
        electron_maxstep = 1000,
        conv_thr         = 1.0d-12,
        mixing_beta      = $mixbeta,
        mixing_mode      = 'plain',
/

&IONS
        ion_dynamics = 'bfgs'
        upscale      = 100.d0
/

&CELL
        cell_dofree    = "all"
        cell_dynamics  = "bfgs"
        cell_factor    =  2.0d0
        press_conv_thr = 0.005d0
        press          = 0.0d0
/

ATOMIC_SPECIES
$pseudopotencial1
$pseudopotencial2

ATOMIC_POSITIONS {crystal}
N             0.0000000000        0.0000000000        0.0000000000  
Si            0.0000000000        0.5000000000        0.5000000000
Si            0.5000000000        0.0000000000        0.5000000000
Si            0.5000000000        0.5000000000        0.0000000000
Si            0.7500000000        0.7500000000        0.2500000000
Si            0.7500000000        0.2500000000        0.7500000000
Si            0.2500000000        0.7500000000        0.7500000000
Si            0.2500000000        0.2500000000        0.2500000000
$k_points

CELL_PARAMETERS {crystal}
1.006765926   0.000000000   0.000000000
0.000000000   1.006765926   0.000000000
0.000000000   0.000000000   1.006765926
EOF

set OMP_NUM_THREADS=1; export OMP_NUM_THREADS=1
/usr/bin/mpirun.mpich -np 16 /home/usr/q-e-qe-6.7.0/bin/pw.x -i  $input$tipo.in > $input$tipo.out 

fi

#_______________________________________________________________________________________________________
# 						SCF CALCULATION
#_______________________________________________________________________________________________________	
# Starting the test to check if the scf calculation has already been performed

tipo="scf"
teste=`grep JOB $input$tipo.out | cut -b4-6`
     check="JOB"
     if [[ $teste == $check ]];
         then
             echo "SCF calculation already performed"
		    tipo="vc-relax"
                    # https://www.theunixschool.com/2012/05/different-ways-to-print-next-few-lines.html?m=1
                    # Optimized network parameters, obtained from the vc-relax calculation
                    # This line uses the awk command to find lines in the file $input$tipo.out that match the pattern /Begin final/ (that is, lines that contain the string "Begin final").
                    # When these lines are found, awk prints the current line number (NR) increased by 5. The result (line number) is stored in the variable x
		    x=`awk '/Begin final/{print NR+5}' $input$tipo.out`
                    y=`awk '/Begin final/{print NR+7}' $input$tipo.out`
                    
                    #The sed command is used to extract a range of lines from the file $input$tipo.out, starting at line x and ending at line y. 
                    # The result is stored in the param_rede variable.
                    param_rede=$(sed -n "$x,$y p" $input$tipo.out)
         
                     # Optimized coordinates, obtained from the vc-relax calculation
		     xx=`awk '/Begin final/{print NR+9}' $input$tipo.out`
                     yy=`awk '/Begin final/{print NR+17}' $input$tipo.out`
                     geometria=$(sed -n "$xx,$yy p" $input$tipo.out)
	  else	
	     echo "SCF calculation will be performed"
             tipo="vc-relax"

             #CALCULO SCF
             calc1=1
             while [ $calc1 -eq 1 ]
                do
                    teste=`grep JOB $input$tipo.out | cut -b4-6`
                    check="JOB"
                 if [ $teste == $check ]
                    then
                    calc1=0
		
	            # https://www.theunixschool.com/2012/05/different-ways-to-print-next-few-lines.html?m=1
                    # Optimized network parameters, obtained from the vc-relax calculation
	            x=`awk '/Begin final/{print NR+5}' $input$tipo.out`
                    y=`awk '/Begin final/{print NR+7}' $input$tipo.out`
                    param_rede=$(sed -n "$x,$y p" $input$tipo.out)
         
               	    # Optimized coordinates, obtained from the vc-relax calculation
		    xx=`awk '/Begin final/{print NR+9}' $input$tipo.out`
                    yy=`awk '/Begin final/{print NR+17}' $input$tipo.out`
                    geometria=$(sed -n "$xx,$yy p" $input$tipo.out)
		    tipo="scf"

cat > $input$tipo.in <<EOF
&CONTROL
        title          = '$input$tipo' ,
        calculation    = '$tipo' ,
        outdir         ='temp_file',
        prefix         ='$input',
        pseudo_dir     ='/home/iprodrigues/pseudopotenciais',
        verbosity      ='high',
        restart_mode   ='from_scratch',
        tprnfor        =.true.,
        tstress        =.true.,
        forc_conv_thr  =1.0d-6
        etot_conv_thr  =1.0d-8
        nstep          =1000
 /

 &SYSTEM
        ibrav       = 0,
	A           = 5.43094100,
        nat         = 8,
        ntyp        = $tipo_atomos,
        ecutwfc     = $ewfc ,
        ecutrho     = $erho ,
        occupations = 'smearing',
        smearing    = 'mv',
        degauss     = 0.005d0,
        lspinorb    = .false.,
        noncolin    = .false.,
        nbnd        = $num_bands,
 /

&ELECTRONS
        electron_maxstep = 1000,
        conv_thr         = 1.0d-12,
        mixing_beta      = $mixbeta,
        mixing_mode      = 'plain',
/

ATOMIC_SPECIES
$pseudopotencial1
$pseudopotencial2

$geometria

$k_points

CELL_PARAMETERS {crystal}
$param_rede

EOF
                     set OMP_NUM_THREADS=1; export OMP_NUM_THREADS=1
                     /usr/bin/mpirun.mpich -np 16 /home/usr/q-e-qe-6.7.0/bin/pw.x -i  $input$tipo.in > $input$tipo.out 
                 else
                     calc1=1

                 fi
                done 
    fi

#_______________________________________________________________________________________________________
# 						BAND STRUCTURE CALCULATION
#_______________________________________________________________________________________________________

# Starting the test to check if the band calculation has already been performed
tipo="bands"
teste=`grep JOB $input$tipo.out | cut -b4-6`
     check="JOB"
    if [[ $teste == $check ]];
         then
             echo "Band calculation already carried out"
	    else	
		     echo "Band calculation will be carried out" 
            tipo="scf"
            #CALCULO DE BANDAS
            calc2=1
            while [ $calc2 -eq 1 ]
               do
                 teste=`grep JOB $input$tipo.out | cut -b4-6`
                 check="JOB"
                 if [ $teste == $check ]
                   then
                     calc2=0
                     tipo=bands
cat > $input$tipo.in <<EOF
&CONTROL
        title          = '$input$tipo' ,
        calculation    = '$tipo' ,
        outdir         ='temp_file',
        prefix         ='$input',
        pseudo_dir     ='/home/iprodrigues/pseudopotenciais',
        verbosity      ='high',
        restart_mode   ='from_scratch',
        tprnfor        =.true.,
        tstress        =.true.,
        forc_conv_thr  =1.0d-6
        etot_conv_thr  =1.0d-8
        nstep          =1000
 /

 &SYSTEM
        ibrav       = 0,
        nat         = 8,
	A           = 5.43094100,		
        ntyp        = $tipo_atomos,
        ecutwfc     = $ewfc ,
        ecutrho     = $erho ,
        occupations = 'smearing',
        smearing    = 'mv',
        degauss     = 0.005d0,
        lspinorb    = .false.,
        noncolin    = .false.,
	nbnd= $num_bands ,
 /

&ELECTRONS
        electron_maxstep = 1000,
        conv_thr         = 1.0d-12,
        mixing_beta      = $mixbeta,
        mixing_mode      = 'plain',
/

ATOMIC_SPECIES
$pseudopotencial1
$pseudopotencial2

$geometria

$k_points_bands

CELL_PARAMETERS {crystal}
$param_rede

EOF
                     set OMP_NUM_THREADS=1; export OMP_NUM_THREADS=1
                     /usr/bin/mpirun.mpich -np 16 /home/usr/q-e-qe-6.7.0/bin/pw.x -i  $input$tipo.in > $input$tipo.out 
                 else
                    calc2=1
                 fi
               done 
    fi

#_______________________________________________________________________________________________________
# 						BAND STRUCTURE AUXILIAR CALCULATION
#_______________________________________________________________________________________________________
#Starting the test to check if the auxiliary band calculation has already been carried out
tipo="bands"
aux="_aux"
teste=`grep JOB $input$tipo$aux.out | cut -b4-6`
     check="JOB"
     if [[ $teste == $check ]];
         then
             echo "Auxiliary band calculation already carried out"
	    else	
	     echo "Auxiliary band calculation will be carried out"

            #CALCULO AUXILAR DE BANDAS
             calc3=1
            while [ $calc3 -eq 1 ]
              do
               teste=`grep JOB $input$tipo.out | cut -b4-6`
               check="JOB"
              if [[ $teste == $check ]];
               then
               tipo="bands"
               aux=_aux
cat > $input$tipo$aux.in <<EOF
&BANDS
       outdir = 'temp_file',
       prefix = '$input',
       filband= '$input$tipo.dat',
/
EOF
               set OMP_NUM_THREADS=1; export OMP_NUM_THREADS=1
               /usr/bin/mpirun.mpich -np 16 /home/usr/q-e-qe-6.7.0/bin/bands.x -i  $input$tipo$aux.in > $input$tipo$aux.out
               calc3=0
               else
               calc3=1
              fi
              done
     fi

#_______________________________________________________________________________________________________
# 						NSCF CALCULATION
#_______________________________________________________________________________________________________

#Starting the test to check if the nscf calculation has already been performed
tipo="nscf"
teste=`grep JOB $input$tipo.out | cut -b4-6`
     check="JOB"
     if [[ $teste == $check ]];
         then
             echo "Calculo nscf já realizado"
	    else	
		     echo "Calculo nscf será realizado"
             tipo=nscf
cat > $input$tipo.in <<EOF
&CONTROL
        title          = '$input$tipo' ,
        calculation    = '$tipo' ,
        outdir         ='temp_file',
        prefix         ='$input',
        pseudo_dir     ='/home/iprodrigues/pseudopotenciais',
        verbosity      ='high',
        restart_mode   ='from_scratch',
        tprnfor        =.true.,
        tstress        =.true.,
        forc_conv_thr  =1.0d-6
        etot_conv_thr  =1.0d-8
        nstep          =1000
 /

 &SYSTEM
        ibrav       = 0,
        A           = 5.43094100,
        nat         = 8,
        ntyp        = $tipo_atomos,
        ecutwfc     = $ewfc ,
        ecutrho     = $erho ,
        occupations = $occupations_dos,
        smearing    = 'mv',
        degauss     = 0.005d0,
        lspinorb    = .false.,
        noncolin    = .false.,
        nbnd        = $num_bands,
 /

&ELECTRONS
        electron_maxstep = 1000,
        conv_thr         = 1.0d-12,
        mixing_beta      = $mixbeta,
        mixing_mode      = 'plain',
/

ATOMIC_SPECIES
$pseudopotencial1
$pseudopotencial2

$geometria

$k_points

CELL_PARAMETERS {crystal}
$param_rede

EOF
             set OMP_NUM_THREADS=1; export OMP_NUM_THREADS=1
             /usr/bin/mpirun.mpich -np 16 /home/usr/q-e-qe-6.7.0/bin/pw.x -i  $input$tipo.in > $input$tipo.out 

     fi 

#_______________________________________________________________________________________________________
# 						DoS AUXILIAR CALCULATION
#_______________________________________________________________________________________________________

#CALCULO AUXILAR DOS
# Starting the test to check if the auxiliary DoS calculation has already been performed
tipo="dos"
aux="_aux"
teste=`grep JOB $input$tipo$aux.out | cut -b4-6`
     check="JOB"
     if [[ $teste == $check ]];
         then
             echo "Auxiliary calculation of dos already carried out"
	      else	
		     echo "Auxiliary calculation of days will be performedo"
s
            tipo="nscf"
            calc5=1
            while [ $calc5 -eq 1 ] #enquanto calc for igual a 1 faça
             do
              teste=`grep JOB $input$tipo.out | cut -b4-6` #procurando pela palavra chave do calculo alterior
              check="JOB"
              if [ $teste == $check ] # se a palavra teste for igual a check entao
              then
                 tipo=dos
                 aux=_aux
cat > $input$tipo$aux.in <<EOF
&DOS
       outdir ='temp_file',
       prefix ='$input',
       Emin   = -40, ! valor de menor energia encontrado nas bandas
       Emax   = 40,  !Valor de maior energia encontrado nas bandas
       DeltaE =0.05,
       fildos = '$input$tipo.dat',

/
EOF
                  set OMP_NUM_THREADS=1; export OMP_NUM_THREADS=1
                  /usr/bin/mpirun.mpich -np 16 /home/usr/q-e-qe-6.7.0/bin/dos.x -i  $input$tipo$aux.in > $input$tipo$aux.out
                  calc5=0
               else
                  calc5=1
              fi
             done
     fi



#CALCULO PDOS
# Iniciando o teste para verificar se o calculo PDoS já foi realizado
tipo="pdos"
teste=`grep JOB $input$tipo.out | cut -b4-6`
     check="JOB"
     if [[ $teste == $check ]];
         then
             echo "Calculo auxiliar de pdos já realizado"
              else
                     echo "Calculo auxiliar de pdos será realizado"

            tipo="nscf"
            calc6=1
            while [ $calc6 -eq 1 ] #enquanto calc for igual a 1 faça
             do
              teste=`grep JOB $input$tipo.out | cut -b4-6` #procurando pela palavra chave do calculo al$
              check="JOB"
              if [ $teste == $check ] # se a palavra teste for igual a check entao
              then
                 tipo=pdos

cat > $input$tipo.in <<EOF
&PROJWFC
       outdir ='temp_file',
       prefix ='$input',
       lwrite_overlaps =.true.,
       ngauss=-1,
       degauss=0.007,
       Emin   = -40.0, ! valor de menor energia encontrado nas bandas
       Emax   = 40.0,  !Valor de maior energia encontrado nas bandas
       DeltaE =0.05,
       filpdos = '$input$tipo.dat',
 /
EOF
                  set OMP_NUM_THREADS=1; export OMP_NUM_THREADS=1
                  /usr/bin/mpirun.mpich -np 16 /home/usr/q-e-qe-6.7.0/bin/projwfc.x -i  $input$tipo.in $
                  calc6=0
               else
                  calc6=1
              fi
             done
    fi

/home/usr/q-e-qe-6.7.0/bin/sumpdos.x  *\($atomo1\)* > atomo_$atomo1-orbital_total.dat
/home/usr/q-e-qe-6.7.0/bin/sumpdos.x  *\($atomo1\)*\(s\) > atomo_$atomo1-orbital_s.dat
/home/usr/q-e-qe-6.7.0/bin/sumpdos.x  *\($atomo1\)*\(p\) > atomo_$atomo1-orbital_p.dat

/home/usr/q-e-qe-6.7.0/bin/sumpdos.x  *\($atomo2\)* > atomo_$atomo2-orbital_total.dat
/home/usr/q-e-qe-6.7.0/bin/sumpdos.x  *\($atomo2\)*\(s\) > atomo_$atomo2-orbital_s.dat
/home/usr/q-e-qe-6.7.0/bin/sumpdos.x  *\($atomo2\)*\(p\) > atomo_$atomo2-orbital_p.dat


mkdir pdos
mv *\($atomo1\)* pdos
mv *\($atomo1\)*\(s\)
mv *\($atomo1\)*\(p\)
mv *\($atomo2\)* pdos
mv *\($atomo2\)*\(s\)
mv *\($atomo2\)*\(p\)


# CALCULO BULK MODULUS

rm -f bulk_$input.dat
cat > $input.in.tmpl <<EOF
&CONTROL
        title = '$input' ,
        calculation = 'scf' ,
        prefix='$input_@ALAT@',
        pseudo_dir='/home/iprodrigues/pseudopotenciais',
        outdir='temp_files2',
        verbosity='low',
        restart_mode='from_scratch',
        tprnfor=.true.,
        tstress=.true.,
        forc_conv_thr=1.0d-6
        etot_conv_thr=1.0d-8
        nstep=1000
 /

 &SYSTEM
        ibrav = 0,
        celldm(1) = @ALAT@,
        nat = 8,
        ntyp        = $tipo_atomos,
        ecutwfc     = $ewfc ,
        ecutrho     = $erho ,
        occupations= 'smearing',
        smearing='mv',
        degauss=0.005d0,
        lspinorb=.false.,
        noncolin=.false.,
 /

&ELECTRONS
        electron_maxstep = 1000,
        conv_thr=1.0d-12,
        mixing_beta      = $mixbeta,
        mixing_mode='plain',
/

ATOMIC_SPECIES
$pseudopotencial1
$pseudopotencial2

$geometria

$k_points

CELL_PARAMETERS {crystal}
$param_rede

EOF

for alat in `seq -w 09.00 0.05 12.00`
do
  sed "s/@ALAT@/$alat/" $input.in.tmpl > $input-$alat.in
  set OMP_NUM_THREADS=1; export OMP_NUM_THREADS=1
  /usr/bin/mpirun.mpich -np 16 /home/usr/q-e-qe-6.7.0/bin/pw.x < $input-$alat.in > $input-$alat.out
  ene=`grep ! $input-$alat.out | cut -b33-50`
  vol=`grep volume $input-$alat.out | cut -b33-46`
  echo $vol  $ene >> bulk_vol_$input.dat
  echo $alat $ene >> bulk_latt_$input.dat
  rm -r temp_files2
  mv  $input-$alat.in $input-$alat.out bulk/
#  rm -rf $input-$alat.in $input-$alat.out
done
