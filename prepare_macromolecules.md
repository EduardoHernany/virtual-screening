recebe caminho do receptor.pdb "--r"

recebe tamanho do gridbox "--gsize x, y, z" 
recebe centro do gridbox "--gcenter x, y, z"

se a configuracao for com redocking recebe o caminho do ligante "--l"

----

caminho pythonsh_path "~/mgltools_x86_64Linux2_1.5.7/bin/pythonsh"
caminho prepare_recptor "~/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
caminho prepare_ligand "~/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
caminho prepare_gpf "~/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py"
caminho ad4_parameters_path "~/x86_64Linux2/AD4_parameters.dat"
camonho autogrid4_path "~/x86_64Linux2/autogrid4"

----

recebeu o receptor. salva na macros_preparadas/recptor.filename/

pythonsh_path prepare_recptor -r recptor_path
#vai gerar um aqruivo como mesmo nome do recptor.pdbqt


pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_1.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=C,A,N,NA,NS,OA,OS,SA,S,H,HD
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_2.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=HS,P,Br,BR,Ca,CA,Cl,CL,F,Fe,FE
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_3.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=I,Mg,MG,Mn,MN,Zn,ZN,He,Li,Be
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_4.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=B,Ne,Na,Al,Si,K,Sc,Ti,V,Co
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_5.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=Ni,Cu,Ga,Ge,As,Se,Kr,Rb,Sr,Y
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_6.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=Zr,Nb,Cr,Tc,Ru,Rh,Pd,Ag,Cd,In
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_7.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=Sn,Sb,Te,Xe,Cs,Ba,La,Ce,Pr,Nd
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_8.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_9.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_10.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th
pythonsh_path prepare_gpf -r recptor.pdbqt -o grid_11.gpf -p gridcenter=x,y,z -p npts=x,y,z -p ligand_types=Pa,U,Np,Pu,Am,Cm,Bk,Cf,E,Fm

sed -i '1i\\parameter_file ad4_parameters_path' path_grid_1.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_2.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_3.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_4.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_5.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_6.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_7.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_8.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_9.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_10.gpf
sed -i '1i\\parameter_file ad4_parameters_path' path_grid_11.gpf

autogrid_path -p grid_1.gpf -l gridbox.glg
autogrid_path -p grid_2.gpf -l gridbox.glg
autogrid_path -p grid_3.gpf -l gridbox.glg
autogrid_path -p grid_4.gpf -l gridbox.glg
autogrid_path -p grid_5.gpf -l gridbox.glg
autogrid_path -p grid_6.gpf -l gridbox.glg
autogrid_path -p grid_7.gpf -l gridbox.glg
autogrid_path -p grid_8.gpf -l gridbox.glg
autogrid_path -p grid_9.gpf -l gridbox.glg
autogrid_path -p grid_10.gpf -l gridbox.glg
autogrid_path -p grid_11.gpf -l gridbox.glg

#nisso vai gerar um arquivo .maps.fld 
#temos que modificar o conteudo do arquivo
#apagar tudo abaixo da linha 23
#from util import textfld
#texto = textfld()
#novo_texto = texto.replace("kakakakaka", receptorfilename)


['/home/eduardo/mgltools_x86_64Linux2_1.5.7/bin/pythonsh '/home/eduardo/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py '-r 'macros_preparadas/1cjb_a/1cjb_a.pdbqt '-o 'macros_preparadas/1cjb_a/grid_11.gpf '-p 'gridcenter=15.5,20.3,10.2 '-p 'npts=60,60,60 '-p 'ligand_types=Pa,U,Np,Pu,Am,Cm,Bk,Cf,E,Fm']


/home/eduardo/mgltools_x86_64Linux2_1.5.7/bin/pythonsh /home/eduardo/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -r /home/eduardo/Documentos/Projetos/Fiocruz/docking/macros_preparadas/1cjb_a/1cjb_a.pdbqt -o /home/eduardo/Documentos/Projetos/Fiocruz/docking/macros_preparadas/1cjb_a/grid_11.gpf -p gridcenter=15.5,20.3,10.2 -p npts=60,60,60 -p ligand_types=Pa,U,Np,Pu,Am,Cm,Bk,Cf,E,Fm