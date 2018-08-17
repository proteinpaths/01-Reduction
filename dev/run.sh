PATHWAY=$1
NATIVE=native-`basename $1`.pdb
echo ">>>" run1.sh NATIVE PATHWAY "<<<"

pr00_main.py $PATHWAY shm/out 100 0.5 5 4

#pathway-plotting-tmscore.R $NATIVE  /home/ppath/01-Reduction/dev/shm/out/pdbs 4
#pathway-plotting-tmscore.R $NATIVE  /home/ppath/01-Reduction/dev/shm/out/pdbsLocal 4
#pathway-plotting-tmscore.R $NATIVE  /home/ppath/01-Reduction/dev/shm/out/pdbsGlobal 4

#mkdir -p shm/out/pdfs
#mv pdbs.pdf pdbsLocal.pdf pdbsGlobal.pdf shm/out/pdfs

