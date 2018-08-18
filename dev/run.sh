PATHWAY=$1
NATIVE=native-`basename $1`.pdb
NCORES=1

echo ">>>" run1.sh NATIVE PATHWAY "<<<"

pr00_main.py $PATHWAY shm/out 40 0.5 8 $NCORES

echo ""
echo ">>>>>>>>>>>>> PLOTTING <<<<<<<<<<<<<<"
echo ""
pathway-plotting-tmscore.R $NATIVE  shm/out/pdbs $NCORES
pathway-plotting-tmscore.R $NATIVE  shm/out/pdbsLocal $NCORES
pathway-plotting-tmscore.R $NATIVE  shm/out/pdbsGlobal $NCORES

mkdir -p shm/out/pdfs
mv pdbs.pdf pdbsLocal.pdf pdbsGlobal.pdf shm/out/pdfs

