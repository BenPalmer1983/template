#!/bin/bash
# Make eampa bash script
#-------------------------------------------------------------------------------
# Load arguments to array
args=("$@")
# Change to working directory (if not already in)
cd ${0%/*}
# Store processor count
procCount="$(nproc 2>&1)"
# Store working directory in variable
workingDir=$(pwd)
# Directories
srcDir=$workingDir"/src"
modDir=$workingDir"/mod"
binDir=$workingDir"/bin"
libDir=$workingDir"/libBP"
mkdir -p $modDir
mkdir -p $binDir
# Flags to refresh/build library
refreshLib=0
buildlib=0
for i in "${args[@]}"
do
	if [ $i == "refreshlib" ]; then
    refreshLib=1
    buildLib=1
  fi
  if [ $i == "buildlib" ]; then
    buildLib=1
  fi
done
# Update/refresh Library files if required (dev only)
libRootDir="/code/lib"
libDestDir="/code/template"
if [[ $refreshLib == 1 ]]; then
  echo "Refreshing library files."
  rm -fR $libDestDir"/libBP"
  rm -f $libDestDir"/makeLib.sh"
  mkdir -p $libDestDir"/libBP"
  cp -R $libRootDir"/libBP/src" $libDestDir"/libBP"
  cp $libRootDir"/makeLib.sh" $libDestDir
fi
# Build library
if [[ $buildLib == 1 ]]; then
  echo "Building library files."
  "$workingDir/makeLib.sh"
fi
# -----------Library Start--------
modDirLib=$libDir"/mod"
libDirLib=$libDir"/lib"
libLink="-L"$libDirLib" -lBP"
# -----------Library End--------
binFile="template.x "
fortLine="mpif90 -g \
-fbounds-check -mtune=native -O3 " # -Wno-unused-function -fcheck=all  -Wall -O3
# Build files
buildFiles="types.f90 \
msubs.f90 \
loadData.f90 \
globals.f90 \
output.f90 \
initialise.f90 \
input.f90 \
finalise.f90 \
template.f90"
# Update compile time in globals.f90
#python make.py 0
# Compile
cd $srcDir
commandLine=$fortLine" -J "$modDir"  -I"$modDirLib" "
commandLine=$commandLine" "$buildFiles" "
commandLine=$commandLine" "$libLink" "
commandLine=$commandLine" -o "$binDir"/"$binFile
echo "----------------------------------------------------------------------------------"
echo "Compile with Fortran:"
echo " "
echo $commandLine
echo "----------------------------------------------------------------------------------"
eval $commandLine
# add export to profile so user can run activity.x
exportLine="export PATH=\"\$PATH:"$binDir"\""
profileFile=$HOME"/.bash_profile"
touch $profileFile
if grep -q "$binDir" "$profileFile";
then
  echo $profileFile" already has the bin directory path"
else
  echo "Adding "$binDir" to "$profileFile
  echo $exportLine >> $profileFile
fi
profileFile=$HOME"/.bashrc"
touch $profileFile
if grep -q "$binDir" "$profileFile";
then
  echo $profileFile" already has the bin directory path"
else
  echo "Adding "$binDir" to "$profileFile
  echo $exportLine >> $profileFile
fi
