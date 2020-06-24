for i in {0..9..1}
do
  ./gambit -rf yaml_files/XENON1T_evidences/bkg_tritium_${i}.yaml > bkg_tritium_${i}.log
  ./gambit -rf yaml_files/XENON1T_evidences/signal_${i}.yaml > signal_${i}.log
done
