#!/bin/bash
for p in "kt" "diff" "quartic"
do
	for i in {1..6}
	do
		filename="job_"$p"_"$i".pbs"
		echo "#!/bin/bash">$filename
		echo "#PBS -l nodes=1:ppn=1,mem=256mb,walltime=08:00:00">>$filename
		echo "#PBS -V">>$filename
		echo "cd /storage/molsim/phuhzr/lattice_switch/jobs">>$filename
		echo "julia run_"$p"_"$i".jl">>$filename
	done
done
