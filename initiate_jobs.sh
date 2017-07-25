#!/bin/bash
for f in /storage/molsim/phuhzr/job_dir/*.pbs; do
	qsub $f
done

