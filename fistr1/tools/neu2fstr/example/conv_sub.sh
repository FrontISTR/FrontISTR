#!/bin/sh


neu2fstr=../../neu2fstr

for i in ${model}
do
	neu="${i}.NEU"
	msh="${i}.msh"
	cnt="${i}.cnt"
	echo ${neu2fstr} ${1} ${analysis} ${neu} ${msh} ${cnt}
	${neu2fstr} ${1} ${analysis} ${neu} ${msh} ${cnt}
done 


for i in ${model_shell}
do
	neu="${i}.NEU"
	msh="${i}.msh"
	cnt="${i}.cnt"
	echo ${neu2fstr} -d ${analysis} ${neu} ${msh} ${cnt}
	${neu2fstr} -d ${analysis} ${neu} ${msh} ${cnt}
done 

