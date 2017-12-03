data_targets = new_data.dat 
image_targets = grafica.png 

Resulta_hw5.pdf : Results_hw5.tex $(image_targets) grafica.png
	pdflatex $< && rm *.aux *.log

grafica.png : Plots.py
	python $< 

$(image_targets) : Plots.py $(data_targets)
	python $<

$(data_targets) : a.out
	./a.out

a.out : CurvaRotacion.c
	gcc -lm -g CurvaRotacion.c
