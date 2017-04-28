plotshock :
	python Shock.py
shock :
	./Shock.x
sedov :
	./Sedov.x
plotsedov :
	python Sedov.py
exec :
	gcc Shock.c -lm -o Shock.x
	gcc Sedov.c -lm -o Sedov.x
clean :
	rm  ushock.dat pshock.dat rhoshock.dat u.pdf p.pdf rho.pdf Shock.x tshock.dat
