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

cleankaboom :
	rm kaboom.x p1sedov.dat p2sedov.dat p3sedov.dat rho1sedov.dat rho2sedov.dat rho3sedov.dat rho10sedov.pdf rho60sedov.pdf rho120sedov.pdf t1sedov.dat t2sedov.dat t3sedov.dat u1sedov.dat u2sedov.dat u3sedov.dat
