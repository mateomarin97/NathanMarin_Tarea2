plotshock :
	python Shock.py
shock :
	./Shock.x
exec :
	gcc Shock.c -lm -o Shock.x
clean :
	rm  ushock.dat pshock.dat rhoshock.dat u.pdf p.pdf rho.pdf Shock.x tshock.dat
