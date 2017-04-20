u.pdf p.pdf rho.pdf : Shock.py ushock.dat pshock.dat rhoshock.dat
	python Shock.py
ushock.dat pshock.dat rhoshock.dat : Shock.x
	./Shock.x
Shock.x : Shock.c
	gcc -lm Shock.c -o Shock.x
clean :
	rm  ushock.dat pshock.dat rhoshock.dat u.pdf p.pdf rho.pdf Shock.x
