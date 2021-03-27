CC = g++
OBJ = Orbit.o ObsData.o Observatory.o Vec_DP.o gaussj.o newt.o fdjac.o ludcmp.o epv00.o cal2jd.o zbrent.o mrqmin.o pvtob.o era00.o CSZ.o Planet.o kalman.o

fig1:	fig1.o $(OBJ)
	$(CC) fig1.o $(OBJ)
fig2:	fig2.o $(OBJ)
	$(CC) fig2.o $(OBJ)