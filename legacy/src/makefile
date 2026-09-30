# Makefile-Datei f�r MHD Riemann-Solver

CC = g++
CFLAGS = -Wall
 
rpx : MHD_Main1.o  MHD_Residuum.o MHD_OutPut.o MHD_Rarefaction.o MHD_fShock.o MHD_sShock.o MHD_fRareMax.o MHD_sShockMax.o DiffGlgn.o Newton.o tools.o
	$(CC) -o rpx MHD_Main1.o  MHD_Residuum.o MHD_OutPut.o MHD_Rarefaction.o MHD_fShock.o MHD_sShock.o MHD_fRareMax.o MHD_sShockMax.o DiffGlgn.o Newton.o tools.o
    
MHD_Main1.o : MHD_Main1.c MHD1.h
	$(CC) $(CFLAGS) -c MHD_Main1.c
     
MHD_Residuum.o : MHD_Residuum.c MHD1.h
	$(CC) $(CFLAGS) -c MHD_Residuum.c

MHD_OutPut.o : MHD_OutPut.c MHD1.h
	$(CC) $(CFLAGS) -c MHD_OutPut.c

MHD_Rarefaction.o : MHD_Rarefaction.c MHD1.h
	$(CC) $(CFLAGS) -c MHD_Rarefaction.c

MHD_fShock.o : MHD_fShock.c MHD1.h 
	$(CC) $(CFLAGS) -c MHD_fShock.c

MHD_sShock.o : MHD_sShock.c MHD1.h 
	$(CC) $(CFLAGS) -c MHD_sShock.c

MHD_fRareMax.o : MHD_fRareMax.c MHD1.h
	$(CC) $(CFLAGS) -c MHD_fRareMax.c

MHD_sShockMax.o : MHD_sShockMax.c MHD1.h
	$(CC) $(CFLAGS) -c MHD_sShockMax.c 

DiffGlgn.o :  DiffGlgn.c dopri5.h
	$(CC) $(CFLAGS) -c DiffGlgn.c	
 
Newton.o : Newton.c 
	$(CC) $(CFLAGS) -c Newton.c

tools.o : tools.c 
	$(CC) $(CFLAGS) -c tools.c

clean :
	rm -f *.o *~
