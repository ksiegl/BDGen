OBJECTS = BD2.o BDA021Z11.o CompUVPASpline.o KillUVPSpline.o \
	MakeUVPSpline.o UVPSint.o UVPSloc.o UVPSval.o c3ran.o ran2.o \
	GaussianDev.o c3ranradius.o
bd : $(OBJECTS)
	g++ -o bd $(OBJECTS)
BD2.o : BD2.c BD2.h BDA021Z11.h ran2.h UVPSpline.h CompUVPASpline.h \
	UVPSval.h
	g++ -c BD2.c
BDA021Z11.o : BDA021Z11.c BD2.h BDA021Z11.h c3ran.h \
	GaussianDev.h c3ranradius.h
	g++ -c BDA021Z11.c
CompUVPASpline.o : CompUVPASpline.c UVPSpline.h CompUVPASpline.h \
	MakeUVPSpline.h KillUVPSpline.h
	g++ -c CompUVPASpline.c
KillUVPSpline.o : KillUVPSpline.c KillUVPSpline.h UVPSpline.h
	g++ -c KillUVPSpline.c
MakeUVPSpline.o : MakeUVPSpline.c UVPSpline.h KillUVPSpline.h
	g++ -c MakeUVPSpline.c
UVPSint.o : UVPSint.c UVPSint.h UVPSpline.h UVPSloc.h
	g++ -c UVPSint.c
UVPSloc.o : UVPSloc.c UVPSloc.h UVPSpline.h
	g++ -c UVPSloc.c
UVPSval.o : UVPSval.c UVPSval.h UVPSpline.h UVPSloc.h
	g++ -c UVPSval.c
c3ran.o : c3ran.c c3ran.h BD2.h
	g++ -c c3ran.c
ran2.o : ran2.c ran2.h
	g++ -c ran2.c
GaussianDev.o : GaussianDev.h
	g++ -c GaussianDev.c
c3ranradius.o : c3ranradius.c c3ranradius.h BD2.h
	g++ -c c3ranradius.c
clean :
	rm $(OBJECTS)
