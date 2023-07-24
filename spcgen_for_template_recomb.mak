spcgen_for_template_recomb: spcgen_for_template_recomb.o nrutil.o cal_lst.o polcof.o polint.o gasdev.o \
	beam_definition.o precess.o spline.o splint.o refraction.o ran1.o fit.o gammq.o gcf.o gammln.o gser.o
	gfortran -o spcgen_for_template_recomb spcgen_for_template_recomb.o nrutil.o cal_lst.o precess.o gasdev.o \
	polcof.o polint.o beam_definition.o spline.o splint.o refraction.o ran1.o fit.o gammq.o gcf.o gammln.o gser.o \
	 -L/usr/local/miriad/darwin_x86_64/lib \
	 -L/usr/local/lib \
	 -L/usr/lib \
	 -lmir \
	 -lm
spcgen_for_template_recomb.o : spcgen_for_template_recomb.c
	gfortran -c spcgen_for_template_recomb.c
nrutil.o : nrutil.c
	gfortran -c nrutil.c
cal_lst.o : cal_lst.c
	  gfortran -c cal_lst.c
polcof.o : polcof.c
	gfortran -c polcof.c
polint.o : polint.c
	gfortran -c polint.c
beam_definition.o : beam_definition.c
	gfortran -c beam_definition.c
fit.o : fit.c
	gfortran -c fit.c
gammq.o : gammq.c
	gfortran -c gammq.c
gcf.o : gcf.c
	gfortran -c gcf.c
gammln.o : gammln.c
	gfortran -c gammln.c
gser.o : gser.c
	gfortran -c gser.c
precess.o : precess.f
	gfortran -c precess.f
spline.o : spline.c
	gfortran -c spline.c
splint.o : splint.c
	gfortran -c splint.c
gasdev.o : gasdev.c
	gfortran -c gasdev.c
ran1.o : ran1.c
	gfortran -c ran1.c
refraction.o : refraction.c
	gfortran -c refraction.c
