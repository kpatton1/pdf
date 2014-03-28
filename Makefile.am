bin_PROGRAMS = pdf

pdf_SOURCES = progs/main.cpp src/GF.cpp src/Halo.cpp src/LSS.cpp src/MF.cpp src/NFW.cpp src/PDF.cpp src/PS.cpp src/TF.cpp
pdf_CPPFLAGS = -I$(top_srcdir)/include $(GSLDEPS_CFLAGS)
pdf_LDADD = $(GSLDEPS_LIBS)