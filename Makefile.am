bin_PROGRAMS = pdf nongaussian test test_cosmo test_mf test_ps_lin test_ps_nonlin test_ps_var integrate integrate_logx integrate_logrithmic integrate_weighted integrate_std

AM_CPPFLAGS = -I$(top_srcdir)/include $(GSLDEPS_CFLAGS)
LDADD = $(GSLDEPS_LIBS)

pdf_SOURCES = progs/main.cpp src/GF.cpp src/Halo.cpp src/LSS.cpp src/MF.cpp src/NFW.cpp src/PDF.cpp src/PS.cpp src/TF.cpp

nongaussian_SOURCES = progs/nongaussian_main.cpp src/CDF.cpp

test_SOURCES = progs/test.cpp

test_mf_SOURCES = progs/test_mf.cpp src/MF.cpp src/GF.cpp src/PS.cpp src/TF.cpp

test_ps_lin_SOURCES = progs/test_ps_lin.cpp src/GF.cpp src/PS.cpp src/TF.cpp
test_ps_nonlin_SOURCES = progs/test_ps_nonlin.cpp src/GF.cpp src/PS.cpp src/TF.cpp

test_ps_var_SOURCES = progs/test_ps_var.cpp src/GF.cpp src/PS.cpp src/TF.cpp

integrate_SOURCES = progs/integrate.cpp
integrate_logx_SOURCES = progs/integrate_logx.cpp
integrate_logrithmic_SOURCES = progs/integrate_logrithmic.cpp
integrate_weighted_SOURCES = progs/integrate_weighted.cpp
integrate_std_SOURCES = progs/integrate_std.cpp

test_cosmo_SOURCES = progs/test_cosmo.cpp src/cosmo.cpp