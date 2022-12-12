
cimport cython

""" Low-level cython kernels for reach proximity testing
"""
# Authors: Darren Engwirda

@cython.boundscheck(False)  # deactivate bnds checking
def check_short(list rnet, rdat):
    """
    Check whether reach RDAT is "too short" wrt. the desired
    spatial proximity metric.

    """

    cdef double xxdn, yydn, zzdn, hhdn
    cdef double xxup, yyup, zzup, hhup
    cdef double xdel, ydel, zdel, dsqr, hbar, hsqr

    hhup = rdat.vert[+0].hval
    xxup = rdat.vert[+0].xpos
    yyup = rdat.vert[+0].ypos
    zzup = rdat.vert[+0].zpos

    hhdn = rdat.vert[-1].hval
    xxdn = rdat.vert[-1].xpos
    yydn = rdat.vert[-1].ypos
    zzdn = rdat.vert[-1].zpos

    hbar = 0.5 * (hhup + hhdn)

    xdel = xxup - xxdn
    ydel = yyup - yydn
    zdel = zzup - zzdn

    dsqr = xdel * xdel + \
           ydel * ydel + \
           zdel * zdel

    hsqr = hbar * hbar
    
    return dsqr >= (5 ** 2) * hsqr


@cython.boundscheck(False)  # deactivate bnds checking
def check_close(list rnet, rdat, list vert, list near):
    """
    Check whether reach RDAT is "too close" wrt. the desired
    spatial proximity metric.
    """

    cdef size_t keep = 1
    cdef double xrdn, yrdn, zrdn, hrdn
    cdef double xrpj, yrpj, zrpj, hrpj
    cdef double xkpj, ykpj, zkpj, hkpj
    cdef double xdel, ydel, zdel, dsqr, hbar, hsqr

    xrdn = rdat.vert[-1].xpos
    yrdn = rdat.vert[-1].ypos
    zrdn = rdat.vert[-1].zpos
    hrdn = rdat.vert[-1].hval

    for rvrt in vert:

        xrpj = rvrt.xpos; yrpj = rvrt.ypos; zrpj = rvrt.zpos
        hrpj = rvrt.hval

        if (rdat.dpos >= +0):

            xdel = xrpj - xrdn
            ydel = yrpj - yrdn
            zdel = zrpj - zrdn

            dsqr = xdel * xdel + \
                   ydel * ydel + \
                   zdel * zdel

            hsqr = hrdn * hrdn

            if (dsqr < 1. / 1. * hsqr): continue

        for kpos in near:
            kdat = rnet[kpos]
            if (kdat.rank < rdat.rank): continue

            for kvrt in reversed(kdat.vert):

                xkpj = kvrt.xpos; ykpj = kvrt.ypos; zkpj = kvrt.zpos
                hkpj = kvrt.hval

                xdel = xkpj - xrdn
                ydel = ykpj - yrdn
                zdel = zkpj - zrdn

                dsqr = xdel * xdel + \
                       ydel * ydel + \
                       zdel * zdel

                hsqr = hrdn * hrdn

                if (dsqr < 1. / 1. * hsqr): continue

                hbar = 0.5 * (hkpj + hrpj)

                xdel = xkpj - xrpj
                ydel = ykpj - yrpj
                zdel = zkpj - zrpj

                dsqr = xdel * xdel + \
                       ydel * ydel + \
                       zdel * zdel

                hsqr = hbar * hbar

                if (dsqr < 7. / 8. * hsqr): keep = 0; break
                
            if (keep == 0): break

        if (keep == 0): break

    return (keep == 1)


@cython.boundscheck(False)  # deactivate bnds checking
def check_seeds(list rnet, rdat, int base, list near):
    """
    Check whether "seed" vertices could be added to the mesh
    without violating the proximity metric.

    """

    cdef size_t okay = 1
    cdef double xbse, ybse, zbse, hbse
    cdef double xkdn, ykdn, zkdn, hkdn    
    cdef double xkup, ykup, zkup, hkup
    cdef double xdel, ydel, zdel, dsqr, hbar, hsqr

    hbse = rdat.vert[base].hval
    xbse = rdat.vert[base].xpos
    ybse = rdat.vert[base].ypos
    zbse = rdat.vert[base].zpos

    hsqr = hbse * hbse

    for kpos in near:
        kdat = rnet[kpos]
        if (kdat.rank >= rdat.rank):
            if (kdat.vert[+0].seed >= 1):

                xkup = kdat.vert[+0].xpos
                ykup = kdat.vert[+0].ypos
                zkup = kdat.vert[+0].zpos

                xdel = xkup - xbse
                ydel = ykup - ybse
                zdel = zkup - zbse

                dsqr = xdel * xdel + \
                       ydel * ydel + \
                       zdel * zdel

                if (dsqr < 7. / 8. * hsqr): okay = 0

            if (kdat.vert[-1].seed >= 1):

                xkdn = kdat.vert[-1].xpos
                ykdn = kdat.vert[-1].ypos
                zkdn = kdat.vert[-1].zpos

                xdel = xkdn - xbse
                ydel = ykdn - ybse
                zdel = zkdn - zbse

                dsqr = xdel * xdel + \
                       ydel * ydel + \
                       zdel * zdel

                if (dsqr < 7. / 8. * hsqr): okay = 0

        if (okay == 0): break

    return (okay == 1)